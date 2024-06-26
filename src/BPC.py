import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BPC:

    #---instantiate a BPC algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.INT_tol = 1e-4
        self.OA_tol = 1e-4
        self.CG_tol = 1e-4
        self.inf = 1e+9
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_RMP = 0.0
        self.rt_pricing = 0.0        
        
        self.OAcuts = []
        self.yvec = [] #---for global interdiction cuts
        
        n = BB_node.BB_node(self.network, 0, 0, self.LB, self.inf, [], [], False)
        self.BB_nodes.append(n)
        
    def getCandidates(self):
        return [n for n in self.BB_nodes if n.active==True]
        
    def getLB(self, candidates):
        return min([i.LB for i in candidates])
        
    def getGap(self):           
        return (self.UB - self.LB)/self.UB
    
    def nodeSelection(self, candidates):
        min_LB = min([n.LB for n in candidates])
        min_LB_nodes = [n for n in candidates if n.LB==min_LB]
        return min_LB_nodes[0]
    
    def branch(self, can):

        fixed00 = list(can.fixed0)
        fixed00.append(can.ybr)
        fixed01 = list(can.fixed1)
        fixed10 = list(can.fixed0)
        fixed11 = list(can.fixed1)
        fixed11.append(can.ybr)        
        
        cnt = len(self.BB_nodes) 
        
        BB_node_id = cnt
        can.children.append(BB_node_id)
        n0 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, self.inf, fixed00, fixed01, False)
        n0.paths = dict(can.paths)
        self.BB_nodes.append(n0)
        
        BB_node_id = cnt+1
        can.children.append(BB_node_id)
        n1 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, self.inf, fixed10, fixed11, False)
        n1.paths = dict(can.paths)
        self.BB_nodes.append(n1)
    
        return
    
    def getOAcut(self):
        cut = {'a':{}, 'b':{}}
        
        for a in self.network.links:
        
            if a.y == 1:                
                cut['a'][a] = round(a.getTravelTime(a.x,'SO'), 3)
                cut['b'][a] = round(- pow(a.x, 2) * a.getDerivativeTravelTime(a.x), 3)
                          
            else:
                cut['a'][a] = 0.0
                cut['b'][a] = 0.0
                
        return cut
    
    def knapsack(self, type, can):
        
        knp = Model()
        knp.y = {a:knp.binary_var() for a in self.network.links2}
        knp.add_constraint(sum(knp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        #---Interdiction cuts
        for yv in self.yvec:
            knp.add_constraint(sum(knp.y[a] + yv[a] - 2*knp.y[a]*yv[a] for a in self.network.links2) >= 1)
            
        #---Branch cuts
        for a in self.network.links2:
            
            if a.id in can.fixed0:
                knp.add_constraint(knp.y[a] == 0)
                               
            if a.id in can.fixed1:
                knp.add_constraint(knp.y[a] == 1)
        
        if type == 'capacity':
            knp.maximize(sum(knp.y[a] * a.C for a in self.network.links2))
        elif type == 'x':
            knp.maximize(sum(knp.y[a] * a.x for a in self.network.links2))                        
        else:
            print('unknown type')
            
        knp.solve(log_output=False)
            
        if knp.solve_details.status == 'infeasible' or knp.solve_details.status == 'integer infeasible':
            print('infeasible instance?')
            return {}
        else:            
            yKNP = {}
            for a in self.network.links2:
                yKNP[a] = knp.y[a].solution_value
                
            return yKNP
        
    def initOAcuts(self, can, nKNP):
        
        if self.params.PRINT_BB_INFO:
            print('Running knapsack heuristic with',nKNP,'round(s)')
        
        yKNP = {}
        for a in self.network.links2:
            yKNP[a] = 1
        
        t0_TAP = time.time()
        tstt = round(self.network.tapas('SO',yKNP), 3) 
        self.ydict.insertSO(yKNP, tstt)
        self.rt_TAP += (time.time() - t0_TAP)
        
        best = self.inf
        for n in range(nKNP):
            #yKNP = self.knapsack('capacity',can)
            yKNP = self.knapsack('x',can)
                        
            t0_TAP = time.time()
            tstt = round(self.network.tapas('SO',yKNP), 3) 
            self.ydict.insertSO(yKNP, tstt)
            self.rt_TAP += (time.time() - t0_TAP)
            self.OAcuts.append(self.getOAcut())                            
            self.nSO += 1        
                    
            self.yvec.append(yKNP)
            
            if tstt < best:
                best = tstt
                self.yopt = yKNP
                
        t0_TAP = time.time()
        can.UB = round(self.network.tapas('UE',self.yopt), 3)
        self.ydict.insertUE(self.yopt, can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB
    
    def checkIntegral(self, y):
        
        #---check integrality
        frac = []
        for a in self.network.links2:
            if y[a] > self.INT_tol and y[a] < 1 - self.INT_tol:
                frac.append(a)
        
        if len(frac) == 0:
            return 'integral',frac
        else:
            return 'fractional',frac
    
    def initPaths(self, can):
        for a in self.network.links2:
            a.y = 0
            
        new = 0
        for r in self.network.origins:        
            self.network.dijkstras(r,'UE')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                
                    p = self.network.trace(r,s)
                    can.paths[r][s].append(p)                
                    #print(r.id, s.id, p.links)                
                    new += 1
                
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('#paths',new)
        
    
    def pricing(self, can):
        
        t0_pricing = time.time()
        
        new = 0
        minrc = self.inf
        
        for a in self.network.links:
            a.dual = can.duals['link'][a]
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                
                    rc = - can.duals['dem'][(r,s)] + s.cost                    
                    
                    if rc <= - self.CG_tol:
                        p = self.network.trace(r,s)
                        can.paths[r][s].append(p)                 
                        new += 1
                    
                    if rc < minrc:
                        minrc = rc
            
        if self.params.PRINT_BB_INFO:
            print('pricing',new,minrc)
            
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
    
    
    def rmp_path(self, can, type):
        
        #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively
           
        t0_RMP = time.time()
        rmp = Model()
        
        rmp.x = {a:rmp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        rmp.h = {p:rmp.continuous_var(lb=0) for p in can.getPaths()}
        rmp.mu = rmp.continuous_var(lb=can.LB)
        
        if type == 'LP':
            rmp.y = {a:rmp.continuous_var(lb=0,ub=1) for a in self.network.links2}
        elif type == 'MILP':
            rmp.y = {a:rmp.binary_var() for a in self.network.links2}
                    
        rmp.add_constraint(sum(rmp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        for a in self.network.links2:
            rmp.add_constraint(rmp.x[a] <= rmp.y[a] * self.network.TD)
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    #rmp.add_constraint(sum(rmp.h[p] for p in can.paths[r][s]) == r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
                    rmp.add_constraint(sum(rmp.h[p] for p in can.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
    
        for a in self.network.links:
            #rmp.add_constraint(sum(rmp.h[p] for p in can.getPaths() if a in p.links) == rmp.x[a], 'link_%d_%d' % (a.start.id,a.end.id))
            rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in can.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))
                        
        for OAcut in self.OAcuts:
            #---OA cuts
            rmp.add_constraint(rmp.mu >= sum(rmp.x[a]*OAcut['a'][a] + OAcut['b'][a] for a in self.network.links))
            
        for yv in self.yvec:        
            #---Interdiction Cuts - useless unless MILP?
            rmp.add_constraint(sum(rmp.y[a] + yv[a] - 2*rmp.y[a]*yv[a] for a in self.network.links2) >= 1)
            
        #---Branch cuts        
        for a in self.network.links2:
            
            if a.id in can.fixed0:
                rmp.add_constraint(rmp.y[a] == 0)
                               
            if a.id in can.fixed1:
                rmp.add_constraint(rmp.y[a] == 1)
        
        rmp.minimize(rmp.mu)
        
        rmp.parameters.threads = self.params.CPLEX_threads
        
        rmp.solve(log_output=False)
        
        if self.params.PRINT_BB_INFO:
            print('nb of paths: %d, cplex time: %.1f, rmp time: %.1f' % (len(can.getPaths()),rmp.solve_details.time,(time.time() - t0_RMP)))
        
        self.rt_RMP += (time.time() - t0_RMP)
            
        if rmp.solve_details.status == 'infeasible' or rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            OFV = rmp.objective_value
            RMP_status = rmp.solve_details.status
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = rmp.y[a].solution_value             
            
            if type == 'LP':
                dual_link = {} 
                for a in self.network.links:
                    dual_link[a] = max(rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0)
                    #if dual_link[a] < 0: 
                    #    print('dual link:',a.id,dual_link[a])
            
                dual_dem = {}
                for r in self.network.origins:
                    for s in self.network.zones:
                        if r.getDemand(s) > 0:
                            dual_dem[(r,s)] = max(rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0)
                            #if dual_dem[(r,s)] < 0: 
                            #    print('dual dem:',r.id,s.id,dual_dem[(r,s)])
                    
                can.duals = {'link':dual_link,'dem':dual_dem}
                
        return RMP_status,OFV,yopt

    
    def CG(self, can):
        
        nCG = 0
        conv = False

        while conv == False:        

            RMP_status,OFV,yRMP = self.rmp_path(can,'LP')
    
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                break
    
            minrc = self.pricing(can)
            
            if self.params.PRINT_BB_INFO:
                npaths = len(can.getPaths())
                print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))
            
            if minrc >= -self.CG_tol:
                if self.params.PRINT_BB_INFO:
                    print('CG converged')
                conv = True
            
            nCG += 1
            
        if conv == True:            
            
            CG_status, can.frac = self.checkIntegral(yRMP)
            
            if self.params.PRINT_BB_INFO:
                print('CG status',CG_status)
            
            return CG_status,OFV,yRMP
        
        else:
            return CG_status,self.inf,yRMP
              

    def OA_lp_path(self,can):
        nOA = 0
        UB_OA = self.inf
        LB_OA = 0.0 #---value to be returned that will serve as the LB of the BB node
        yOA = {}
        
        t0_OA = time.time()
        
        conv = False
        while conv == False:            
                        
            CG_status,CG_OFV,yCG = self.CG(can)
            
            can.LB = CG_OFV
        
            if CG_status != 'integral':
                
                conv = True
                
                if CG_status == 'infeasible':
                    if nOA == 0:
                        print('==========WARNING: LP infeasible at first iteration===================')
                        #---lp must be feasible at first iteration since interdiction cuts are reset at each BB node
                        #   and budget constraint violations are checked before executing OA_link
                    
                    #---search is complete and UB_OA is the optimal OFV                    
                    LB_OA = max(can.LB,UB_OA)
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> convergence by feasibility (due to interdiction cuts)')
                    return nOA,LB_OA,yOA,CG_status
                
                else:
                    #---solution is fractional: stop OA an return LP_OFV as BB node LB                   
                    LB_OA = CG_OFV                    
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> fractional solution')                      
                    
                    if self.params.KNP:
                        if self.params.PRINT_BB_INFO:
                            print('-----------------------------------> solving KNP')                      
                        
                        yKNP = self.knapsack('x', can)
                        
                        if len(yKNP) > 0:
                            #---knp is feasible
                            t0_TAP = time.time()
                            tstt = round(self.network.tapas('SO',yKNP), 3)                
                            self.ydict.insertSO(yKNP, tstt)
                            self.rt_TAP += time.time() - t0_TAP
                            self.OAcuts.append(self.getOAcut())
                            self.nSO += 1
                    
                
            else:
                #---integral solution, generate OA cut and keep going
                        
                if CG_OFV >= self.UB:
                    #---search can be stopped and lp obj can be used as the LB of the BB node
                    conv = True
                    LB_OA = max(can.LB,CG_OFV)
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> convergence by bounding in OA_lp_path')
                    return nOA,LB_OA,yOA,CG_status
            
                #---add interdiction cut: in a single search tree yvec can be global
                self.yvec.append(yCG)
                
                #---solve SO-TAP to get OA cut
                if self.ydict.hasSO(yCG) == True:
                    tstt = self.ydict.getSO(yCG)
                    print('\n\nWARNING: with global interdiction cuts, yCG should not have been feasible!!!\n\n')
                
                else:
                    t0_TAP = time.time()
                    tstt = round(self.network.tapas('SO',yCG), 3)                
                    self.ydict.insertSO(yCG, tstt)
                    self.rt_TAP += time.time() - t0_TAP
                    self.OAcuts.append(self.getOAcut())
                    self.nSO += 1
            
                if tstt < UB_OA:
                    UB_OA = tstt
                    yOA = yCG
                
            for a in self.network.links2:
                can.score[a.id] = round(a.x * a.getTravelTime(a.x,'SO'), 3)
            
            gap = (UB_OA-CG_OFV)/UB_OA
            if self.params.PRINT_BB_INFO:
                print('-----------------------------------> %d\t%.1f\t%.1f\t%.2f%%' % (nOA,CG_OFV,UB_OA,100*gap))
            
            if gap <= self.OA_tol:
                #---search can be stopped and the min between lp obj and UB_OA can be used as the LB of the BB node
                conv = True
                LB_OA = max(can.LB,min(CG_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by optimality gap in OA_lp_path')                    
                
            if (time.time() - t0_OA) >= self.params.BB_timelimit/2:
                #---search is stopped and the min between lp obj and UB_OA can be used as the LB of the BB node
                LB_OA = max(can.LB,min(CG_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> time limit exceeded in OA_lp_path')
                break
            
            nOA += 1
            
        return nOA,LB_OA,yOA,CG_status

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---BPC---')        
        
        self.network.resetTapas()
 
        t0 = time.time()
        
        #---initialize OA cuts
        nInitCuts = len(self.network.links2)
        self.initOAcuts(self.BB_nodes[0],nInitCuts)        
    
        #---initialize paths
        self.initPaths(self.BB_nodes[0])
    
        conv = False
        while conv == False:                    
             
            can = self.nodeSelection(self.getCandidates())
            status = can.check()
        
            if self.params.PRINT_BB_INFO:
                print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
     
            prune = False
            integral = False
                     
            if status == 'infeasible':
                if self.params.PRINT_BB_INFO:
                    print('--> prune by feasibility')
                prune = True                 
                can.LB = self.inf
                can.UB = self.inf        
                 
            elif status == 'fixed' or status == 'stop':
                if self.params.PRINT_BB_INFO:
                    print('--> prune by check',status)
                prune = True
                integral = True
                 
                yUB = {}
                for a in self.network.links2:
                    if a.id in can.fixed1:
                        yUB[a] = 1
                    elif a.id in can.fixed0:
                        yUB[a] = 0 
                    else:
                        yUB[a] = 0
                        can.fixed0.append(a.id)        
                 
            else:
                #---LB is obtained from LP relaxation of OA master problem
                nOA,can.LB,yOA,CG_status = self.OA_lp_path(can)
                
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True 
                
                elif CG_status == 'integral':
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by integrality')
                        print('yOA',yOA)
                    prune = True
                    integral = True
                    yUB = yOA
                 
            if integral == True:
                
                #---solve UE TAP to get UB
                if self.ydict.hasUE(yUB) == True:
                    can.UB = self.ydict.getUE(yUB)
                
                else:
                    t0_TAP = time.time()
                    can.UB = round(self.network.tapas('UE',yUB), 3)
                    self.ydict.insertUE(yUB, can.UB)
                    self.rt_TAP += time.time() - t0_TAP
                    self.nUE += 1
                
                if can.UB < self.UB:            
                    self.UB = can.UB
                    self.yopt = yUB
                    if self.params.PRINT_BB_INFO:
                        print('--> update UB')
                    
                    for n in self.BB_nodes:                    
                        if n.active == True and n.LB >= self.UB:
                            n.active = False

            if prune == False:
                
                frac = [a.id for a in can.frac]                
                free_sorted = sorted(frac, key = lambda ele: can.score[ele], reverse = True)
                can.ybr = free_sorted[0]
                
                if self.params.PRINT_BB_INFO:                    
                    for a in self.network.links2:
                        if a.id == can.ybr:
                            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can.ybr))
                
                self.branch(can)
         
            can.active = False
            candidates = self.getCandidates()
               
            if len(candidates) == 0:
                conv = True
                self.LB = self.UB
                self.gap = 0.0
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by inspection')
                break
                
            else:
                self.LB = self.getLB(candidates)
                self.gap = self.getGap()
            
            if self.params.PRINT_BB_INFO:
                print('--> can (after): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))            
            
            if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
                print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (self.nBB,self.nSO,self.nUE,self.LB,self.UB,100*self.gap))
            
            if self.gap <= self.params.BB_tol:
                conv = True
                self.LB = self.UB
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by optimality gap')
                break
            
            if (time.time() - t0) >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print('--> time limit exceeded')
                break
            
            self.nBB += 1
 
        self.rt = time.time() - t0

        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,self.rt,self.nBB,self.nSO,self.nUE,self.UB,100*self.gap))
            print(self.rt_TAP)
            print(self.rt_RMP)
            print(self.rt_pricing)
            print(self.yopt)
        