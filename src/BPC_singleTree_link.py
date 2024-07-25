import time
import math
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BPC_singleTree_link:

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
        self.UB_SO_DNDP = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.t0 = 0.0
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_RMP = 0.0
        self.rt_pricing = 0.0        
        
        self.OAcuts = []
        self.yvec = [] #---for global interdiction cuts
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        
        n = BB_node.BB_node(self.network, 0, 0, self.LB, self.inf, [], [], False)
        self.BB_nodes.append(n)
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.network.zones:
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths        
        
    def getCandidates(self):
        return [n for n in self.BB_nodes if n.active==True]
        
    def getLB(self, candidates):
        return min([i.LB for i in candidates])
        
    def getGap(self):           
        return (self.UB - self.LB)/self.UB
    
    def nodeSelection(self, type, candidates):
        
        if type == 'bestBound':
            min_LB = min([n.LB for n in candidates])
            min_LB_nodes = [n for n in candidates if n.LB==min_LB]
            return min_LB_nodes[0]
            
        elif type == 'depth':
            max_depth = max([len(n.fixed0) + len(n.fixed1) for n in candidates])
            max_depth_nodes = [n for n in candidates if (len(n.fixed0) + len(n.fixed1) == max_depth)]
            return max_depth_nodes[0]
        
    def branch_unfixed(self, can):
                                     
        fixed00 = list(can.fixed0)
        fixed00.append(can.ybr)
        fixed01 = list(can.fixed1)
        fixed10 = list(can.fixed0)
        fixed11 = list(can.fixed1)
        fixed11.append(can.ybr)
        
        cnt = len(self.BB_nodes) 
        
        for a in self.network.links2:
            if can.ybr == a.id:
                if self.params.PRINT_BB_INFO:
                    print('a.id found',a.id,' y =',can.y[a])                
                break
            
        if can.y[a] == 0:
            solved0 = True
            solved1 = False
            UB0 = can.UB
            UB1 = self.inf
        else:
            solved0 = False
            solved1 = True
            UB0 = self.inf
            UB1 = can.UB          
        
        BB_node_id = cnt
        can.children.append(BB_node_id)
        n0 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, UB0, fixed00, fixed01, solved0)
        self.BB_nodes.append(n0)
        
        BB_node_id = cnt+1
        can.children.append(BB_node_id)
        n1 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, UB1, fixed10, fixed11, solved1)
        self.BB_nodes.append(n1)

        if solved0 == True:
            n0.y = can.y
            n0.score = can.score
        else:
            n1.y = can.y
            n1.score = can.score
    
        return    
    
    def branch_frac(self, can):

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
        self.BB_nodes.append(n0)
        
        BB_node_id = cnt+1
        can.children.append(BB_node_id)
        n1 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, self.inf, fixed10, fixed11, False)
        self.BB_nodes.append(n1)
    
        return      
    
    def getOAcut(self):
        cut = {'a':{}, 'b':{}}
        
        for a in self.network.links:
        
            if a.y == 1:                
                cut['a'][a] = a.getTravelTime(a.x,'SO')
                cut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
                          
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
                yKNP[a] = round(knp.y[a].solution_value)
                
            return yKNP
        
    def initOAcuts(self, can, nKNP):
        
        if self.params.PRINT_BB_BASIC:
            print('Running knapsack heuristic with',nKNP,'round(s)')
        
        yKNP = {}
        for a in self.network.links2:
            yKNP[a] = 1
        
        t0_TAP = time.time()
        tstt = self.network.tapas('SO',yKNP)
        self.ydict.insertSO(yKNP, tstt)
        self.OAcuts.append(self.getOAcut()) 
        self.rt_TAP += (time.time() - t0_TAP)
        self.nSO += 1
        
        self.xvec = self.network.getFlowMap() 
        
        for a in self.network.links2:
            can.score[a.id] = a.x * a.getTravelTime(a.x,'SO') 
        
        for n in range(nKNP):
            yKNP = self.knapsack('x',can)
                       
            if not self.params.useAONcuts: 
                t0_TAP = time.time()
                tstt = self.network.tapas('SO',yKNP)
                self.ydict.insertSO(yKNP, tstt)
                self.rt_TAP += (time.time() - t0_TAP)
                self.OAcuts.append(self.getOAcut())                            
                self.nSO += 1        
            
            if self.params.useAONcuts:
                t0_TAP = time.time()
                tstt = self.network.setAON('SO', yKNP)
                self.rt_TAP += (time.time() - t0_TAP)
                self.OAcuts.append(self.getOAcut())   
             
            if self.params.useInterdictionCuts:       
                self.yvec.append(yKNP)
            
            if tstt < self.UB_SO_DNDP:
                self.UB_SO_DNDP = tstt
                self.yopt = yKNP
                
                
                
        t0_TAP = time.time()
        can.UB = self.network.tapas('UE',self.yopt)
        self.ydict.insertUE(self.yopt, can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB
        
        '''
        allY = {}
        for a in self.network.links2:
            allY[a] = 1
            
        sotstt = self.network.tapas('SO',allY) 
        self.ydict.insertSO(allY, sotstt)
        self.OAcuts.append(self.getOAcut())
        self.nSO += 1
        '''
 
    
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
    
    def initPaths(self):
        for a in self.network.links2:
            a.y = 0
            
        new = 0
        for r in self.network.origins:        
            self.network.dijkstras(r,'UE')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                
                    p = self.network.trace(r,s)
                    self.paths[r][s].append(p)
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
                        self.paths[r][s].append(p)
                        new += 1
                    
                    if rc < minrc:
                        minrc = rc
            
        #if self.params.PRINT_BB_INFO:
        #    print('pricing',new,minrc)
            
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
    
    def rmp_path(self, can, type):
        
        #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively
           
        t0_RMP = time.time()
        rmp = Model()
        
        rmp.x = {a:rmp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        rmp.h = {p:rmp.continuous_var(lb=0) for p in self.getPaths()}
        rmp.mu = {a:rmp.continuous_var(lb=0) for a in self.network.links}
        
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
                    rmp.add_constraint(sum(rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))                    
    
        for a in self.network.links:
            rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))
                        
        for OAcut in self.OAcuts:
            #---OA cuts
            for a in self.network.links:
                rmp.add_constraint(rmp.mu[a] >= rmp.x[a]*OAcut['a'][a] + OAcut['b'][a])
        
        if self.params.useInterdictionCuts:
            for yv in self.yvec:        
                #---Interdiction Cuts - useless unless MILP?
                rmp.add_constraint(sum(rmp.y[a] + yv[a] - 2*rmp.y[a]*yv[a] for a in self.network.links2) >= 1)
        
        #---Branch cuts        
        for a in self.network.links2:
            
            if a.id in can.fixed0:
                rmp.add_constraint(rmp.y[a] == 0)
                               
            if a.id in can.fixed1:
                rmp.add_constraint(rmp.y[a] == 1)
        
        rmp.minimize(sum(rmp.mu[a] for a in self.network.links))
        
        rmp.parameters.threads = self.params.CPLEX_threads
        rmp.parameters.timelimit = self.params.BB_timelimit
        
        rmp.solve(log_output=False)
        
        #if self.params.PRINT_BB_INFO:
        #    print('nb of paths: %d, cplex time: %.1f, rmp time: %.1f' % (len(can.getPaths()),rmp.solve_details.time,(time.time() - t0_RMP)))
        
        self.rt_RMP += (time.time() - t0_RMP)
            
        if rmp.solve_details.status == 'infeasible' or rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            OFV = rmp.objective_value
            #print("RMP obj", OFV)
            RMP_status = rmp.solve_details.status
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = rmp.y[a].solution_value
                maxscore = 0
                for OAcut in self.OAcuts:
                    maxscore = max(rmp.x[a].solution_value * OAcut['a'][a] + OAcut['b'][a],maxscore)                    
                can.score[a.id] = rmp.x[a].solution_value * maxscore                

            
            if type == 'LP':
                dual_link = {} 
                for a in self.network.links:                    
                    dual_link[a] = max(rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0)                    
            
                dual_dem = {}
                for r in self.network.origins:
                    for s in self.network.zones:
                        if r.getDemand(s) > 0:
                            dual_dem[(r,s)] = max(rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0)
                    
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
            
            #if self.params.PRINT_BB_INFO:
            #    npaths = len(self.getPaths())
            #    print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))
            
            if minrc >= -self.CG_tol:
                #if self.params.PRINT_BB_INFO:
                #    print('CG converged')
                conv = True
            
            nCG += 1
            
        if conv == True:            
            
            CG_status, can.frac = self.checkIntegral(yRMP)
            
            if self.params.PRINT_BB_INFO:
                if CG_status != 'fractional':
                    print('CG status',CG_status)
                    
            #npaths = len(self.getPaths())
            #print('CG: %s\t%d\t%d\t%.1f\t%.2f' % (CG_status,nCG,npaths,OFV,minrc))        
            
            return CG_status,OFV,yRMP
        
        else:
            return CG_status,self.inf,yRMP
  

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---')
        
        self.network.resetTapas()
 
        self.t0 = time.time()
        
        #---initialize OA cuts
        nInitCuts = round(len(self.network.links2))
        self.initOAcuts(self.BB_nodes[0],nInitCuts)        
    
        #---initialize paths
        self.initPaths()        
    
        conv = False
        while conv == False:                    
             
            can = self.nodeSelection('bestBound',self.getCandidates())
            #can = self.nodeSelection('depth',self.getCandidates())
            status = can.check()
        
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
     
            prune = False
            runSO = False
            runUE = False
                     
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
                #runSO = True #---waste of time if too many OA cuts?
                runUE = True
                 
                for a in self.network.links2:
                    if a.id in can.fixed1:
                        can.y[a] = 1
                    elif a.id in can.fixed0:
                        can.y[a] = 0 
                    else:
                        can.y[a] = 0
                        can.fixed0.append(a.id)
                        
            else:
                
                #---LB is obtained from LP relaxation of OA MP
                CG_status,can.LB,yCG = self.CG(can)
               
                if CG_status == 'infeasible':
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by feasibility')                    
                    prune = True
                    
                elif can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True                    
                   
                elif CG_status == 'integral':
                    #if self.params.PRINT_BB_INFO:
                    #    print('--> prune by integrality')

                    #---CG solution is integral and better than UB ==> runSO
                    #prune = True
                    runSO = True
                    runUE = True #optional needs to be there if interdiction cuts
                    
                    if self.params.useInterdictionCuts:
                        runUE = True
                        
                    for a in self.network.links2:
                        can.y[a] = round(yCG[a])
                        

                          
            if runSO:

                if not self.params.useAONcuts:
                    if self.ydict.hasSO(can.y) == True:
                        sotstt = self.ydict.getSO(can.y)
                        if self.params.PRINT_BB_INFO:
                            print('--> has SO') 


                             

                    else:
                        t0_TAP = time.time()
                        sotstt = self.network.tapas('SO',can.y)             
                        self.ydict.insertSO(can.y, sotstt)
                        self.OAcuts.append(self.getOAcut())
                        self.rt_TAP += time.time() - t0_TAP                        
                        self.nSO += 1
                    
                if self.params.useAONcuts:
                    t0_TAP = time.time()
                    self.network.setFlows(self.xvec)
                    sotstt = self.network.setAON('SO', can.y) # change to AON
                    self.rt_TAP += (time.time() - t0_TAP)
                    self.OAcuts.append(self.getOAcut())  
                    
                    # !!! can be memoized      
                    
                    #print('solved SO',sotstt)
                    
                if sotstt < self.UB_SO_DNDP:
                    self.UB_SO_DNDP = sotstt
                    if self.params.PRINT_BB_INFO:
                        print('--> update UB SO DNDP, prune nodes?')                        
                 
            if runUE:
                
                if self.ydict.hasUE(can.y) == True:
                    can.UB = self.ydict.getSO(can.y)
                    if self.params.PRINT_BB_INFO:
                        print('--> has UE')                    
                
                #---solve UE TAP to get UB 
                t0_TAP = time.time()
                #can.UB = self.network.tapas('UE',can.y)
                can.UB = self.network.tapas_ubstop('UE',can.y, self.UB)
                self.ydict.insertUE(can.y, can.UB)
                self.rt_TAP += time.time() - t0_TAP
                self.nUE += 1     
                
                if self.params.PRINT_BB_INFO:
                    print('TSTT UE',can.UB)
                
                if can.UB < self.UB:            
                    self.UB = can.UB
                    self.yopt = can.y
                    if self.params.PRINT_BB_INFO:
                        print('--> update UB')
                    
                    for n in self.BB_nodes:                    
                        if n.active == True and n.LB >= self.UB:
                            n.active = False
                            
            #---add global interdiction cuts
            if self.params.useInterdictionCuts and (runSO or runUE):
                self.yvec.append(can.y)

            if prune == False:
                fixed = can.fixed0 + can.fixed1
                free = [a.id for a in self.network.links2 if a.id not in fixed]
                
                #---first look for unfixed fractional variable                
                frac = [a.id for a in can.frac if a.id in free]
                if len(frac) > 0:
                    #print('--> branching on unfixed and fractional variable')
                    fsorted = sorted(frac, key = lambda ele: can.score[ele], reverse = True)
                    can.ybr = fsorted[0]
                    self.branch_frac(can)
                
                #---else look for unfixed variable
                else:
                    print('--> branching on unfixed variable')
                    fsorted = sorted(free, key = lambda ele: can.score[ele], reverse = True)                    
                    can.ybr = fsorted[0]
                    self.branch_unfixed(can)
                    
                #if self.params.PRINT_BB_INFO:                    
                #    for a in self.network.links2:
                #        if a.id == can.ybr:
                #            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can.ybr))                
                
            can.active = False
            candidates = self.getCandidates()
               
            if len(candidates) == 0:
                conv = True
                self.LB = self.UB
                self.gap = 0.0
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by inspection')
                
            else:
                self.LB = self.getLB(candidates)
                self.gap = self.getGap()
            
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (after): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))            
            
            if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
                print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (self.nBB,self.nSO,self.nUE,self.LB,self.UB,100*self.gap))
            
            if self.gap <= self.params.BB_tol:
                conv = True
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by optimality gap')
            
            if (time.time() - self.t0) >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print('--> time limit exceeded')
                break
            
            self.nBB += 1
 
        self.rt = time.time() - self.t0

        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,self.rt,self.nBB,self.nSO,self.nUE,self.UB,100*self.gap))
            print(self.rt_TAP)
            print(self.rt_RMP)
            print(self.rt_pricing)
            print(self.yopt)
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---BPC_singleTree end---')
        