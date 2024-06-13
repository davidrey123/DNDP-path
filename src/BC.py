import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BC:

    #---instantiate FS NETS's branch-and-bound algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.inf = 1e+9
        self.INT_tol = 1e-4
        self.OA_tol = 1e-4
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt_TAP = 0.0
        self.rt_LP = 0.0        
        
        self.OAcuts = []
        
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
        
        for yv in can.yvec:
            knp.add_constraint(sum(knp.y[a] + yv[a] - 2*knp.y[a]*yv[a] for a in self.network.links2) >= 1)        
        
        if type == 'capacity':
            knp.maximize(sum(knp.y[a] * a.C for a in self.network.links2))
        elif type == 'x':
            knp.maximize(sum(knp.y[a] * a.x for a in self.network.links2))                        
        else:
            print('unknown type')
            
        knp.solve(log_output=False)
            
        if knp.solve_details.status == 'infeasible' or knp.solve_details.status == 'integer infeasible':
            print('infeasible instance')
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
                    
            can.yvec.append(yKNP)      
            
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
    
    def lp_link(self,can):
        
        #---to do: recode such that lp is setup once and only cuts are added at each OA iteration
                
        t0_LP = time.time()
        lp = Model()    
        
        lp.x = {a:lp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        lp.xc = {(a,s):lp.continuous_var() for a in self.network.links for s in self.network.zones}
        lp.y = {a:lp.continuous_var(lb=0,ub=1) for a in self.network.links2}
        lp.mu = lp.continuous_var(lb=can.LB)
        
        lp.add_constraint(sum(lp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        for a in self.network.links2:
            lp.add_constraint(lp.x[a] <= lp.y[a] * self.network.TD)
            
        for a in self.network.links:
            lp.add_constraint(sum(lp.xc[(a,s)] for s in self.network.zones) == lp.x[a])
                
        for i in self.network.nodes:                    
            for s in self.network.zones:            
                
                if i.id == s.id:
                    dem = - sum(r.getDemand(s) for r in self.network.zones)                
                elif isinstance(i, type(s)) == True:
                    dem = i.getDemand(s)
                else:
                    dem = 0
                
                lp.add_constraint(sum(lp.xc[(a,s)] for a in i.outgoing) - sum(lp.xc[(a,s)] for a in i.incoming) == dem)
                        
        for OAcut in self.OAcuts:
            #---OA cuts
            lp.add_constraint(lp.mu >= sum(lp.x[a]*OAcut['a'][a] + OAcut['b'][a] for a in self.network.links))
            
        for yv in can.yvec:        
            #---Interdiction Cuts
            lp.add_constraint(sum(lp.y[a] + yv[a] - 2*lp.y[a]*yv[a] for a in self.network.links2) >= 1)
            
        #---Branch cuts    
        for a in self.network.links2:
            
            if a.id in can.fixed0:
                lp.add_constraint(lp.y[a] == 0)
                               
            if a.id in can.fixed1:
                lp.add_constraint(lp.y[a] == 1)
        
        lp.minimize(lp.mu)
        
        lp.solve(log_output=False)
                
        #print(lp.solve_details.time,(time.time() - t0_lp))
        self.rt_LP += (time.time() - t0_LP)
            
        if lp.solve_details.status == 'infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            LP_OFV = lp.objective_value
                         
            yLP = {}
            for a in self.network.links2:
                yLP[a] = lp.y[a].solution_value
                
            LP_status, can.frac = self.checkIntegral(yLP)
            
            #print('XXX',LP_OFV,LP_status,len(self.OAcuts))            
                
            return LP_status,LP_OFV,yLP

    
    def OA_lp_link(self,can):
        nOA = 0
        UB_OA = self.inf
        LB_OA = 0.0 #---value to be returned that will serve as the LB of the BB node
        yOA = {}
        
        t0_OA = time.time()
        
        conv = False
        while conv == False:            
                        
            LP_status,LP_OFV,yLP = self.lp_link(can)
            
            can.LB = LP_OFV
        
            if LP_status != 'integral':
                
                conv = True
                
                if LP_status == 'infeasible':
                    if nOA == 0:
                        print('==========WARNING: LP infeasible at first iteration===================')
                        #---lp must be feasible at first iteration since interdiction cuts are reset at each BB node
                        #   and budget constraint violations are checked before executing OA_link
                    
                    #---search is complete and UB_OA is the optimal OFV                    
                    LB_OA = max(can.LB,UB_OA)
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> convergence by feasibility (due to interdiction cuts)')
                    return nOA,LB_OA,yOA,LP_status
                
                else:
                    #---solution is fractional: stop OA an return LP_OFV as BB node LB                   
                    LB_OA = LP_OFV                    
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> fractional solution')                        
                
                
            else:
                #---integral solution, generate OA cut and keep going
                        
                if LP_OFV >= self.UB:
                    #---search can be stopped and lp obj can be used as the LB of the BB node
                    conv = True
                    LB_OA = max(can.LB,LP_OFV)
                    if self.params.PRINT_BB_INFO:
                        print('-----------------------------------> convergence by bounding in OA_lp_link')
                    return nOA,LB_OA,yOA,LP_status
            
                #---add interdiction cut
                can.yvec.append(yLP)
                
                #---solve SO-TAP to get OA cut
                if self.ydict.hasSO(yLP) == True:
                    tstt = self.ydict.getSO(yLP)
                
                else:
                    t0_TAP = time.time()
                    tstt = round(self.network.tapas('SO',yLP), 3)                
                    self.ydict.insertSO(yLP, tstt)
                    self.rt_TAP += time.time() - t0_TAP
                    self.OAcuts.append(self.getOAcut())
                    self.nSO += 1
            
                if tstt < UB_OA:
                    #if self.params.PRINT_BB_INFO:
                    #    print('-----------------------------------> update UB in OA_lp_link')
                    UB_OA = tstt
                    yOA = yLP                    
                
                
            for a in self.network.links2:
                can.score[a.id] = round(a.x * a.getTravelTime(a.x,'SO'), 3)
            
            gap = (UB_OA-LP_OFV)/UB_OA
            if self.params.PRINT_BB_INFO:
                print('-----------------------------------> %d\t%.1f\t%.1f\t%.2f%%' % (nOA,LP_OFV,UB_OA,100*gap))
            
            if gap <= self.OA_tol:
                #---search can be stopped and the min between lp obj and UB_OA can be used as the LB of the BB node
                conv = True
                LB_OA = max(can.LB,min(LP_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by optimality gap in OA_link')                    
                
            if (time.time() - t0_OA) >= self.params.BB_timelimit/2:
                #---search is stopped and the min between lp obj and UB_OA can be used as the LB of the BB node
                LB_OA = max(can.LB,min(LP_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> time limit exceeded in OA_link')
                break
            
            nOA += 1
        
            
        return nOA,LB_OA,yOA,LP_status
        

    def BB(self):
        
        print('---BC---')        
        
        self.network.resetTapas()
 
        t0 = time.time()
        
        #---initialize OA cuts
        self.initOAcuts(self.BB_nodes[0],10)        
    
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
                nOA, can.LB, yOA, LP_status = self.OA_lp_link(can)
                
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True 
                
                elif LP_status == 'integral':
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
                gap = 0.0
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by inspection')
                break
                
            else:
                self.LB = self.getLB(candidates)
                gap = self.getGap()
            
            if self.params.PRINT_BB_INFO:
                print('--> can (after): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))            
            
            print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (self.nBB,self.nSO,self.nUE,self.LB,self.UB,100*gap))
            
            if gap <= self.params.BB_tol:
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
 
        rt = time.time() - t0
 
        print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,rt,self.nBB,self.nSO,self.nUE,self.UB,100*gap))
        print(self.rt_TAP)
        print(self.rt_LP)
        print(self.yopt)
        
        return
    