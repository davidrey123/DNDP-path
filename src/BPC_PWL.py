import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BPC_PWL:

    #---instantiate a BPC algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.INT_tol = 1e-4
        self.CG_tol = 1e-4
        self.inf = 1e+9
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.UBlin = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.t0 = 0.0        
        self.M = self.inf
        self.rootNodeLB = 0
        
        self.nBB = 0
        self.nPWL = 100
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_RMP = 0.0
        self.rt_pricing = 0.0
        self.rt_rootNode = 0.0
        
        self.rmp = None
        self.tap = None
        self.cntUnscaledInf = 0 
        self.nInitYvec = len(self.network.links2)
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        
        n = BB_node.BB_node(self.network, 0, 0, self.LB, self.inf, [], [], False)
        self.BB_nodes.append(n)
        
        self.V = set([i for i in range(self.nPWL+1)])  
        self.Vp = self.V.difference({0})
        self.alpha = {(a,v):float() for a in self.network.links for v in self.V}    
        for a in self.network.links:
            cnt = 0
            step = self.network.TD/(len(self.V)-1)
            for v in self.V:
                self.alpha[(a,v)] = cnt*step
                cnt += 1     
                
        self.c = {a:a.t_ff*a.alpha/(a.C**a.beta) for a in self.network.links}                       
        
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
        return (self.UBlin - self.LB)/self.UBlin
    
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
    
    def knapsack(self, yvec):
        
        knp = Model()
        knp.y = {a:knp.binary_var() for a in self.network.links2}
        knp.add_constraint(sum(knp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        #---Interdiction cuts
        for yv in yvec:
            knp.add_constraint(sum(knp.y[a] + yv[a] - 2*knp.y[a]*yv[a] for a in self.network.links2) >= 1)
        
        knp.maximize(sum(knp.y[a] * a.x for a in self.network.links2))
            
        knp.solve(log_output=False)
            
        if knp.solve_details.status == 'infeasible' or knp.solve_details.status == 'integer infeasible':
            print('infeasible instance?')
            return {}
        else:            
            yKNP = {}
            for a in self.network.links2:
                yKNP[a] = round(knp.y[a].solution_value)
                
            return yKNP
        
    def kBestKNP(self):
        
        tsttBest = self.inf
        yBest = None
        yvec = []
        for n in range(self.nInitYvec):
            yKNP = self.knapsack(yvec)
                       
            tstt = self.network.tapas('SO_OA_cuts',yKNP)
            self.ydict.insertSO(yKNP, tstt)  
            
            if tstt < tsttBest:
                tsttBest = tstt
                yBest = yKNP                
                
            #---interdict yKNP to get different y
            yvec.append(yKNP)
        
        return yBest
    
    def LocalSearchKNP(self):
        
        yKNP = self.knapsack([])                       
        tstt = self.network.tapas('SO_OA_cuts',yKNP)
        self.ydict.insertSO(yKNP, tstt)        
        
        yBest = yKNP
        tsttBest = tstt
        
        for a in self.network.links2:
            
            yLS = dict(yKNP)
            
            if yLS[a] == 1:
                yLS[a] = 0
                
            else:
                yLS[a] = 1            
                
            tstt = self.network.tapas('SO_OA_cuts',yLS)
            self.ydict.insertSO(yLS, tstt)
            
            yCost = sum(b.cost*yLS[b] for b in yLS)
            
            if yCost <= self.network.B:
            
                if tstt < tsttBest:
                    tsttBest = tstt
                    yBest = yLS
        
        return yBest
    
    def LocalSearchY1(self):
        
        yBest = self.knapsack([])                       
        tstt = self.network.tapas('SO_OA_cuts',yBest)
        self.ydict.insertSO(yBest, tstt)   
        
        y1 = {a:1 for a in self.network.links2}
        
        for a in self.network.links2:
            
            yLS = dict(y1)            
            yLS[a] = 0
                
            tstt = self.network.tapas('SO_OA_cuts',yLS)
            self.ydict.insertSO(yLS, tstt)  
                    
        return yBest    
        
    def init(self, type, can):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('Initialization',self.params.initOAheuristic)
            
        #---solve SO-TAP(y1) to initialize OA cuts and branching scores
        y1 = {a:1 for a in self.network.links2}
        
        self.network.tapas('SO_OA_cuts',y1)        
        
        for a in self.network.links2:
            can.score[a.id] = a.x * a.getTravelTime(a.x,'SO')         
        
        #---select heuristic to explore y vectors and initialize OA cuts
        if type == 'kBestKNP':
            yBest = self.kBestKNP()
            
        elif type == 'LocalSearchKNP':
            yBest = self.LocalSearchKNP()
            
        elif type == 'LocalSearchY1':
            yBest = self.LocalSearchY1()            
                
        #---solve UE-TAP(yBest) to get UB
        t0_TAP = time.time()
        #can.UB = self.network.tapas('UE',yBest)
        can.UB = self.CGTAP(yBest)
        self.ydict.insertUE(yBest,can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UBlin = can.UB
        self.yopt = yBest    
        
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
                    new += 1
                
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('#paths',new)
    
    def pricing(self, can):
        
        t0_pricing = time.time()
        
        minrc = self.inf
        
        for a in self.network.links:
            a.dual = can.duals['link'][a]
            if abs(a.dual) <= self.CG_tol:
                a.dual = 0.0
                
            #print('a,a.dual',a.start.id,a.end.id,a.dual)
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                    
                    rc = - can.duals['dem'][(r,s)] + s.cost
                    
                    if rc < - self.CG_tol:
                        
                        p = self.network.trace(r,s)
                        self.paths[r][s].append(p)
                        
                        #---add new path var                                                
                        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
                        self.tap.h[p] = self.tap.continuous_var(lb=0)
                        
                        #---update constraints
                        self.rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).lhs.add_term(self.rmp.h[p], 1)
                        self.tap.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).lhs.add_term(self.tap.h[p], 1)
                        
                        for a in self.network.links:
                            if a in p.links:
                                self.rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).lhs.add_term(self.rmp.h[p], -1)
                                self.tap.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).lhs.add_term(self.tap.h[p], -1)
                    
                    if rc < minrc:
                        minrc = rc
                
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
    
    def createRMP(self):
                   
        rmp = Model()
        
        rmp.x = {a:rmp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        rmp.h = {p:rmp.continuous_var(lb=0) for p in self.getPaths()}
        rmp.y = {a:rmp.continuous_var(lb=0,ub=1) for a in self.network.links2}
        
        rmp.add_constraint(sum(rmp.y[a] * a.cost for a in self.network.links2) <= self.network.B)                   
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    rmp.add_constraint(sum(rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
    
        for a in self.network.links:
            rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))
               
        for a in self.network.links2:            
            rmp.add_constraint(rmp.x[a] <= rmp.y[a] * self.network.TD)

        rmp.ll = {(a,v): rmp.continuous_var(lb=0) for a in self.network.links for v in self.V}
        rmp.lr = {(a,v): rmp.continuous_var(lb=0) for a in self.network.links for v in self.V}                
        for a in self.network.links:
            rmp.add_constraint(rmp.x[a] == sum(rmp.ll[(a,v)]*self.alpha[(a,v-1)] + rmp.lr[(a,v)]*self.alpha[(a,v)] for v in self.Vp))      
            rmp.add_constraint(sum(rmp.ll[(a,v)] + rmp.lr[(a,v)] for v in self.V) == 1)  
        
        rmp.minimize(sum(a.t_ff*rmp.x[a] + self.c[a]*sum(rmp.ll[(a,v)]*(self.alpha[(a,v-1)]**(a.beta+1)) + rmp.lr[(a,v)]*(self.alpha[(a,v)]**(a.beta+1)) for v in self.Vp) for a in self.network.links))
            
        rmp.parameters.threads = self.params.CPLEX_threads
        rmp.parameters.timelimit = self.params.BB_timelimit
        rmp.parameters.read.scale = -1 #---turns off data scaling in cplex. Using default (0) occasionally yields unscaled infeasibilities
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d' % (rmp.number_of_variables,rmp.number_of_constraints))
        
        self.rmp = rmp
        
    def addBranchCuts(self,can):
                       
        Bcuts0 = self.rmp.add_constraints([self.rmp.y[a] == 0 for a in self.network.links2 if a.id in can.fixed0])
        Bcuts1 = self.rmp.add_constraints([self.rmp.y[a] == 1 for a in self.network.links2 if a.id in can.fixed1])
        
        return Bcuts0,Bcuts1
    
    def removeBranchCuts(self,Bcuts0,Bcuts1):
                       
        self.rmp.remove_constraints(Bcuts0)
        self.rmp.remove_constraints(Bcuts1)
        
    def solveRMP(self,can):
        
        t0_RMP = time.time()
        
        self.rmp.solve(log_output=False)
        
        #print(self.rmp.solve_details.status)      
        if self.rmp.solve_details.status == 'optimal with unscaled infeasibilities':
            self.rmp.parameters.read.scale = -1
            self.rmp.solve(log_output=False)
            self.rmp.parameters.read.scale = 0            
            self.cntUnscaledInf += 1
            
        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            OFV = self.rmp.objective_value
            RMP_status = self.rmp.solve_details.status
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = self.rmp.y[a].solution_value
                a.y = round(self.rmp.y[a].solution_value,self.params.rd)                
                can.score[a.id] = self.rmp.x[a].solution_value #---used for selecting branching variable

            dual_link = {} 
            for a in self.network.links:
                #---set link flows to generate OA cuts and get duals
                a.x = round(self.rmp.x[a].solution_value,self.params.rd) 
                dual_link[a] = round(max(self.rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0),self.params.rd)
        
            dual_dem = {}
            for r in self.network.origins:
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        dual_dem[(r,s)] = round(max(self.rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0),self.params.rd)
                
            can.duals = {'link':dual_link,'dem':dual_dem}
            
        self.rt_RMP += (time.time() - t0_RMP)
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, nOAcuts: %d, nIcuts: %d, cplexTime: %.1f, RMPTime: %.1f' % (self.rmp.number_of_variables,self.rmp.number_of_constraints,self.nOA,self.rmp.solve_details.time,(time.time() - t0_RMP)))
        
        return RMP_status,OFV,yopt
    
    def CG(self, can):
        
        nCG = 0
        conv = False

        while conv == False:        
                
            RMP_status,OFV,yRMP = self.solveRMP(can)
                
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                break
    
            minrc = self.pricing(can)
            
            if self.params.PRINT_BB_INFO:
                npaths = len(self.getPaths())
                print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))
                            
            if minrc >= -self.CG_tol:
                conv = True            
            
            nCG += 1
            
        if conv == True:            
            
            CG_status, can.frac = self.checkIntegral(yRMP)
            
            if self.params.PRINT_BB_INFO:
                if CG_status != 'fractional':
                    print('CG status',CG_status)      
            
            return CG_status,OFV,yRMP
        
        else:
            return CG_status,self.inf,yRMP

    def BB(self):        
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---'+' with nPWL = '+str(self.nPWL))
        
        self.network.resetTapas()
 
        self.t0 = time.time()
    
        #---initialize paths
        self.initPaths()
        
        #---initialize RMP
        self.createRMP()
        self.createTAP()
            
        #---initialize OA cuts
        self.init(self.params.initOAheuristic,self.BB_nodes[0])
            
        conv = False
        while conv == False:                    
             
            can = self.nodeSelection('bestBound',self.getCandidates())
            status = can.check()
                
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
                
            prune = False
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
                
                if can.solved == False:
                
                    #---add Branch cuts
                    Bcuts0,Bcuts1 = self.addBranchCuts(can)
                    
                    CG_status,OFV,yCG = self.CG(can)
                    can.LB = max(OFV,can.LB)
                   
                    if CG_status == 'infeasible':
                        if self.params.PRINT_BB_INFO:
                            print('--> CG infeasible - prune by feasibility')
                        prune = True
                        
                    elif can.LB >= self.UBlin:
                        if self.params.PRINT_BB_INFO:
                            print('--> prune by bounding')
                        prune = True
                       
                    if CG_status == 'integral':
                        if self.params.PRINT_BB_INFO:
                            print('--> CG integral')
                            
                        for a in self.network.links2:
                            can.y[a] = round(yCG[a])
                            
                        if self.params.runUEifCGIntegral:
                            runUE = True
                            
                    #---remove Branch cuts
                    self.removeBranchCuts(Bcuts0,Bcuts1)
                        
                else:
                    if self.params.PRINT_BB_INFO:
                        print('already solved')
                    if can.LB >= self.UBlin:
                        if self.params.PRINT_BB_INFO:
                            print('--> prune by bounding2')
                        prune = True
                 
            if runUE:
                
                t0_TAP = time.time()
                
                if self.ydict.hasUE(can.y) == True:
                    can.UB = self.ydict.getUE(can.y)
                    if self.params.PRINT_BB_INFO:
                        print('--> has UE')                    
                
                else:
                    #---solve UE TAP lin to get UB                     
                    can.UB = self.CGTAP(can.y)
                    self.ydict.insertUE(can.y, can.UB)
                    self.nUE += 1  
                        
                self.rt_TAP += time.time() - t0_TAP
                
                if self.params.PRINT_BB_INFO:
                    print('TSTT UE - UB: %.1f, time: %.1f' % (can.UB,time.time() - t0_TAP))
                
                if can.UB < self.UBlin:            
                    self.UBlin = can.UB
                    self.yopt = can.y
                    if self.params.PRINT_BB_INFO:
                        print('--> update UB lin')
                    
                    for n in self.BB_nodes:                    
                        if n.active == True and n.LB >= self.UBlin:
                            n.active = False   
            
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
                    #print('--> branching on unfixed variable')
                    fsorted = sorted(free, key = lambda ele: can.score[ele], reverse = True)                    
                    can.ybr = fsorted[0]
                    self.branch_unfixed(can)
                    
                if self.params.PRINT_BB_INFO:

                    for a in self.network.links2:
                        if a.id == can.ybr:
                            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can.ybr))

            if self.nBB == 0:
                self.rootNodeLB = can.LB
                self.rt_rootNode = time.time() - self.t0                
                
            can.active = False
            candidates = self.getCandidates()
               
            if len(candidates) == 0:
                conv = True
                self.LB = self.UBlin
                self.gap = 0.0
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by inspection')
                
            else:
                self.LB = self.getLB(candidates)
                self.gap = self.getGap()     
            
            if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
                npaths = len(self.getPaths())
                elapsed_time = time.time() - self.t0
                print('%d\t%d\t%d\t%.1f\t%.1f\t%.2f%%\t%.1f' % (self.nBB,npaths,self.nUE,self.LB,self.UBlin,100*self.gap,elapsed_time))
            
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

        self.UB = self.network.tapas('UE',self.yopt)

        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            npaths = len(self.getPaths())
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%\t%.1f' % (conv,self.rt,self.nBB,npaths,self.nUE,self.UBlin,100*self.gap,self.UB))
            print("Time TAP: %.2f" % self.rt_TAP)
            print("Time RMP: %.2f" % self.rt_RMP)
            print("Time Prc: %.2f" % self.rt_pricing)
            print(self.yopt)
            if self.cntUnscaledInf > 0:
                print(self.cntUnscaledInf)
                
        if self.params.PRINT_LOG:
            logfile.close()
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')
        
    def createTAP(self):
        tap = Model()
        
        tap.x = {a:tap.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        tap.h = {p:tap.continuous_var(lb=0) for p in self.getPaths()}                       
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    tap.add_constraint(sum(tap.h[p] for p in self.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
    
        for a in self.network.links:
            tap.add_constraint(tap.x[a] - sum(tap.h[p] for p in self.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))   
    
        tap.ll = {(a,v): tap.continuous_var(lb=0) for a in self.network.links for v in self.V}
        tap.lr = {(a,v): tap.continuous_var(lb=0) for a in self.network.links for v in self.V}                
        for a in self.network.links:
            tap.add_constraint(tap.x[a] == sum(tap.ll[(a,v)]*self.alpha[(a,v-1)] + tap.lr[(a,v)]*self.alpha[(a,v)] for v in self.Vp))      
            tap.add_constraint(sum(tap.ll[(a,v)] + tap.lr[(a,v)] for v in self.V) == 1)  

        tap.minimize(sum(a.t_ff*tap.x[a] + (self.c[a]/(a.beta+1))*sum(tap.ll[(a,v)]*(self.alpha[(a,v-1)]**(a.beta+1)) + tap.lr[(a,v)]*(self.alpha[(a,v)]**(a.beta+1)) for v in self.Vp) for a in self.network.links))                        
        
        tap.parameters.threads = 1
        tap.parameters.simplex.display = 0        
        
        self.tap = tap
        
    def solveTAP(self,y):
        
        t0_TAP = time.time()
        
        self.tap.solve(log_output=False)   
        OFV = self.tap.objective_value

        dual_link = {} 
        for a in self.network.links:
            a.x = round(self.tap.x[a].solution_value,self.params.rd) 
            dual_link[a] = round(max(self.tap.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0),self.params.rd)
    
        dual_dem = {}
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    dual_dem[(r,s)] = round(max(self.tap.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0),self.params.rd)
            
        dualsTAP = {'link':dual_link,'dem':dual_dem}
        
        xopt = {a:self.tap.solution.get_value(self.tap.x[a]) for a in self.network.links}
        llopt = {(a,v):self.tap.solution.get_value(self.tap.ll[(a,v)]) for a in self.network.links for v in self.V}
        lropt = {(a,v):self.tap.solution.get_value(self.tap.lr[(a,v)]) for a in self.network.links for v in self.V}
        UE_TSTT = 0
        topt = {}
        for a in self.network.links:
            topt[a] = (a.t_ff + self.c[a]*sum(llopt[(a,v)]*(self.alpha[(a,v-1)]**(a.beta)) + lropt[(a,v)]*(self.alpha[(a,v)]**(a.beta)) for v in self.Vp))        
            UE_TSTT += topt[a]*xopt[a]
            
        self.rt_TAP += (time.time() - t0_TAP)            
        
        return UE_TSTT,dualsTAP       
        
    def CGTAP(self, y):
        
        ycuts = self.tap.add_constraints([self.tap.x[a] == 0 for a in self.network.links2 if y[a] <= 1e-4])
        
        nCG = 0
        conv = False

        while conv == False:        
                
            OFV,dualsTAP = self.solveTAP(y)
    
            minrc = self.pricingTAP(dualsTAP)
            
            if minrc >= -self.CG_tol:
                conv = True
                
            if nCG >= 50:
                print('WARNING nCG TAP lin >= 50')
                conv = True                
            nCG += 1
            
        self.tap.remove_constraints(ycuts)            
            
        return OFV
        
    def pricingTAP(self, duals):
        
        t0_pricing = time.time()
        
        minrc = self.inf
        
        for a in self.network.links:
            a.dual = duals['link'][a]
            if abs(a.dual) <= self.CG_tol:
                a.dual = 0.0
                
            #print('a,a.dual',a.start.id,a.end.id,a.dual)
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                    
                    rc = - duals['dem'][(r,s)] + s.cost
                    
                    if rc < - self.CG_tol:
                        
                        p = self.network.trace(r,s)
                        self.paths[r][s].append(p) 
                        
                        #---add new path var
                        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
                        self.tap.h[p] = self.tap.continuous_var(lb=0)
                        
                        #---update constraints
                        self.rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).lhs.add_term(self.rmp.h[p], 1)
                        self.tap.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).lhs.add_term(self.tap.h[p], 1)
                        
                        for a in self.network.links:
                            if a in p.links:
                                self.rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).lhs.add_term(self.rmp.h[p], -1)
                                self.tap.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).lhs.add_term(self.tap.h[p], -1)
                    
                    if rc < minrc:
                        minrc = rc
                
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc        