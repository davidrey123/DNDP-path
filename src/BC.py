import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BC:

    #---instantiate a BC algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.INT_tol = 1e-4
        self.inf = 1e+9
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.t0 = 0.0
        self.M = self.inf
        self.rootNodeLB = 0
        
        self.nBB = 0
        self.nOA = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_LP = 0.0
        self.rt_rootNode = 0.0        
        
        self.LP = None
        self.cntUnscaledInf = 0 
        self.nIcuts = 0
        self.nOAcuts = 0
        self.nOABcuts = 0
        self.nVFcuts = 0
        self.nInitYvec = len(self.network.links2)
        
        n = BB_node.BB_node(self.network, 0, 0, self.LB, self.inf, [], [], False)
        self.BB_nodes.append(n)
        
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
    
    def getOAcuts(self):
        
        #---determine OA of TSTT function based on link flows and add link-based OA cuts to lp
        
        for a in self.network.links:
            OAcut = {}
            
            if a.y == 1:
                
                addOAcut = True
                for cut in a.OAcuts: 
                    
                    #---check if a.x or a.getTravelTime(a.x,'SO') is sufficiently different from existing OA cuts        
                    if abs(a.x - cut['x']) <= self.params.OAcut_tol*a.C:
                        addOAcut = False
                        break
                
                if addOAcut == True:
                    
                    OAcut['x'] = a.x
                    OAcut['a'] = a.getTravelTime(a.x,'SO')
                    OAcut['b'] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
                                                  
                    a.OAcuts.append(OAcut)
                    
                    #---add OA cut to lp
                    self.lp.add_constraint(self.lp.mu[a] >= self.lp.x[a]*OAcut['a'] + OAcut['b'])
                    self.nOAcuts += 1
                
    def getOABcuts(self):
        
        #---determine OA of Beckmann function based on (UE) link flows and add link-based OAB cuts to lp
        
        for a in self.network.links:
            OABcut = {}
            
            if a.y == 1:
                
                addOABcut = True
                for cut in a.OABcuts: 
                    
                    #---check if a.x is sufficiently different from existing OAB cuts        
                    if abs(a.x - cut['x']) <= self.params.OAcut_tol*a.C:
                        addOABcut = False
                        break
                
                if addOABcut == True:
                    
                    OABcut['x'] = a.x
                    OABcut['a'] = a.getTravelTime(a.x,'UE')
                    OABcut['b'] = a.getPrimitiveTravelTime(a.x) - a.x * a.getTravelTime(a.x,'UE')
                                                  
                    a.OABcuts.append(OABcut)
                    
                    #---add OAB cut to lp
                    self.lp.add_constraint(self.lp.muB[a] >= self.lp.x[a]*OABcut['a'] + OABcut['b'])
                    self.nOABcuts += 1
                               
                    
    def getVFcut1(self, beck):
        
        self.rmp.add_constraint(sum(self.rmp.muB[a] for a in self.network.links) <= beck + sum(self.M*(1 - self.rmp.y[a]) for a in self.network.links2 if a.x > 1e-4))
        self.nVFcuts += 1
    
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
            self.getOAcuts()                            
            self.nOA += 1        
            
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
        self.getOAcuts()                            
        self.nOA += 1
        
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
            self.getOAcuts()                            
            self.nOA += 1
            
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
        self.getOAcuts()                            
        self.nOA += 1
        
        y1 = {a:1 for a in self.network.links2}
        
        for a in self.network.links2:
            
            yLS = dict(y1)            
            yLS[a] = 0
                
            tstt = self.network.tapas('SO_OA_cuts',yLS)
            self.ydict.insertSO(yLS, tstt)       
            self.getOAcuts()                            
            self.nOA += 1
                    
        return yBest    
        
    def initOAcuts(self, type, can):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('Initializing OA cuts with',self.params.initOAheuristic)
            
        t0_OA = time.time()
        
        #---solve SO-TAP(y1) to initialize OA cuts and branching scores
        y1 = {a:1 for a in self.network.links2}
        
        self.network.tapas('SO_OA_cuts',y1)
        self.getOAcuts()
        self.nOA += 1
        
        for a in self.network.links2:
            can.score[a.id] = a.x * a.getTravelTime(a.x,'SO')         
        
        #---select heuristic to explore y vectors and initialize OA cuts
        if type == 'kBestKNP':
            yBest = self.kBestKNP()
            
        elif type == 'LocalSearchKNP':
            yBest = self.LocalSearchKNP()
            
        elif type == 'LocalSearchY1':
            yBest = self.LocalSearchY1()            
            
        self.rt_OA = time.time() - t0_OA
                
        #---solve UE-TAP(yBest) to get UB
        t0_TAP = time.time()
        can.UB = self.network.tapas('UE',yBest)
        self.ydict.insertUE(yBest,can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB
        
        if self.params.useValueFunctionCuts1 or self.params.useValueFunctionCuts2:                        
            beck = self.network.getBeckmannOFV()            
            self.getOABcuts()
        
            if self.params.useValueFunctionCuts1:
                self.getVFcut1(beck)
                
            if self.params.useValueFunctionCuts2:
                self.ydict.insertBeck(yBest, beck)
        
        #---add best y for which we ran UE
        if self.params.useInterdictionCuts:
            self.lp.add_constraint(sum(self.lp.y[a] + yBest[a] - 2*self.lp.y[a]*yBest[a] for a in self.network.links2) >= 1)
            self.nIcuts += 1
            
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
    
    def createLP(self):
                   
        lp = Model()
        
        lp.x = {a:lp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        lp.xc = {(a,s):lp.continuous_var(lb=0) for a in self.network.links for s in self.network.zones}
        lp.mu = {a:lp.continuous_var(lb=0) for a in self.network.links}
        lp.y = {a:lp.continuous_var(lb=0,ub=1) for a in self.network.links2}
        
        if self.params.useValueFunctionCuts1 or self.params.useValueFunctionCuts2:
            lp.muB = {a:lp.continuous_var(lb=0) for a in self.network.links}
  
        lp.add_constraint(sum(lp.y[a] * a.cost for a in self.network.links2) <= self.network.B)                   
            
        for a in self.network.links:
            lp.add_constraint(sum(lp.xc[(a,s)] for s in self.network.zones) == lp.x[a])
            #lp.add_constraint(lp.mu[a] >= lp.x[a] * a.t_ff)
                
        for i in self.network.nodes:          
            for s in self.network.zones:            
                
                if i.id == s.id:
                    dem = - sum(r.getDemand(s) for r in self.network.zones)                
                elif isinstance(i, type(s)) == True:
                    dem = i.getDemand(s)
                else:
                    dem = 0
                
                lp.add_constraint(sum(lp.xc[(a,s)] for a in i.outgoing) - sum(lp.xc[(a,s)] for a in i.incoming) == dem)
               
        for a in self.network.links2:            
            lp.add_constraint(lp.x[a] <= lp.y[a] * self.network.TD)
        
        lp.minimize(sum(lp.mu[a] for a in self.network.links))
        
        lp.parameters.threads = self.params.CPLEX_threads
        lp.parameters.timelimit = self.params.BB_timelimit
        lp.parameters.read.scale = -1 #---turns off data scaling in cplex. Using default (0) occasionally yields unscaled infeasibilities
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d' % (lp.number_of_variables,lp.number_of_constraints))
        
        self.lp = lp
        
    def addBranchCuts(self,can):
                       
        Bcuts0 = self.lp.add_constraints([self.lp.y[a] == 0 for a in self.network.links2 if a.id in can.fixed0])
        Bcuts1 = self.lp.add_constraints([self.lp.y[a] == 1 for a in self.network.links2 if a.id in can.fixed1])
        
        return Bcuts0,Bcuts1
    
    def removeBranchCuts(self,Bcuts0,Bcuts1):
                       
        self.lp.remove_constraints(Bcuts0)
        self.lp.remove_constraints(Bcuts1)                
        
    def solveLP(self,can):
        
        t0_lp = time.time()
        
        self.lp.solve(log_output=False)
        
        #print(self.lp.solve_details.status)      
        if self.lp.solve_details.status == 'optimal with unscaled infeasibilities':
            self.lp.parameters.read.scale = -1
            self.lp.solve(log_output=False)
            self.lp.parameters.read.scale = 0            
            self.cntUnscaledInf += 1
            
        if self.lp.solve_details.status == 'infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            OFV = self.lp.objective_value
            LP_status = self.lp.solve_details.status
            
            #if self.params.PRINT_BB_INFO:
            #    print('OFV: %.1f, lp_status: %s' % (OFV,lp_status))
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = self.lp.y[a].solution_value
                maxscore = 0
                for OAcut in a.OAcuts:
                    maxscore = max(self.lp.x[a].solution_value * OAcut['a'] + OAcut['b'],maxscore)                    
                can.score[a.id] = self.lp.x[a].solution_value * maxscore
            
            #---set link flows to generate OA cuts
            for a in self.network.links:
                a.x = round(self.lp.x[a].solution_value,self.params.rd) 
                
            #---check LP solution integrality
            LP_status, can.frac = self.checkIntegral(yopt)
            
            if self.params.PRINT_BB_INFO:
                if LP_status != 'fractional':
                    print('LP status',LP_status)                
            
        self.rt_LP += (time.time() - t0_lp)
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, nOAcuts: %d, nIcuts: %d, cplexTime: %.1f, lpTime: %.1f' % (self.lp.number_of_variables,self.lp.number_of_constraints,self.nOAcuts,self.nIcuts,self.lp.solve_details.time,(time.time() - t0_lp)))
        
        return LP_status,OFV,yopt

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---')
        
        self.network.resetTapas()
 
        self.t0 = time.time()
        
        #---initialize lp
        self.createLP()
    
        #---initialize OA cuts
        self.initOAcuts(self.params.initOAheuristic,self.BB_nodes[0])
        
        if self.params.useValueFunctionCuts1 or self.params.useValueFunctionCuts2:
            
            y0 = {a:0 for a in self.network.links2}
            tstt = self.network.tapas('UE',y0)
            self.ydict.insertUE(y0, tstt)
            beck = self.network.getBeckmannOFV()
            self.ydict.insertBeck(y0, beck)
            self.getOABcuts()
            self.M = beck          
    
        conv = False
        while conv == False:                    
             
            can = self.nodeSelection('bestBound',self.getCandidates())
            status = can.check()
        
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
     
            prune = False
            runOA = False
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
                    
                    #---LB is obtained from LP relaxation of OA MP
                    LP_status,can.LB,yLP = self.solveLP(can)
                   
                    if LP_status == 'infeasible':
                        if self.params.PRINT_BB_INFO:
                            print('--> LP infeasible - prune by feasibility')
                        prune = True
                        
                    elif can.LB >= self.UB:
                        if self.params.PRINT_BB_INFO:
                            print('--> prune by bounding')
                        prune = True
                        
                    else:
                        #---LP solution is feasible ==> runOA
                        runOA = True                    
                       
                    if LP_status == 'integral':
                        if self.params.PRINT_BB_INFO:
                            print('--> LP integral')
                            
                        for a in self.network.links2:
                            can.y[a] = round(yLP[a])                        
                        
                        if self.params.runUEifCGIntegral:
                            runUE = True
                            
                    #---remove Branch cuts
                    self.removeBranchCuts(Bcuts0,Bcuts1)     
                
                else:
                    if self.params.PRINT_BB_BASIC:
                        print('already solved')
                    if can.LB >= self.UB:
                        if self.params.PRINT_BB_BASIC:
                            print('--> prune by bounding2')
                        if self.params.PRINT_BB_INFO:
                            print('--> prune by bounding2')
                        prune = True                    
                          
            if runOA:

                t0_OA = time.time()
                
                if self.params.solveSO:

                    #---look for can.y in ydict (hash), solve SO-TAP if not found
                    if self.ydict.hasSO(can.y) == True:
                        sotstt = self.ydict.getSO(can.y)
                        if self.params.PRINT_BB_INFO:
                            print('--> has SO') 
    
                    else:                        
                        sotstt = self.network.tapas('SO_OA_cuts',can.y)                        
                        self.ydict.insertSO(can.y, sotstt)
                        self.getOAcuts()

                else:
                    
                    #---get OA cuts from LP x solution
                    self.getOAcuts()

                self.rt_OA += time.time() - t0_OA
                self.nOA += 1                                        
                 
            if runUE:
                
                t0_TAP = time.time()
                
                if self.ydict.hasUE(can.y) == True:
                    can.UB = self.ydict.getUE(can.y)
                    if self.params.PRINT_BB_INFO:
                        print('--> has UE')                    
                
                else:
                    #---solve UE TAP to get UB                    
                    can.UB = self.network.tapas('UE',can.y)
                    self.ydict.insertUE(can.y, can.UB)
                    self.nUE += 1  
                    
                    if self.params.useValueFunctionCuts1 or self.params.useValueFunctionCuts2:                        
                        beck = self.network.getBeckmannOFV()
                        self.ydict.insertBeck(can.y, beck)
                        self.getOABcuts()
                        
                        if self.params.useValueFunctionCuts1:
                            self.getVFcut1(beck)
                        
                self.rt_TAP += time.time() - t0_TAP                        
                
                if self.params.PRINT_BB_INFO:
                    print('TSTT UE - UB: %.1f, time: %.1f' % (can.UB,time.time() - t0_TAP))
                
                if can.UB < self.UB:            
                    self.UB = can.UB
                    self.yopt = can.y
                    if self.params.PRINT_BB_INFO:
                        print('--> update UB')
                    
                    for n in self.BB_nodes:                    
                        if n.active == True and n.LB >= self.UB:
                            n.active = False
                            
            #---add interdiction cuts
            if self.params.useInterdictionCuts and runUE:
                self.lp.add_constraint(sum(self.lp.y[a] + can.y[a] - 2*self.lp.y[a]*can.y[a] for a in self.network.links2) >= 1)
                self.nIcuts += 1

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
                    print(fixed,frac)
                    print(can.score)
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
                print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (self.nBB,self.nOA,self.nUE,self.LB,self.UB,100*self.gap))
            
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
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,self.rt,self.nBB,self.nOA,self.nUE,self.UB,100*self.gap))
            print(self.rt_TAP)
            print(self.rt_OA)            
            print(self.rt_LP)
            print(self.yopt)
            print(self.cntUnscaledInf)
            
            if self.params.useValueFunctionCuts1:
                print(self.nOABcuts,self.nVFcuts)
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')
        