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
        self.UB_SO_DNDP = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.t0 = 0.0
        
        if self.params.useValueFunctionCuts:
            #self.M = {a:a.getPrimitiveTravelTime(self.network.TD) for a in self.network.links2}
            
            y0 = {}
            for a in self.network.links2:
                y0[a] = 0
                
            self.network.tapas('UE',y0)
            OFV = self.network.getBeckmannOFV()
            self.M = {a:OFV for a in self.network.links2}
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_LP = 0.0   
        
        self.LP = None
        self.cntUnscaledInf = 0 
        self.yvec = [] #---for global interdiction cuts --- is it needed to store if directly adding to lp?
        self.nOAcuts = 0
        self.nOABcuts = 0
        self.nVFcuts = 0
        self.nInitKNP = len(self.network.links2)
        
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
                    
                    #---check if a.x is sufficiently different from existing OABB cuts        
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
                               
                    
    def getVFcut(self):
        
        self.lp.add_constraint(sum(self.lp.muB[a] for a in self.network.links) <= self.network.getBeckmannOFV() + sum(self.M[a]*(1 - self.lp.y[a]) for a in self.network.links2 if a.x > 1e-4)) 
        self.nVFcuts += 1
    
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
        
    def initOAcuts(self, can):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('Running knapsack heuristic with',self.nInitKNP,'round(s)')
        
        yKNP = {}
        for a in self.network.links2:
            yKNP[a] = 1
        
        t0_TAP = time.time()
        tstt = self.network.tapas('SO_OA_cuts',yKNP)
        self.ydict.insertSO(yKNP, tstt)
        self.getOAcuts()
        self.rt_TAP += (time.time() - t0_TAP)
        self.nSO += 1
       
        for a in self.network.links2:
            can.score[a.id] = a.x * a.getTravelTime(a.x,'SO') 
        
        yBest = None
        for n in range(self.nInitKNP):
            yKNP = self.knapsack('x',can)
                       
            t0_TAP = time.time()
            tstt = self.network.tapas('SO_OA_cuts',yKNP)
            self.ydict.insertSO(yKNP, tstt)
            self.rt_TAP += (time.time() - t0_TAP)
            self.getOAcuts()                            
            self.nSO += 1        
            
            if tstt < self.UB_SO_DNDP:
                self.UB_SO_DNDP = tstt
                yBest = yKNP                
                
            #---interdict to get different y
            self.yvec.append(yKNP)
                
        t0_TAP = time.time()
        can.UB = self.network.tapas('UE',yBest)
        self.ydict.insertUE(yBest, can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB
        
        if self.params.useValueFunctionCuts:
            self.getOABcuts()
            self.getVFcut()
        
        #---reset yvec then add best y for which we ran UE
        self.yvec = []
        if self.params.useInterdictionCuts:       
            self.yvec.append(yBest) #---superfluous?
            self.lp.add_constraint(sum(self.lp.y[a] + yBest[a] - 2*self.lp.y[a]*yBest[a] for a in self.network.links2) >= 1)
            
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
        
        if self.params.useValueFunctionCuts:
            lp.muB = {a:lp.continuous_var(lb=0) for a in self.network.links}
  
        lp.add_constraint(sum(lp.y[a] * a.cost for a in self.network.links2) <= self.network.B)                   
            
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
                
            #---check LP solution integrality
            LP_status, can.frac = self.checkIntegral(yopt)
            
            if self.params.PRINT_BB_INFO:
                if LP_status != 'fractional':
                    print('LP status',LP_status)                
            
        self.rt_LP += (time.time() - t0_lp)
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, nOAcuts: %d, nIcuts: %d, cplexTime: %.1f, lpTime: %.1f' % (self.lp.number_of_variables,self.lp.number_of_constraints,self.nOAcuts,len(self.yvec),self.lp.solve_details.time,(time.time() - t0_lp)))
        
        return LP_status,OFV,yopt

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---')
        
        self.network.resetTapas()
 
        self.t0 = time.time()
        
        #---initialize lp
        self.createLP()
    
        #---initialize OA cuts
        self.initOAcuts(self.BB_nodes[0])
    
        conv = False
        while conv == False:                    
             
            can = self.nodeSelection('bestBound',self.getCandidates())
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
                runUE = True
                
                #---condition to be tuned: idea is to get OA cuts early in the search then stop                
                #if self.nSO <= 5*nInitCuts: # **2 *5
                
                runSO = True
                 
                for a in self.network.links2:
                    if a.id in can.fixed1:
                        can.y[a] = 1
                    elif a.id in can.fixed0:
                        can.y[a] = 0 
                    else:
                        can.y[a] = 0
                        can.fixed0.append(a.id)
                        
            else:
                
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
                   
                elif LP_status == 'integral':
                    if self.params.PRINT_BB_INFO:
                        print('--> LP integral')

                    #---CG solution is integral and better than UB ==> runSO
                    runSO = True
                    
                    if self.params.runUEifCGIntegral:
                        runUE = True
                        
                    for a in self.network.links2:
                        can.y[a] = round(yLP[a])
                        
                #---remove Branch cuts
                self.removeBranchCuts(Bcuts0,Bcuts1)                        
                          
            if runSO:

                if self.ydict.hasSO(can.y) == True:
                    sotstt = self.ydict.getSO(can.y)
                    if self.params.PRINT_BB_INFO:
                        print('--> has SO') 

                else:
                    t0_TAP = time.time()
                    sotstt = self.network.tapas('SO_OA_cuts',can.y)                        
                    self.ydict.insertSO(can.y, sotstt)
                    self.getOAcuts()
                    self.rt_TAP += time.time() - t0_TAP
                    self.nSO += 1                                         
                 
            if runUE:
                
                if self.ydict.hasUE(can.y) == True:
                    can.UB = self.ydict.getUE(can.y)
                    if self.params.PRINT_BB_INFO:
                        print('--> has UE')                    
                
                else:
                    #---solve UE TAP to get UB 
                    t0_TAP = time.time()
                    can.UB = self.network.tapas('UE',can.y)
                    self.ydict.insertUE(can.y, can.UB)
                    self.rt_TAP += time.time() - t0_TAP
                    self.nUE += 1  
                    
                    if self.params.useValueFunctionCuts:
                        self.getOABcuts()
                        self.getVFcut()                    
                
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
                self.yvec.append(can.y)
                self.lp.add_constraint(sum(self.lp.y[a] + can.y[a] - 2*self.lp.y[a]*can.y[a] for a in self.network.links2) >= 1)

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
                    #print(fixed,frac)
                    #print(can.score)
                    for a in self.network.links2:
                        if a.id == can.ybr:
                            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can.ybr))                
                
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
            print(self.rt_LP)
            print(self.yopt)
            print(self.cntUnscaledInf)
            
            if self.params.useValueFunctionCuts:
                print(self.nOABcuts,self.nVFcuts)
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')
        