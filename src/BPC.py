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
        
        self.rmp = None
        self.cntUnscaledInf = 0 
        self.yvec = [] #---for global interdiction cuts --- is it needed to store if directly adding to RMP?
        self.nOAcuts = 0
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
    
    def getOAcuts(self):
        
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
                    
                    #---add OA cut to RMP
                    self.rmp.add_constraint(self.rmp.mu[a] >= self.rmp.x[a]*OAcut['a'] + OAcut['b'])
                    self.nOAcuts += 1
                
    
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
        tstt = self.network.tapas('SO_OA_cuts',yKNP)
        self.ydict.insertSO(yKNP, tstt)
        self.getOAcuts()
        self.rt_TAP += (time.time() - t0_TAP)
        self.nSO += 1
        
        self.xvec = self.network.getFlowMap() 
        
        for a in self.network.links2:
            can.score[a.id] = a.x * a.getTravelTime(a.x,'SO') 
        
        yBest = None
        for n in range(nKNP):
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
            #print(n,yKNP)
                
        t0_TAP = time.time()
        can.UB = self.network.tapas('UE',yBest)
        self.ydict.insertUE(yBest, can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB
        
        #---reset yvec then add best y for which we ran UE
        self.yvec = []
        if self.params.useInterdictionCuts:       
            self.yvec.append(yBest) #---superfluous?
            self.rmp.add_constraint(sum(self.rmp.y[a] + yBest[a] - 2*self.rmp.y[a]*yBest[a] for a in self.network.links2) >= 1)
            
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
                        self.paths[r][s].append(p) #---is it needed to store paths if directly adding to RMP?
                        
                        #---add new path var to RMP                                                
                        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
                        
                        #---update RMP constraints
                        self.rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).lhs.add_term(self.rmp.h[p], 1)
                        
                        for a in self.network.links:
                            if a in p.links:
                                self.rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).lhs.add_term(self.rmp.h[p], -1)
                        
                        new += 1
                    
                    if rc < minrc:
                        minrc = rc
            
        #if self.params.PRINT_BB_INFO:
        #    print('pricing',new,minrc)        
            
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
    
    def createRMP(self):
                   
        rmp = Model()
        
        rmp.x = {a:rmp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        rmp.h = {p:rmp.continuous_var(lb=0) for p in self.getPaths()}
        rmp.mu = {a:rmp.continuous_var(lb=0) for a in self.network.links}
        rmp.y = {a:rmp.continuous_var(lb=0,ub=1) for a in self.network.links2}
  
        rmp.add_constraint(sum(rmp.y[a] * a.cost for a in self.network.links2) <= self.network.B)                   
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    rmp.add_constraint(sum(rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
                    
                    #---disaggregate linking constraints: not working as expected....very slow?
                    #for a in self.network.links2:
                    #    rmp.add_constraint(sum(rmp.h[p] for p in self.paths[r][s] if a in p.links) <= rmp.y[a] * r.getDemand(s))
    
        for a in self.network.links:
            rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))           
               
        for a in self.network.links2:            
            rmp.add_constraint(rmp.x[a] <= rmp.y[a] * self.network.TD)
        
        rmp.minimize(sum(rmp.mu[a] for a in self.network.links))
        
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
            
            #if self.params.PRINT_BB_INFO:
            #    print('OFV: %.1f, RMP_status: %s' % (OFV,RMP_status))
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = self.rmp.y[a].solution_value
                maxscore = 0
                for OAcut in a.OAcuts:
                    maxscore = max(self.rmp.x[a].solution_value * OAcut['a'] + OAcut['b'],maxscore)                    
                can.score[a.id] = self.rmp.x[a].solution_value * maxscore

            dual_link = {} 
            for a in self.network.links:                    
                dual_link[a] = max(self.rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0)                    
        
            dual_dem = {}
            for r in self.network.origins:
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        dual_dem[(r,s)] = max(self.rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0)
                
            can.duals = {'link':dual_link,'dem':dual_dem}
            
        self.rt_RMP += (time.time() - t0_RMP)
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, nOAcuts: %d, nIcuts: %d, cplexTime: %.1f, RMPTime: %.1f' % (self.rmp.number_of_variables,self.rmp.number_of_constraints,self.nOAcuts,len(self.yvec),self.rmp.solve_details.time,(time.time() - t0_RMP)))
        
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
    
        #---initialize paths
        self.initPaths()
        
        #---initialize RMP
        self.createRMP()
    
        #---initialize OA cuts
        nInitCuts = round(len(self.network.links2))
        self.initOAcuts(self.BB_nodes[0],nInitCuts)    
    
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
                CG_status,can.LB,yCG = self.CG(can)            
               
                if CG_status == 'infeasible':
                    if self.params.PRINT_BB_INFO:
                        print('--> CG infeasible - prune by feasibility')
                    prune = True
                    
                elif can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True                    
                   
                elif CG_status == 'integral':
                    if self.params.PRINT_BB_INFO:
                        print('--> CG integral')

                    #---CG solution is integral and better than UB ==> runSO
                    runSO = True
                    
                    if self.params.runUEifCGIntegral:
                        runUE = True
                        
                    for a in self.network.links2:
                        can.y[a] = round(yCG[a])
                        
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
                self.rmp.add_constraint(sum(self.rmp.y[a] + can.y[a] - 2*self.rmp.y[a]*can.y[a] for a in self.network.links2) >= 1)

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
            print(self.rt_RMP)
            print(self.rt_pricing)
            print(self.yopt)
            print(self.cntUnscaledInf)
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')
        