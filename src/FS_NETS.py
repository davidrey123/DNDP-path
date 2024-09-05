import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class FS_NETS:

    #---instantiate FS NETS's branch-and-bound algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.inf = 1e+9        
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.UB_SO_DNDP = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.OA_tol = self.params.BB_tol
        self.t0 = 0.0        
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_MILP = 0.0        
        
        self.milp = None
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
            n0.score = can.score
            n0.y = can.y
        else:
            n1.score = can.score       
            n1.y = can.y
        
        return
    
    def addBranchCuts(self,can):
        
        #---Branch cuts
        Bcuts0 = self.milp.add_constraints([self.milp.y[a] == 0 for a in self.network.links2 if a.id in can.fixed0])
        Bcuts1 = self.milp.add_constraints([self.milp.y[a] == 1 for a in self.network.links2 if a.id in can.fixed1])
              
        return Bcuts0,Bcuts1
    
    def removeBranchCuts(self,Bcuts0,Bcuts1):
                       
        self.milp.remove_constraints(Bcuts0)
        self.milp.remove_constraints(Bcuts1)
    
    def getOAcut(self):
        OAcut = {'a':{}, 'b':{}}
        
        for a in self.network.links:
        
            if a.y == 1:                
                OAcut['a'][a] = a.getTravelTime(a.x,'SO')
                OAcut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
            
            else:
                OAcut['a'][a] = 0.0
                OAcut['b'][a] = 0.0
                
        #---OA cut to MILP
        self.milp.add_constraint(self.milp.mu >= sum(self.milp.x[a]*OAcut['a'][a] + OAcut['b'][a] for a in self.network.links))
    
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
        
    def initOAcuts(self, can):
        
        if self.params.PRINT_BB_INFO:
            print('Running knapsack heuristic (designed for 1 round only)')
        
        yKNP = {}
        for a in self.network.links2:
            yKNP[a] = 1
        
        t0_TAP = time.time()
        tstt = self.network.tapas('SO',yKNP) 
        self.ydict.insertSO(yKNP, tstt)
        self.rt_TAP += (time.time() - t0_TAP)
        self.getOAcut()
        self.nSO += 1
        
        yKNP = self.knapsack('x',can)
                    
        t0_TAP = time.time()
        tstt = self.network.tapas('SO',yKNP) 
        self.ydict.insertSO(yKNP, tstt)
        self.rt_TAP += (time.time() - t0_TAP)
        self.getOAcut()                           
        self.nSO += 1
                
        can.yvec.append(yKNP)
        
        if tstt < self.UB_SO_DNDP:
            self.UB_SO_DNDP = tstt
            self.yopt = yKNP
                
        t0_TAP = time.time()
        can.UB = self.network.tapas('UE',self.yopt)
        self.ydict.insertUE(self.yopt, can.UB)
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        self.UB = can.UB           
    
    def createMILP(self):        
        
        milp = Model()    
        
        milp.x = {a:milp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        milp.xc = {(a,s):milp.continuous_var() for a in self.network.links for s in self.network.zones}
        milp.y = {a:milp.binary_var() for a in self.network.links2}
        milp.mu = milp.continuous_var(lb=0,name='mu')
        
        milp.add_constraint(sum(milp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        for a in self.network.links2:
            milp.add_constraint(milp.x[a] <= milp.y[a] * self.network.TD)
            
        for a in self.network.links:
            milp.add_constraint(sum(milp.xc[(a,s)] for s in self.network.zones) == milp.x[a])
                
        for i in self.network.nodes:                    
            for s in self.network.zones:            
                
                if i.id == s.id:
                    dem = - sum(r.getDemand(s) for r in self.network.zones)                
                elif isinstance(i, type(s)) == True:
                    dem = i.getDemand(s)
                else:
                    dem = 0
                
                milp.add_constraint(sum(milp.xc[(a,s)] for a in i.outgoing) - sum(milp.xc[(a,s)] for a in i.incoming) == dem)
                
        milp.minimize(milp.mu)
                
        milp.parameters.mip.tolerances.mipgap = self.params.BB_tol
        milp.parameters.threads = self.params.CPLEX_threads
        milp.parameters.timelimit = self.params.BB_timelimit
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d' % (milp.number_of_variables,milp.number_of_constraints))
        
        self.milp = milp        
        
    def solveMILP(self,can):
        
        t0_MILP = time.time()
        
        self.milp.solve(log_output=False)
                
        if self.milp.solve_details.status == 'infeasible' or self.milp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            MILP_OFV = self.milp.objective_value
            MILP_status = self.milp.solve_details.status
             
            yopt = {}
            for a in self.network.links2:
                yopt[a] = round(self.milp.y[a].solution_value)
                
        self.rt_MILP += (time.time() - t0_MILP)
                
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, cplexTime: %.1f, MILPTime: %.1f' % (self.milp.number_of_variables,self.milp.number_of_constraints,self.milp.solve_details.time,(time.time() - t0_MILP)))
                
        return MILP_status,MILP_OFV,yopt    

    
    def OA_link(self,can):
        
        for n in self.BB_nodes:
            if n.id == can.parent:
                LB_OA = n.LB #---value to be returned that will serve as the LB of the BB node
                break
        
        nOA = 0
        UB_OA = self.inf
        yOA = {}
        Icuts = []
               
        
        conv_OA = False
        while conv_OA == False:
            
            MILP_status,MILP_OFV,yMILP = self.solveMILP(can)
            
            if MILP_status == 'infeasible':
                if nOA == 0:
                    print('==========WARNING: MILP infeasible at first iteration===================')
                    #---milp must be feasible at first iteration since interdiction cuts are reset at each BB node
                    #   and budget constraint violations are checked before executing OA_link
                
                #---search is complete and UB_OA is the optimal OFV
                conv_OA = True
                LB_OA = max(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by feasibility (due to interdiction cuts)')
                return nOA,LB_OA,yOA,MILP_status
                        
            if MILP_OFV >= self.UB:
                #---search can be stopped and milp obj can be used as the LB of the BB node
                conv_OA = True
                LB_OA = max(LB_OA,MILP_OFV)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by bounding in OA_link')
                return nOA,LB_OA,yOA,MILP_status
            
            #---update LB of mu
            #can.LB = MILP_OFV #---can.LB is temporarily updated during the OA loop but will be overwritten after exiting OA_link
            self.milp.get_var_by_name('mu').lb = MILP_OFV
            
            #---update interdiction cuts
            #can.yvec.append(yMILP) #---superfluous if directly adding cut to MILP?
            Icut = self.milp.add_constraint(sum(self.milp.y[a] + yMILP[a] - 2*self.milp.y[a]*yMILP[a] for a in self.network.links2) >= 1)
            Icuts.append(Icut)            
            
            if self.ydict.hasSO(yMILP) == True:
                tstt = self.ydict.getSO(yMILP)
                if self.params.PRINT_BB_INFO:
                    print('hasSO')
            
            else:
                t0_TAP = time.time()
                tstt = self.network.tapas('SO',yMILP)              
                self.ydict.insertSO(yMILP, tstt)
                self.rt_TAP += time.time() - t0_TAP
                self.getOAcut()
                self.nSO += 1
            
            if tstt < UB_OA:
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> update UB in OA_link')
                UB_OA = tstt
                yOA = yMILP
                
                for a in self.network.links2:
                    can.score[a.id] = a.x * a.getTravelTime(a.x,'SO')
            
            gap = (UB_OA-MILP_OFV)/UB_OA
            if self.params.PRINT_BB_INFO:
                print('-----------------------------------> %d\t%.1f\t%.1f\t%.2f%%' % (nOA,MILP_OFV,UB_OA,100*gap))
            
            if gap <= self.OA_tol:
                #---search can be stopped and the min between milp obj and UB_OA can be used as the LB of the BB node
                conv_OA = True
                LB_OA = max(LB_OA,min(MILP_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by optimality gap in OA_link')                    
                
            if (time.time() - self.t0) >= self.params.BB_timelimit:
                #---search is stopped and the min between milp obj and UB_OA can be used as the LB of the BB node
                LB_OA = max(LB_OA,min(MILP_OFV,UB_OA))
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> time limit exceeded in OA_link')
                break
            
            nOA += 1
            
        #---reset LB of mu
        self.milp.get_var_by_name('mu').lb = LB_OA
        
        #---remove local interdiction cuts
        for Icut in Icuts:
            self.milp.remove_constraint(Icut)        
            
        return nOA,LB_OA,yOA,MILP_status
        

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---')   
        
        self.network.resetTapas()
 
        self.t0 = time.time()
        
        #---initialize MILP
        self.createMILP()
        
        #---initialize OA cuts
        self.initOAcuts(self.BB_nodes[0])
    
        conv = False
        while conv == False:                   
             
            can = self.nodeSelection(self.getCandidates())
            status = can.check()
            
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
     
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
                if can.solved == False:
                    
                    #---add Branch cuts
                    Bcuts0,Bcuts1 = self.addBranchCuts(can)
                    
                    #---LB is obtained from OA algorithm
                    nOA,can.LB,can.y,OA_status = self.OA_link(can)                    
                    
                    #---solve UE-TAP at root node to get initial UB
                    if self.nBB == 0:
                        t0_TAP = time.time()
                        can.UB = self.network.tapas('UE',can.y)
                        self.rt_TAP += time.time() - t0_TAP
                        self.nUE += 1
                        
                        if can.UB < self.UB:            
                            self.UB = can.UB
                            self.yopt = can.y
                            if self.params.PRINT_BB_INFO:
                                print('--> update UB')
                                
                    #---remove Branch cuts
                    self.removeBranchCuts(Bcuts0,Bcuts1)                    
                 
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True
                     
                else:
                    prune = False
     
            if integral == True:
                if self.ydict.hasUE(yUB) == True:
                    can.UB = self.ydict.getUE(yUB)
                
                else:
                    t0_TAP = time.time()
                    can.UB = self.network.tapas('UE',yUB)
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
                fixed = can.fixed0 + can.fixed1
                free = [a.id for a in self.network.links2 if a.id not in fixed]            
                free_sorted = sorted(free, key = lambda ele: can.score[ele], reverse = True)
                can.ybr = free_sorted[0]
                
                #if self.params.PRINT_BB_INFO:                    
                #    for a in self.network.links2:
                #        print(a,can.score[a.id])
                #        if a.id == can.ybr:
                #            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can.ybr))
                
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
            
            #if self.params.PRINT_BB_INFO:
            #    print('--> can (after): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))            
            
            if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
                print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (self.nBB,self.nSO,self.nUE,self.LB,self.UB,100*self.gap))
            
            if self.gap <= self.params.BB_tol:
                conv = True
                self.LB = self.UB
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by optimality gap')
                break
            
            if (time.time() - self.t0) >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print('--> time limit exceeded')
                break
            
            self.nBB += 1
 
        self.rt = time.time() - self.t0
 
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,self.rt,self.nBB,self.nSO,self.nUE,self.UB,100*self.gap))
            print(self.rt_TAP)
            print(self.rt_MILP)
            print(self.yopt)

        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')