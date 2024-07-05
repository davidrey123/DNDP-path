import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class BPC_nestedTree:

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
        self.UB_OABPC = self.inf #---UB that is used in OABPC after outer root node is solved
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
            n0.y = can.y
            n0.score = can.score
        else:
            n1.y = can.y
            n1.score = can.score
    
        return    
    
    def branch_OABPC(self, can, OABPC_nodes):

        fixed00 = list(can.fixed0)
        fixed00.append(can.ybr)
        fixed01 = list(can.fixed1)
        fixed10 = list(can.fixed0)
        fixed11 = list(can.fixed1)
        fixed11.append(can.ybr)        
        
        cnt = len(OABPC_nodes) 
        
        OABPC_node_id = cnt
        can.children.append(OABPC_node_id)
        n0 = BB_node.BB_node(self.network, OABPC_node_id, can.id, can.LB, self.inf, fixed00, fixed01, False)
        OABPC_nodes.append(n0)
        
        OABPC_node_id = cnt+1
        can.children.append(OABPC_node_id)
        n1 = BB_node.BB_node(self.network, OABPC_node_id, can.id, can.LB, self.inf, fixed10, fixed11, False)
        OABPC_nodes.append(n1)
    
        return OABPC_nodes      
    
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
                yKNP[a] = round(knp.y[a].solution_value)
                
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
        
        for n in range(nKNP):
            yKNP = self.knapsack('x',can)
                        
            t0_TAP = time.time()
            tstt = round(self.network.tapas('SO',yKNP), 3) 
            self.ydict.insertSO(yKNP, tstt)
            self.rt_TAP += (time.time() - t0_TAP)
            self.OAcuts.append(self.getOAcut())                            
            self.nSO += 1        
                    
            self.yvec.append(yKNP)
            
            if tstt < self.UB_SO_DNDP:
                self.UB_SO_DNDP = tstt
                self.yopt = yKNP
                
                for a in self.network.links2:
                    can.score[a.id] = a.x #---update to x*t(x)
                
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
                    rmp.add_constraint(sum(rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))
    
        for a in self.network.links:
            rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in self.getPaths() if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))
                        
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
        rmp.parameters.timelimit = self.params.BB_timelimit
        
        rmp.solve(log_output=False)
        
        #if self.params.PRINT_BB_INFO:
        #    print('nb of paths: %d, cplex time: %.1f, rmp time: %.1f' % (len(can.getPaths()),rmp.solve_details.time,(time.time() - t0_RMP)))
        
        self.rt_RMP += (time.time() - t0_RMP)
            
        if rmp.solve_details.status == 'infeasible' or rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            OFV = rmp.objective_value
            RMP_status = rmp.solve_details.status
            
            yopt = {}
            for a in self.network.links2:
                yopt[a] = rmp.y[a].solution_value
                can.score[a.id] = rmp.x[a].solution_value
            
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
        
    def OABPC(self, can):
                
       #---single search tree OA-BPC to solve SO-DNDP

       n = BB_node.BB_node(self.network, 0, 0, can.LB, self.UB_OABPC, can.fixed0, can.fixed1, False)
       OABPC_nodes = [n]
       nOABPC = 0
       yOABPC = self.yopt #---check
       
       t0_OABPC = time.time()
        
       conv_OABPC = False
       while conv_OABPC == False:                    
            
           candidates_OABPC = [n for n in OABPC_nodes if n.active==True]
           #can_OABPC = self.nodeSelection('depth',candidates_OABPC)
           can_OABPC = self.nodeSelection('bestBound',candidates_OABPC)
           status = can_OABPC.check()
               
           prune = False
           runSO = False        
                
           if status == 'infeasible':
               if self.params.PRINT_BB_INFO:
                   print('----> prune by feasibility')
               prune = True                 
               can_OABPC.LB = self.inf
               can_OABPC.UB = self.inf        
                
           elif status == 'fixed' or status == 'stop':
               if self.params.PRINT_BB_INFO:
                   print('----> prune by check',status)
               prune = True
               #runSO = True #---waste of time if too many OA cuts?
                
               for a in self.network.links2:
                   if a.id in can_OABPC.fixed1:
                       can_OABPC.y[a] = 1
                   elif a.id in can_OABPC.fixed0:
                       can_OABPC.y[a] = 0 
                   else:
                       can_OABPC.y[a] = 0
                       can_OABPC.fixed0.append(a.id)
                
           else:
               #---LB is obtained from LP relaxation of OA MP
               CG_status,can_OABPC.LB,yRMP = self.CG(can_OABPC)
               
               if CG_status == 'infeasible':
                   prune = True
                   can_OABPC.LB = self.inf
                   can_OABPC.UB = self.inf
                   
               elif can_OABPC.LB >= self.UB_OABPC:
                   if self.params.PRINT_BB_INFO:
                       print('----> prune by bounding')
                   prune = True                   
                   
               elif CG_status == 'integral':
                   if self.params.PRINT_BB_INFO:
                       print('----> prune by integrality')

                   prune = True
                   runSO = True

                   for a in self.network.links2:
                       can_OABPC.y[a] = round(yRMP[a])
                                   
           if runSO:
               
               #---add interdiction cut: global since single search tree 
               self.yvec.append(can_OABPC.y)
                
               #---solve SO-TAP to get OA cut
               if self.ydict.hasSO(can_OABPC.y) == True:
                   can_OABPC.UB = self.ydict.getSO(can_OABPC.y)
                   #print('----> has SO')                    
                
               else:
                   t0_TAP = time.time()
                   can_OABPC.UB = round(self.network.tapas('SO',can_OABPC.y), 3)                
                   self.ydict.insertSO(can_OABPC.y, can_OABPC.UB)
                   self.rt_TAP += time.time() - t0_TAP
                   self.OAcuts.append(self.getOAcut())
                   self.nSO += 1

               #---update UB_OABPC                
               if can_OABPC.UB < self.UB_OABPC:
                   self.UB_OABPC = can_OABPC.UB
                   yOABPC = can_OABPC.y
                   
                   if self.params.PRINT_BB_INFO:
                       print('----> update UB OABPC')
                   
                   for n in OABPC_nodes:                    
                       if n.active == True and n.LB >= self.UB_OABPC:
                           n.active = False

           if prune == False:
               
               frac = [a.id for a in can_OABPC.frac]                
               free_sorted = sorted(frac, key = lambda ele: can_OABPC.score[ele], reverse = True)
               can_OABPC.ybr = free_sorted[0]
               
               #if self.params.PRINT_BB_INFO:                    
               #    for a in self.network.links2:
               #        if a.id == can_OABPC.ybr:
               #            print('--> branch on link %s (id: %d)' % ((a.start.id, a.end.id), can_OABPC.ybr))
               
               OABPC_nodes = self.branch_OABPC(can_OABPC,OABPC_nodes)
        
           can_OABPC.active = False
           candidates_OABPC = [n for n in OABPC_nodes if n.active==True]
              
           if len(candidates_OABPC) == 0:
               conv_OABPC = True
               LB_OABPC = self.UB_OABPC 
               gap_OABPC = 0.0
               OABPC_status = 'optimal'
               if self.params.PRINT_BB_INFO:
                   print('----> convergence by inspection')
               
           else:
               LB_OABPC = self.getLB(candidates_OABPC)
               gap_OABPC = (self.UB_OABPC - LB_OABPC)/self.UB_OABPC                      
           
           if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
               print('OABPC> %d\t%d\t%.1f\t%.1f\t%.2f%%' % (nOABPC,self.nSO,LB_OABPC,self.UB_OABPC,100*gap_OABPC))
           
           if gap_OABPC <= self.params.OABPC_tol:
               conv_OABPC = True               
               OABPC_status = 'optimal'
               if self.params.PRINT_BB_INFO:
                   print('----> convergence by optimality gap')
           
           if (time.time() - self.t0) >= self.params.BB_timelimit:
               OABPC_status = 'timelimit'
               if self.params.PRINT_BB_INFO:
                   print('----> time limit exceeded')
               break
           
           nOABPC += 1

       #rt_OABPC = time.time() - t0_OABPC       
       #print('OABPC it: %d, runtime: %.1f, status: %s' % (nOABPC,rt_OABPC,OABPC_status))
       
       return OABPC_status,LB_OABPC,yOABPC

    def BB(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---BPC_nestedTree---')        
        
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
            status = can.check()
        
            if self.params.PRINT_BB_INFO:
                print('--> can (before): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))
     
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
                
                #---after outer root node is solved, used UB_DNDP as UB
                if self.nBB > 0:
                    self.UB_OABPC = self.UB
                
                #---LB is obtained from a single search tree OABPC algorithm
                OABPC_status,can.LB,can.y = self.OABPC(can)
                
                if self.nBB == 0:
                                        
                    self.UB_SO_DNDP = can.LB                    
                    #print('\nSO-DNDP: root node solved: %.1f\n' % self.UB_SO_DNDP)                    
                    runUE = True   
                                        
                if can.LB >= self.UB:                    
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True     
                    
                for a in self.network.links2:
                    can.score[a.id] = round(a.x * a.getTravelTime(a.x,'SO'), 3)                    
                 
            if runUE:
                
                #---solve UE TAP to get UB 
                t0_TAP = time.time()
                can.UB = round(self.network.tapas('UE',can.y), 3)
                self.rt_TAP += time.time() - t0_TAP
                self.nUE += 1                
                
                if can.UB < self.UB:            
                    self.UB = can.UB
                    self.yopt = can.y
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
                
            else:
                self.LB = self.getLB(candidates)
                self.gap = self.getGap()
            
            if self.params.PRINT_BB_INFO:
                print('--> can (after): %d\t%d\t%.1f\t%.1f\t%s\t%s' % (can.id, can.parent, can.LB, can.UB, can.solved, status))            
            
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
            print('---BPC_nestedTree end---')        