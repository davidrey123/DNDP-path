import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class FS_NETS:

    #---instantiate Leblanc's branch-and-bound algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.inf = 1e+9
        self.OA_tol = 1e-2
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        
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
        n1 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, can.UB, fixed10, fixed11, True)
        n1.score = can.score
        self.BB_nodes.append(n1)
    
        return
    
    def getOAcut(self):
        cut = {'a':{}, 'b':{}}
        
        for a in self.network.links:
        
            if a.y == 1:                
                cut['a'][a] = a.getTravelTime(a.x,'SO')
                cut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)                
                
                #if self.params.PRINT_BB_INFO:
                #    print('%s\t\t%.1f\t\t%.3f\t\t%.1f\t\t%.1f\t\t%.1f' 
                #          % ((a.start.id,a.end.id),a.x,a.getDerivativeTravelTime(a.x),a.getTravelTime(a.x,'SO'),cut['a'][(a.start.id,a.end.id)],cut['b'][(a.start.id,a.end.id)]))
            
            else:
                cut['a'][a] = 0
                cut['b'][a] = 0
                
        return cut
    
    def milp_link(self,can):
        
        #---setup RMP 
        milp = Model()    
        
        milp.x = {a:milp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        milp.xc = {(a,s):milp.continuous_var() for a in self.network.links for s in self.network.zones}
        milp.y = {a:milp.binary_var() for a in self.network.links2}
        milp.mu = milp.continuous_var(lb=can.LB)
        
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
                        
        for OAcut in self.OAcuts:
            #---OA cuts
            milp.add_constraint(milp.mu >= sum(milp.x[a]*OAcut['a'][a] + OAcut['b'][a] for a in self.network.links))
            
        for yv in can.yvec:        
            #---Interdiction Cuts
            milp.add_constraint(sum(milp.y[a] + yv[a] - 2*milp.y[a]*yv[a] for a in self.network.links2) >= 1)
            
        #---Branch cuts    
        for a in self.network.links2:
            
            if a.id in can.fixed0:
                milp.add_constraint(milp.y[a] == 0)
                               
            if a.id in can.fixed1:
                milp.add_constraint(milp.y[a] == 1)
        
        milp.minimize(milp.mu)
        
        milp.solve(log_output=False)
        
        if milp.solve_details.status == 'infeasible' or milp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            LB = milp.objective_value
            status = milp.solve_details.status
            #rt_MILP = milp.solve_details.time
            #gap_MILP = milp.solve_details.gap
             
            yopt = {}
            for a in self.network.links2:
                yopt[a] = int(milp.y[a].solution_value)
                
        return status,LB,yopt    

    
    def OA_link(self,can):
        nOA = 0
        nSO_OA = 0
        UB_OA = self.inf
        LB_OA = 0
        yOA = {}
        
        t0 = time.time()
        
        conv = False
        while conv == False:
            
            milp_status,milp_obj,milp_y = self.milp_link(can)
            
            if milp_status == 'infeasible':
                if nOA == 0:
                    print('==========WARNING: MILP infeasible at first iteration===================')
                    
                conv = True
                LB_OA = min(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by feasibility (due to interdiction cuts)')
                return nOA,LB_OA,yOA,milp_status
            
            LB_OA = milp_obj
            can.LB = milp_obj
            
            if LB_OA >= self.UB:
                conv = True
                LB_OA = min(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by bounding in OA_link')
                break
            
            can.yvec.append(milp_y)           
            
            if self.ydict.hasSO(milp_y) == True:
                tstt = self.ydict.getSO(milp_y)
                #print('*** hasSO ***', tstt)
            
            else:            
                tstt = self.network.tapas('SO',milp_y)
                self.ydict.insertSO(milp_y, tstt)
                self.OAcuts.append(self.getOAcut())
                nSO_OA += 1
            
            if tstt < UB_OA:
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> update UB in OA_link')
                UB_OA = tstt
                yOA = milp_y
                
                for a in self.network.links2:
                    can.score[a.id] = a.x * a.getTravelTime(a.x,'SO')                
            
            gap = (UB_OA-LB_OA)/UB_OA
            if self.params.PRINT_BB_INFO:
                print('-----------------------------------> %d\t%.1f\t%.1f\t%.2f%%' % (nOA,LB_OA,UB_OA,100*gap))
            
            if gap <= self.OA_tol:
                conv = True
                LB_OA = min(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> convergence by optimality gap in OA_link')                    
                
            if (time.time() - t0) >= self.params.BB_timelimit:
                LB_OA = min(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('-----------------------------------> time limit exceeded in OA_link')
                break
            
            nOA += 1
            
        return nOA,nSO_OA,LB_OA,yOA,milp_status
        

    def BB(self):
        
        print('---FS_NETS---')
        
        nBB = 0
        nSO = 0
        nUE = 0
        yopt = None
        
        self.network.resetTapas()
 
        t0 = time.time()
    
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
                if can.solved == False:
                    #---closed link child
                    y = {}
                    for a in self.network.links2:
                        if a.id in can.fixed1:
                            y[a] = 1
                        elif a.id in can.fixed0:
                            y[a] = 0 
                        else:
                            y[a] = 1                                        
                    
                    #---LB is obtained from OA algorithm
                    nOA, nSO_OA, can.LB, yOA, OA_status = self.OA_link(can)
                    nSO += nSO_OA
                                               
                    #if OA_status != 'infeasible' and len(yOA) > 0:
                    if self.params.PRINT_BB_INFO:
                        print(OA_status)
                        for a in self.network.links2:
                            print('--> yOA/score: %d\t%s\t%d\t%.1f' % (a.id, (a.start.id,a.end.id), yOA[a], can.score[a.id]))
                    
                    #---solve UE-TAP at root node to get initial UB
                    if nBB == 0:
                        can.UB = self.network.tapas('UE',yOA)
                        nUE += 1
                        
                        if can.UB < self.UB:            
                            self.UB = can.UB
                            yopt = yOA
                            if self.params.PRINT_BB_INFO:
                                print('--> update UB')                        
                    
                 
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True
                     
                else:
                    prune = False
     
            if integral == True:
                                
                can.UB = self.network.tapas('UE',yUB)
                nUE += 1
                
                if can.UB < self.UB:            
                    self.UB = can.UB
                    yopt = yUB
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
            
            print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (nBB,nSO,nUE,self.LB,self.UB,100*gap))
            
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
            
            nBB += 1
 
        rt = time.time() - t0
 
        print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,rt,nBB,nSO,nUE,self.UB,100*gap))
        print(yopt)
        
        return
    