import time
import copy
from src import BB_node
from src import Params
from docplex.mp.model import Model

class FS_NETS:

    # instantiate Leblanc's branch-and-bound algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.inf = 1e+9
        self.tol = 1e-2        
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.params = Params.Params()
        
        self.yvec = {}
        self.cuts = {}
        
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
                             
        fixed00 = copy.deepcopy(can.fixed0)
        fixed00.append(can.ybr)
        fixed01 = copy.deepcopy(can.fixed1)
        fixed10 = copy.deepcopy(can.fixed0)
        fixed11 = copy.deepcopy(can.fixed1)
        fixed11.append(can.ybr)
        
        cnt = len(self.BB_nodes) 
        
        BB_node_id = cnt+1
        can.children.append(BB_node_id)
        n0 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, self.inf, fixed00, fixed01, False)
        self.BB_nodes.append(n0)
        
        BB_node_id = cnt+2
        can.children.append(BB_node_id)
        n1 = BB_node.BB_node(self.network, BB_node_id, can.id, can.LB, can.UB, fixed10, fixed11, True)
        n1.score = can.score
        self.BB_nodes.append(n1)
    
        return
    
    def getOAcut(self):
        cut = {'a':{}, 'b':{}}
        
        for a in self.network.links:
        
            if a.y == 1:
                
                #cut['a'][(a.start.id,a.end.id)] = a.x * a.getDerivativeTravelTime(a.x) + a.getTravelTime(a.x,'UE')
                cut['a'][(a.start.id,a.end.id)] = a.getTravelTime(a.x,'SO')
                cut['b'][(a.start.id,a.end.id)] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)                
                
                #if self.params.PRINT_BB_INFO:
                #    print('%s\t\t%.1f\t\t%.3f\t\t%.1f\t\t%.1f\t\t%.1f' 
                #          % ((a.start.id,a.end.id),a.x,a.getDerivativeTravelTime(a.x),a.getTravelTime(a.x,'SO'),cut['a'][(a.start.id,a.end.id)],cut['b'][(a.start.id,a.end.id)]))
            
            else:
                cut['a'][(a.start.id,a.end.id)] = 0
                cut['b'][(a.start.id,a.end.id)] = 0
                
        return cut
    
    def getYvec(self,y):
        #---transforms a dict of link objects into a dict of pairs of ints
        yid = {}
        for a in self.network.links2:
            yid[(a.start.id,a.end.id)] = y[a]
        
        return yid
    
    def milp_link(self,can):
        
        B = self.network.B
        Q = self.network.TD        
        D = [i.id for i in self.network.zones]
        V = [i.id for i in self.network.nodes]
        A = [(a.start.id,a.end.id) for a in self.network.links]        
        A2 = [(a.start.id,a.end.id) for a in self.network.links2]
        g = {(a.start.id,a.end.id):a.cost for a in self.network.links2}
        
        d = {(i,s):0 for i in V for s in D}
        for r in self.network.zones:
            for s in self.network.zones:
                d[r.id,s.id] = r.getDemand(s)  
        for s in D:
            d[s,s] = - sum(d[j,s] for j in D)
        
        l2n = {}
        for a in self.network.links:
            l2n[a.id] = (a.start.id,a.end.id)            
        
        #---setup RMP 
        milp = Model()    
        
        milp.x = {(i,j):milp.continuous_var(lb=0,ub=Q) for (i,j) in A}
        milp.xc = {(i,j,s):milp.continuous_var() for (i,j) in A for s in D}
        milp.y = {(i,j):milp.binary_var() for (i,j) in A2}
        milp.mu = milp.continuous_var(lb=can.LB)
        
        milp.add_constraint(sum(milp.y[i,j]*g[i,j] for (i,j) in A2) <= B)
        
        for (i,j) in A2:        
            milp.add_constraint(milp.x[i,j] <= milp.y[i,j]*Q)
            
        for (i,j) in A:
            milp.add_constraint(sum(milp.xc[i,j,s] for s in D) == milp.x[i,j])           
                
        for i in V:
            for s in D:
                milp.add_constraint(sum(milp.xc[i,j,s] for j in V if (i,j) in A) 
                                    - sum(milp.xc[j,i,s] for j in V if (j,i) in A) == d[i,s])       
                        
        for k in self.cuts:
            #---OA cuts
            milp.add_constraint(milp.mu >= sum(milp.x[a]*self.cuts[k]['a'][a] + self.cuts[k]['b'][a] for a in A))
            
        for k in self.yvec:        
            #---Interdiction Cuts
            milp.add_constraint(sum(milp.y[a] + self.yvec[k][a] - 2*milp.y[a]*self.yvec[k][a] for a in A2) >= 1)
            
        #---Branch cuts    
        for a in can.fixed0:
            milp.add_constraint(milp.y[l2n[a]] == 0)
                               
        for a in can.fixed1:
            milp.add_constraint(milp.y[l2n[a]] == 1)        
        
        milp.minimize(milp.mu)
        
        milp.solve(log_output=False)
        
        if milp.solve_details.status == 'infeasible' or milp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf,{}
        
        else:
            LB = milp.objective_value
            status = milp.solve_details.status
            #rt_MILP = milp.solve_details.time
            #gap_MILP = milp.solve_details.gap
             
            yopt = {} #---this is a dict of link objects
            for a in self.network.links2:
                yopt[a] = milp.y[(a.start.id,a.end.id)].solution_value
                
        return status,LB,yopt

    
    def OA_link(self,can):
        nOA = 0
        UB_OA = self.inf
        LB_OA = 0
        
        t0 = time.time()
        
        conv = False
        while conv == False:
            
            milp_status,milp_obj,milp_y = self.milp_link(can)
            
            if milp_status == 'infeasible':
                if self.params.PRINT_BB_INFO:
                    print('MILP infeasible')
                return nOA,self.inf,{}
            
            LB_OA = milp_obj
            self.yvec[len(self.yvec)] = self.getYvec(milp_y)
            
            tstt = self.network.tapas('SO',milp_y)
            self.cuts[len(self.cuts)] = self.getOAcut()
            
            if tstt < UB_OA:
                if self.params.PRINT_BB_INFO:
                    print('*** update UB in OA_link ***')
                UB_OA = tstt
                yOA = milp_y
            
            gap = (UB_OA-LB_OA)/UB_OA
            if self.params.PRINT_BB_INFO:
                print('-----------------------------------> %d\t%.1f\t%.1f\t%.2f%%' % (nOA,LB_OA,UB_OA,100*gap))
            
            if gap <= self.tol:
                conv = True
                LB_OA = min(LB_OA,UB_OA)
                gap = 0
                if self.params.PRINT_BB_INFO:
                    print('Convergence by optimality in OA_link')                    
                
            if (time.time() - t0) >= self.params.BB_timelimit:
                LB_OA = min(LB_OA,UB_OA)
                if self.params.PRINT_BB_INFO:
                    print('Time limit exceeded in OA_link')
                break
            
            nOA += 1
            
        return nOA,LB_OA,yOA
        

    def BB(self):
        nBB = 0
        nSO = 0
        nUE = 0
        yopt = None
 
        t0 = time.time()
    
        conv = False
        while conv == False:
             
            can = self.nodeSelection(self.getCandidates())
            status = can.check()
     
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
                    nOA, can.LB, yOA = self.OA_link(can)
                    nSO += nOA
                    
                    for a in self.network.links2:
                        can.score[a.id] = a.x * a.getTravelTime(a.x,'SO')                    
                    
                    #---solve UE-TAP at root node to get initial UB
                    if nBB == 0:
                        can.UB = self.network.tapas('UE',yOA)
                        nUE += 1
                        
                        if can.UB < self.UB:            
                            self.UB = can.UB
                            yopt = yOA
                            if self.params.PRINT_BB_INFO:
                                print('*** update UB ***')                        
                    
                 
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by optimality')
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
                        print('*** update UB ***')
                    
                    for n in self.BB_nodes:                    
                        if n.active == True and n.LB >= self.UB:
                            n.active = False      

            if prune == False:
                fixed = can.fixed0 + can.fixed1
                free = [a.id for a in self.network.links2 if a.id not in fixed]            
                free_sorted = sorted(free, key = lambda ele: can.score[ele], reverse = True)
                can.ybr = free_sorted[0]
                
                self.branch(can)
         
            can.active = False
            candidates = self.getCandidates()
               
            if len(candidates) == 0:
                conv = True
                self.LB = self.UB
                gap = 0
                if self.params.PRINT_BB_INFO:
                    print('Convergence by inspection')
                break
                
            else:
                self.LB = self.getLB(candidates)            
                gap = self.getGap()
            
            print('==> %d\t%d\t%d\t%.1f\t%.1f\t%.2f%%' % (nBB,nSO,nUE,self.LB,self.UB,100*gap))
            
            if gap <= self.tol:
                conv = True
                self.LB = self.UB
                gap = 0
                if self.params.PRINT_BB_INFO:
                    print('Convergence by optimality')
                break
            
            if (time.time() - t0) >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print('Time limit exceeded')
                break
            
            nBB += 1
 
        rt = time.time() - t0
 
        print(conv,rt,nBB,gap,self.UB,yopt)
        
        return
    