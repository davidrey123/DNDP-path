import time
from src import BB_node
from src import Params

class Leblanc:

    # instantiate Leblanc's branch-and-bound algorithm
    def __init__(self, network):
        self.network = network
        self.BB_nodes = []
        self.inf = 1e+9
        self.nit = 0
        self.LB = 0
        self.UB = self.inf
        self.yopt = None
        self.params = Params.Params()
        
        self.nBB = 0
        self.nSO = 0
        self.nUE = 0
        self.rt_TAP = 0.0
        
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

    def BB(self):
        
        print('---Leblanc---')    
        
        self.network.resetTapas()
 
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
                       
                    #---LB is obtained from SO-TAP with unfixed links opened
                    t0_TAP = time.time()
                    can.LB = round(self.network.tapas('SO',y), 3)
                    self.rt_TAP += time.time() - t0_TAP
                    self.nSO += 1
                    
                    for a in self.network.links2:
                        can.score[a.id] = round(a.x * a.getTravelTime(a.x,'SO'), 3)
                        
                        #if self.params.PRINT_BB_INFO:
                        #    print('--> y/score: %d\t%s\t%d\t%.1f' % (a.id, (a.start.id,a.end.id), y[a], can.score[a.id]))
                 
                if can.LB >= self.UB:
                    if self.params.PRINT_BB_INFO:
                        print('--> prune by bounding')
                    prune = True
                     
                else:
                    prune = False
     
            if integral == True:
                                
                t0_TAP = time.time()
                can.UB = round(self.network.tapas('UE',yUB), 3)
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
        print(self.yopt)
        print(self.rt_TAP)
        return
    