import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model
import polytope as pc
import numpy as np


class HY_CNDP:
    
    def __init__(self, network):
        self.network = network
        
        
        self.g = {a:a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        
        for a in self.network.links:
            if a.cost > 0:
                self.varlinks.append(a)
                
        
    
    def solve(self):
        time1 = time.time()
        S0, A, b = self.constructS0()
        self.solveSLP(S0, A, b, 1)
        time1 = time.time()
        
        print("time", time1)
    
    
    def solveSLP(self, S0, A, b, muk):
        Sk = S0
        
        # step 1 choose tolerance factor delta k > 0
        delta_k = 1
        
        # let polytope S_o as defined in (6)
        # calculate penalty term H with all-or-nothing flow assignment
        while True:
            extreme = pc.extreme(Sk)

            best_var = None
            best_obj = 1e100
            best_tdotv = 0
            best_t1 = 0

            for i in range(0, len(extreme)):
                SLP_obj, tdotv, t1 = self.evaluate(extreme[i])

                if SLP_obj < best_obj:
                    best_var = extreme[i]
                    best_obj = SLP_obj
                    best_todtv = tdotv
                    best_t1 = t1

            # if tdotv - t1 <= delta^k, terminate with z^k
            if best_tdotv - best_t1 <= delta_k:
                return best_obj
                
            # compute multicutting plane
            # add linear constraints to 
            
    def evaluate(self, vars):
        for r in range(0, len(self.network.origins)):
            for i in range(0, len(self.network.links)):
                index = r * len(self.network.links) + i
                self.network.links[i].x += vars[index]
        
        for i in range(0, len(self.varlinks)):
            self.varlinks[i].add_cap = vars[i + len(self.network.links)]
        
        t1 = vars[-1]
        
        tdotv = sum(a.getTravelTimeC(a.x, a.add_cap) * a.x for a in self.network.links)
        
        output = tdotv
        return output, tdotv, t1
        
    def constructS0(self):
        
        
        A_demo = np.array([[1.0, 0.0],
              [0.0, 1.0],
              [-1.0, -0.0],
              [-0.0, -1.0]])

        b_demo = np.array([2.0, 1.0, 0.0, 0.0])
        
        xvarslen = len(self.network.links)*len(self.network.origins) # x per origin
        xconslen = 2*len(self.network.origins)*len(self.network.nodes) 
        num_cons = xconslen+len(self.varlinks)+1 # +1 for t_1
        cols = xvarslen + len(self.varlinks) + 1 # 1 for t1 limit
        #num_cons = len(self.varlinks)+1
        #cols = len(self.varlinks)+1
        
        
        rows = num_cons + cols
        
        # v first
        
        A = np.zeros((rows, cols))
        b = np.zeros(rows)
        
        print(rows, cols, len(A), len(A[0]))
        
        
        
        # link flows
        for ri in range(0, len(self.network.origins)):
            r = self.network.origins[ri]

            for j in self.network.nodes:
            
                rowi = 2* (ri*len(self.network.nodes) + (j.id-1))
                
                row1 = A[rowi]
                row2 = A[rowi+1]
            
                dem = 0
                if j.id == r.id:
                    dem = - sum(r.getDemand(s) for s in self.network.zones)                
                elif isinstance(j, type(r)) == True:
                    dem = r.getDemand(j)

                for ij in j.incoming:
                    varidx = ri*len(self.network.links)+ij.id
                    row1[varidx] = 1
                    #row2[varidx] = -1
                for jk in j.outgoing:
                    varidx = ri*len(self.network.links)+jk.id 
                    row1[varidx] = -1
                    #row2[varidx] = 1

                b[rowi] = dem
                #b[rowi+1] = -dem
                break
            break
                
        
        # then y
        for i in range(0, len(self.varlinks)):
            rowi = i + xconslen 
            coli = i + xvarslen
            
            rowy = A[rowi]
            rowy[coli] = 1
            b[rowi] = self.varlinks[i].C/2
        
        # t_1 <= M
        A[rows-1][cols-1] = 1
        b[cols-1] = 1e20
        
        for i in range (0, cols):
            rowi = num_cons + i
            A[rowi][i] = -1
            
        

        # Ax <= b
        p = pc.Polytope(A, b)
        
        
        
        #print(len(pc.extreme(p)))
        
        return p, A, b