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
        time1 = time.time() - time1
        
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
                
            print("best", best_var, best_obj, SLP_obj)
            
            # if tdotv - t1 <= delta^k, terminate with z^k
            if best_tdotv - best_t1 <= delta_k:
                return best_obj
                
            # compute multicutting plane
            # add linear constraints to 
            
    def evaluate(self, vars):
        for a in self.network.links:
            a.x = 0
            
        for r in range(0, len(self.network.origins)):
            for i in range(0, len(self.network.links)):
                index = r * len(self.network.links) + i
                self.network.links[i].x += vars[index]
        
        for i in range(0, len(self.varlinks)):
            self.varlinks[i].add_cap = vars[i + len(self.network.links)]
        
        t1 = vars[-1]
        
        tdotv = sum(a.getTravelTimeC(a.x, a.add_cap, "UE") * a.x for a in self.network.links)
        
        output = tdotv
        return output, tdotv, t1
        
    def constructS0(self):
        
        
        A_demo = np.array([[1.0, 0.0],
              [0.0, 1.0],
              [-1.0, -0.0],
              [-0.0, -1.0]])

        b_demo = np.array([2.0, 1.0, 0.0, 0.0])
        
        xvarslen = len(self.network.links)*len(self.network.origins) # x per origin
        print(len(self.network.links), len(self.network.origins))
        
        xconslen = len(self.network.origins)*len(self.network.nodes) 
        num_cons = xconslen+len(self.varlinks)+2  # +1 for t_1
        cols = xvarslen + len(self.varlinks) + 1 # 1 for t1 limit
        
        
        rows = num_cons + cols + xvarslen
        
        # v first
        
        A = np.zeros((rows, cols))
        b = np.zeros(rows)
        
        print(rows, cols, len(A), len(A[0]))
        
        
        rowidx = 0
        
        
        # link flows
        eps = 1e-4
        
        count = 0
        
        for ri in range(0, len(self.network.origins)):
            r = self.network.origins[ri]

            for j in self.network.nodes:
            
                
                
                row1 = A[rowidx]
                #row2 = A[rowidx+1]
                
                
            
                dem = 0
                if j.id == r.id:
                    dem = - sum(r.getDemand(s) for s in self.network.zones)                
                elif isinstance(j, type(r)) == True:
                    dem = r.getDemand(j)

                

                    
                #print(rowidx, "node", j, "origin", r.id)
                
                if dem <= 0:
                    for ij in j.incoming:
                        #print("\tin", ij.id, ij)
                        varidx = ri*len(self.network.links)+ij.id
                        row1[varidx] = 1
                        #row2[varidx] = -1
                    for jk in j.outgoing:
                        #print("\tout", jk.id, jk)
                        varidx = ri*len(self.network.links)+jk.id 
                        row1[varidx] = -1
                        #row2[varidx] = 1
                        
                    b[rowidx] = dem
                    rowidx += 1
                elif dem > 0:    
                    for ij in j.incoming:
                        print("\tin", ij.id, ij)
                        varidx = ri*len(self.network.links)+ij.id
                        row1[varidx] = -1
                        #row2[varidx] = -1
                    for jk in j.outgoing:
                        print("\tout", jk.id, jk)
                        varidx = ri*len(self.network.links)+jk.id 
                        row1[varidx] = 1
                        #row2[varidx] = 1
                    b[rowidx] = -dem
                    rowidx += 1
                
                #b[rowidx+1] = -dem
                
                
                
                
                
         
        
        # then y
        for i in range(0, len(self.varlinks)):
            
            coli = i + xvarslen
            
            rowy = A[rowidx]
            rowy[coli] = 1
            b[rowidx] = self.varlinks[i].C/2
            rowidx += 1
        
        # t_1 <= M
        A[rowidx][cols-1] = 1
        b[rowidx] = 1e20
        
        rowidx += 1
        
        A[rowidx][cols-1] = -1
        b[rowidx] = 0
        
        rowidx += 1
        
        for i in range (0, cols):
            A[rowidx][i] = -1
            rowidx += 1
            
        
        # upper bounds
        
        for i in range(0, len(self.network.origins)):
            for j in range(0, len(self.network.links)):
                
                #print("check ub", num_cons, cols, i, j)
                origin = self.network.origins[i]
                col = i*len(self.network.links)+self.network.links[j].id
                A[rowidx][col] = 1
                b[rowidx] = origin.totaldemand
        
                rowidx += 1
               
        # Ax <= b
        p = pc.Polytope(A, b)
        
        print(A)
        print(b)
        
        #print([4000, 4000, 0, 0, 0, 0] in p)
        
        ext = pc.extreme(p)
        
        #for e in ext:
        #    print(e)
        
        return p, A, b