import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class OA_CNDP:
    
    def __init__(self, network):
        self.network = network
        self.g = {a:0.01 for a in self.network.links2}
        
    
    def solve(self):
       
        timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e15
        lb = 0
        gap = 1
        cutoff = 0.01
        
        last_xhat = None
        last_yhat = None
        
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB
            yhat, obj_l = self.solveRMP()
            lb = obj_l
            # solve TAP -> x, UB
            xhat, obj_f = self.TAP(yhat)
            ub = min(ub, obj_f)
            # add VF cut
            self.addVFCut(xhat, yhat)
            # add TSTT cut
            self.addTSTTCut(xhat, yhat)
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (ub - lb)/lb
            else:
                gap = 1
            
            print(iteration, lb, ub, gap, elapsed)
            print("\t", yhat)
            
            if elapsed > timelimit:
                break
                

                
            
    def addVFCut(self, xhat, yhat):
        yhat_ext = dict()
        
        for a in self.network.links:
            if a in self.network.links2:
                yhat_ext[a] = yhat[a]
            else:
                yhat_ext[a] = 0
                
        self.rmp.add_constraint(sum( (self.rmp.x[a] - xhat[a]) * a.getTravelTimeC(xhat[a], yhat_ext[a], "UE") for a in self.network.links) <= 0)
            
    def addTSTTCut(self, xhat, yhat):
        for a in self.network.links:
            if a in self.network.links2:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], yhat[a], "UE")
                secondterm = a.getDerivativeTravelTimeCx(xhat[a], yhat[a])
                thirdterm = a.getDerivativeTravelTimeCy(xhat[a], yhat[a])
            
                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]) + thirdterm * (self.rmp.y[a] - yhat[a]))
            else:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], 0, "UE")
                secondterm = a.getDerivativeTravelTimeCx(xhat[a], 0)
                thirdterm = a.getDerivativeTravelTimeCy(xhat[a], 0)
                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]))
        
    
      
    def initRMP(self):   
        self.rmp = Model()
        self.rmp.mu = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.C/2) for a in self.network.links2}
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.xc = {(a,r):self.rmp.continuous_var(lb=0, ub=r.totaldemand) for a in self.network.links for r in self.network.origins}
        
        for a in self.network.links:
            self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for r in self.network.origins) == self.rmp.x[a])
            
        for i in self.network.nodes:                    
            for r in self.network.origins:            

                if i.id == r.id:
                    dem = - sum(r.getDemand(s) for s in self.network.zones)                
                elif isinstance(i, type(r)) == True:
                    dem = r.getDemand(i)
                else:
                    dem = 0

                self.rmp.add_constraint(sum(self.rmp.xc[(a,r)] for a in i.incoming) - sum(self.rmp.xc[(a,r)] for a in i.outgoing) == dem)
                
        self.rmp.minimize(sum(self.rmp.mu[a] for a in self.network.links) + sum(self.g[a] * self.rmp.y[a] for a in self.network.links2))
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        yhat = {a:self.rmp.y[a].solution_value for a in self.network.links2}

        obj_l = self.rmp.objective_value
        
        return yhat, obj_l
        
    def TAP(self, y):
    
        for a in self.network.links2:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.network.getTSTT("UE")
        
        return xhat, obj_f
 