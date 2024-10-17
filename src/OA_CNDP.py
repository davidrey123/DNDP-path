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
        last_x_l = None
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB
            x_l, yhat, obj_l = self.solveRMP()
            lb = obj_l
            B_l = self.calcBeckmann(x_l, yhat)
            
            # solve TAP -> x, UB
            xhat, obj_f = self.TAP(yhat)
            B_f = self.calcBeckmann(xhat, yhat)
            ub = min(ub, obj_f)
            # add VF cut
            self.addVFCut(x_l, xhat, yhat)
            # add TSTT cut
            self.addTSTTCut(xhat, yhat)
            
            OA_l = 0
            if last_x_l is not None:
                OA_l = self.calcOABeckmann(x_l, yhat, last_x_l, last_xhat, last_yhat)
            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (ub - lb)/lb
            else:
                gap = 1
            
            print(iteration, lb, ub, gap, elapsed, B_l-B_f, OA_l)
            print("\t", yhat)
            
            #for a in self.network.links:
            #    print("\t", a, x_l[a], xhat[a])
            
            if elapsed > timelimit:
                break
                
            last_xhat = xhat
            last_yhat = yhat
            last_x_l = x_l
                

    def calcBeckmann(self, xhat, yhat):
        total = 0
        
        for a in self.network.links:
            if a in self.network.links2:
                total += a.getPrimitiveTravelTimeC(xhat[a], yhat[a])
            else:
                total += a.getPrimitiveTravelTimeC(xhat[a], 0)
        return total
            
    def addVFCut(self, x_l, xhat, yhat):
        yhat_ext = dict()
        
        B1 = self.calcBeckmann(x_l, yhat)
        B2 = self.calcBeckmann(xhat, yhat)
        
        for a in self.network.links:
            if a in self.network.links2:
                yhat_ext[a] = yhat[a]
            else:
                yhat_ext[a] = 0
                
        self.rmp.add_constraint(B1-B2 + sum( (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], yhat_ext[a], "UE") for a in self.network.links) + sum( (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a])) for a in self.network.links2) <= 0)
    
    def calcOABeckmann(self, x, y, x_l, xhat, yhat):
        total = self.calcBeckmann(x_l, yhat) - self.calcBeckmann(xhat, yhat)
        for a in self.network.links:
            y_rel = 0
            if a in self.network.links2:
                y_rel = yhat[a]
            total += (x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], y_rel, "UE")
            
        for a in self.network.links2:
            total += (y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a]))
            
        return total
      
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
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        obj_l = self.rmp.objective_value
        
        return x_l, yhat, obj_l
        
    def TAP(self, y):
    
        for a in self.network.links2:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.network.getTSTT("UE")
        
        return xhat, obj_f
 