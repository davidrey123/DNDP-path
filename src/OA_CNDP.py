import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class OA_CNDP:
    
    def __init__(self, network):
        self.network = network
        
        
        self.g = {a:a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = self.network.links
        
    
    def solve(self):
       
        timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e15
        lb = 0
        gap = 1
        cutoff = 0.01
        
        last_xhat = {a:0 for a in self.network.links}
        last_yhat = {a:0 for a in self.varlinks}
        last_x_l = {a:0 for a in self.network.links}
        
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
            

            
            elapsed = time.time() - starttime
            if lb > 0:
                gap = (ub - lb)/lb
            else:
                gap = 1
            
            print(iteration, lb, ub, gap, elapsed, B_f, B_l - B_f)
            
            #for a in self.varlinks:
            #    print("\t", a, yhat[a], a.C/2)
            
 
            if elapsed > timelimit:
                break
                
            last_xhat = xhat
            last_yhat = yhat
            last_x_l = x_l
    
    def calcOFV(self):
        output = self.network.getTSTT("UE")
        
        for a in self.varlinks:
            output += self.g[a] * a.add_cap
        
        return output

    def calcBeckmann(self, xhat, yhat):
        total = 0
        
        for a in self.network.links:
            if a in self.varlinks:
                total += a.getPrimitiveTravelTimeC(xhat[a], yhat[a])
            else:
                total += a.getPrimitiveTravelTimeC(xhat[a], 0)
                
        return total
            
    def addVFCut(self, x_l, xhat, yhat):
        yhat_ext = dict()
        
        B1 = self.calcBeckmann(x_l, yhat)
        B2 = self.calcBeckmann(xhat, yhat)
        
        for a in self.network.links:
            if a in self.varlinks:
                yhat_ext[a] = yhat[a]
            else:
                yhat_ext[a] = 0
                
        self.rmp.add_constraint(B1-B2 + sum( (self.rmp.x[a] - x_l[a]) * a.getTravelTimeC(x_l[a], yhat_ext[a], "UE") for a in self.network.links) + sum( (self.rmp.y[a] - yhat[a]) * (a.intdtdy(x_l[a], yhat[a]) - a.intdtdy(xhat[a], yhat[a])) for a in self.varlinks) <= 0)
    

      
    def addTSTTCut(self, xhat, yhat):
        for a in self.network.links:
            if a in self.varlinks:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], yhat[a], "UE")
                secondterm = a.getTravelTimeC(xhat[a], yhat[a], "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], yhat[a])
                thirdterm = xhat[a] * a.getDerivativeTravelTimeCy(xhat[a], yhat[a])
            
                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]) + thirdterm * (self.rmp.y[a] - yhat[a]))
            else:
                firstterm = xhat[a] * a.getTravelTimeC(xhat[a], 0, "UE")
                secondterm = a.getTravelTimeC(xhat[a], yhat[a], "UE") + xhat[a] * a.getDerivativeTravelTimeCx(xhat[a], yhat[a])

                self.rmp.add_constraint(self.rmp.mu[a] >= firstterm + secondterm * (self.rmp.x[a] - xhat[a]))
        
    
      
    def initRMP(self):   
        self.rmp = Model()
        self.rmp.mu = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.max_add_cap) for a in self.varlinks}
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
        
        
        self.rmp.minimize(sum(self.rmp.mu[a] for a in self.network.links) + sum(self.g[a] * self.rmp.y[a] for a in self.varlinks))
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        yhat = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
        obj_l = self.rmp.objective_value
        
        return x_l, yhat, obj_l
        
    def TAP(self, y):
    
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV()
        
        return xhat, obj_f
 