import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class OA_CNDP_CG:
    
    def __init__(self, network):
        self.network = network
        self.CG_tol = 1e-4
        self.inf = 1e+9
        
        self.g = {a:a.cost for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = self.network.links
        
        self.paths = {r:{s:[] for s in self.network.zones} for r in self.network.origins}
        self.link_cons = dict()
        self.dem_cons = dict()
        
        self.rt_pricing = 0.0
        self.params = Params.Params()
    
    def solve(self):
       
        timelimit = 3600
        iteration = 0
        starttime = time.time()
        ub = 1e15
        lb = 0
        gap = 1
        cutoff = 0.01
        tap_time = 0
        
        last_xhat = {a:0 for a in self.network.links}
        last_yhat = {a:0 for a in self.varlinks}
        last_x_l = {a:0 for a in self.network.links}
        
        self.initRMP()
        
        while gap > cutoff:
            iteration += 1
            
            # solve RMP -> y, LB
            CG_status, obj_l, x_l, yhat = self.CG()
            
            
            
            
            lb = obj_l
            B_l = self.calcBeckmann(x_l, yhat)
            
            # solve TAP -> x, UB
            t1 = time.time()
            xhat, obj_f = self.TAP(yhat)
            t1 = time.time()-t1
            tap_time += t1
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
            
            print(iteration, lb, ub, gap, elapsed, tap_time, B_f, B_l - B_f)
            
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
    
        # init paths
        
        for r in self.network.origins:
            if r.totaldemand > 0:
                self.network.dijkstras(r, "UE")

                self.paths[r] = dict()
                
                for s in self.network.zones:
                    if r.getDemand(s) > 0:
                        p = self.network.trace(r, s)
                        self.paths[r][s] = list()
                        self.paths[r][s].append(p)
                    
                
        self.rmp = Model()
        self.rmp.mu = {a:self.rmp.continuous_var(lb=0,ub=1e10) for a in self.network.links}
        #self.rmp.eta = {a:self.rmp.continuous_var(lb=-1e10,ub=1e10) for a in self.network.links}
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.C/2) for a in self.varlinks}
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.h = {p:self.rmp.continuous_var(lb=0) for p in self.getPaths()}
        

        for a in self.network.links:
            self.link_cons[a] = self.rmp.add_constraint(self.rmp.x[a] - sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) >= 0)
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    self.dem_cons[(r,s)] = self.rmp.add_constraint(sum(self.rmp.h[p] for p in self.paths[r][s]) >= r.getDemand(s))
        
        self.rmp.minimize(sum(self.rmp.mu[a] for a in self.network.links) + sum(self.g[a] * self.rmp.y[a] for a in self.varlinks))
        
    def pricing(self, link_duals, dem_duals):
        
        t0_pricing = time.time()
        
        new = 0
        minrc = self.inf
        
        for a in self.network.links:
            a.dual = link_duals[a]
        
        for r in self.network.origins:
            self.network.dijkstras(r,'RC')
            
            for s in self.network.zones:
                
                if r.getDemand(s) > 0:
                
                    rc = - dem_duals[(r,s)] + s.cost                    
                    
                    if rc <= - self.CG_tol:
                        p = self.network.trace(r,s)
                        self.paths[r][s].append(p) #---is it needed to store paths if directly adding to RMP?
                        
                        #---add new path var to RMP                                                
                        self.rmp.h[p] = self.rmp.continuous_var(lb=0)
                        
                        #---update RMP constraints
                        self.dem_cons[(r, s)].lhs.add_term(self.rmp.h[p], 1)
                        
                        for a in p.links:
                            self.link_cons[a].lhs.add_term(self.rmp.h[p], -1)
                        
                        new += 1
                    
                    if rc < minrc:
                        minrc = rc
                
        self.rt_pricing += (time.time() - t0_pricing)        
        return minrc
        
    def CG(self):
        conv = False
        nCG = 0
        OFV = self.inf
        
        while conv == False:
            RMP_status, OFV, link_duals, dem_duals = self.solveRMP()
            minrc = self.pricing(link_duals, dem_duals)
            
            if RMP_status == 'infeasible':
                CG_status = 'infeasible'
                return CG_status, self.inf, dict(), dict()
            else:
                CG_status = "solved"

            if self.params.PRINT_BB_INFO:
                npaths = len(self.getPaths())
                print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))
                
            if minrc >= -self.CG_tol:
                conv = True
                
            nCG += 1
        
        yhat = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        x_l = {a:self.rmp.x[a].solution_value for a in self.network.links}
 
        
        return CG_status, OFV, x_l, yhat
        
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        
        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf, dict(), dict()
        RMP_status = self.rmp.solve_details.status
        
        OFV = self.rmp.objective_value
        link_duals = {a: self.link_cons[a].dual_value for a in self.network.links}
        
        dem_duals = dict()
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    dem_duals[(r,s)] = self.dem_cons[(r,s)].dual_value
        
        return RMP_status, OFV, link_duals, dem_duals
        
    def TAP(self, y):
    
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV()
        
        return xhat, obj_f
        
        
    def getPaths(self):
        all_paths = []
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    all_paths.append(p)
        return all_paths     
 