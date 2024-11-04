import time
from src import BB_node
from src import Params
from src import YDict
from src import Path
from docplex.mp.model import Model
from src import Heap

class CNDP_MILP:
    
    def __init__(self, network, num_x_pieces, num_y_pieces, k_paths, inflate_costs):
        self.network = network
        self.CG_tol = 1e-4
        self.inf = 1e+9
        self.cs_error = 0.01

        self.last_xhat = None
        self.last_obj_f = 1000000000
        
        self.sameycuts = []
        self.sameyvars = []
        
        self.g = {a:a.cost * inflate_costs for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
            
        self.varlinks = []
        
        
        
        for a in self.network.links:
            if a.cost > 0:
                #print(a, self.g[a])
                self.varlinks.append(a)
            
        
        print("yvars", len(self.varlinks))
        
        self.paths = {r:{s:[] for s in self.network.zones if r.getDemand(s) > 0} for r in self.network.origins}
        self.link_cons = dict()
        self.dem_cons = dict()

        
        self.xbounds = dict()
        self.ybounds = dict()
        
        self.rt_pricing = 0.0
        self.params = Params.Params()
        self.pathcount = 0
        
        self.num_x_pieces = num_x_pieces
        self.num_y_pieces = num_y_pieces
        self.k_paths = k_paths
        
    def solve(self):
        
        
        # initialize paths
        self.createPaths(self.k_paths)
        
        if self.params.PRINT_BB_INFO:
            print("init RMP")
        
        t1 = time.time()
        
        self.initRMP()
        
        if self.params.PRINT_BB_INFO:
            print("solving RMP")
            
        RMP_status, ofv, x, y = self.solveRMP()
        
        
        
        
        
       
        
        
        print(RMP_status, ofv)
        
        x_f, obj = self.TAP(y)
        
        
        t1 = time.time()-t1
        
        '''
        for a in self.network.links:
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = y[a]
            print(a, x[a], x_f[a])
            print("\t", y_ext, self.g[a]*y_ext)
            print("\t", self.rmp.tt[a].solution_value, a.getTravelTime(a.x, "UE"))
        
        print(obj, t1)
        '''
        
        print(t1)
        
        
        
        return obj, t1
    
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
            
    
    def createPaths(self, k):
  
        if self.params.PRINT_BB_INFO:
            print("starting path enumeration")
        
        # solve UE to have  path costs corresponding to feasible link flows
        self.network.tapas("UE", None)
        
        t1 = time.time()
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
        
                    self.printAllPaths(r, s)
                    
                    self.paths[r][s].sort(key=lambda path: path.cost)
                    ln = len(self.paths[r][s])
                    self.paths[r][s] = self.paths[r][s][0:min(k, ln)]
                   
        
        t1 = time.time() - t1
        
        if self.params.PRINT_BB_INFO:
            print("finished path enumeration", len(self.getPaths()), t1)
   
    
    

    
    def initRMP(self):   
    
        # init paths
        
        
        
                    
        epsilon = 1e-4
                
        self.rmp = Model()
        
        
        self.rmp.parameters.read.scale = -1
        #self.rmp.parameters.timelimit = self.params.BB_timelimit
        
        self.rmp.y = {a:self.rmp.continuous_var(lb=0, ub=a.max_add_cap) for a in self.varlinks}
       
        #for a in self.varlinks:
        #    self.rmp.add_constraint(self.rmp.y[a] == 0)
            
        self.rmp.x = {a:self.rmp.continuous_var(lb=0, ub=self.network.TD) for a in self.network.links}
        self.rmp.h = dict()
        self.rmp.c = dict()
        self.rmp.sigma = dict()
        self.rmp.zeta = {a: {i : self.rmp.binary_var() for i in range(0, self.num_x_pieces)} for a in self.network.links}
        self.rmp.kappa = {a: {i: self.rmp.binary_var() for i in range(0, self.num_y_pieces)} for a in self.varlinks}
        
        totaldemand = self.network.getTotalTrips()
        
        self.rmp.tt = {a: self.rmp.continuous_var(lb=0) for a in self.network.links}
        
        # variable for min travel time from r to s
        self.rmp.mu = {(r,s): self.rmp.continuous_var(lb=0) for r in self.network.origins for s in self.paths[r].keys()}
        
        
        
        
        for r in self.network.origins:
            for s in self.paths[r].keys():
                for p in self.paths[r][s]:
                    
                    
                    self.rmp.sigma[p] = self.rmp.binary_var() # CS binary variable
                    self.rmp.c[p] = self.rmp.continuous_var(lb=0)
                    self.rmp.h[p] = self.rmp.continuous_var(lb=0, ub=r.getDemand(s))
                    
                    
                    self.rmp.add_constraint(self.rmp.c[p] == sum(self.rmp.tt[a] for a in p.links))
                    self.rmp.add_constraint(self.rmp.mu[(r,s)] <= self.rmp.c[p])
                    
                    
                    # CS
                    M = sum(a.getTravelTimeC(totaldemand, 0, "UE") for a in p.links)
                    
                    self.rmp.add_constraint(epsilon * (1-self.rmp.sigma[p]) <= self.rmp.c[p] - self.rmp.mu[(r,s)])
                    
                    
                    self.rmp.add_constraint(self.rmp.h[p] <= r.getDemand(s) * self.rmp.sigma[p])
                    
                    self.rmp.add_constraint(self.rmp.mu[(r,s)] + epsilon >= self.rmp.c[p] - M *(1 - self.rmp.sigma[p]))
                    
                    #phi = 1 if mu <= c+epsilon
                    
                
                # VI at least one path is used per OD
                
                self.rmp.add_constraint(sum(self.rmp.sigma[p] for p in self.paths[r][s]) >= 1)
        
        
        for a in self.network.links:
            self.link_cons[a] = self.rmp.add_constraint(self.rmp.x[a] - sum(self.rmp.h[p] for p in self.getPaths() if a in p.links) == 0)
            
        for r in self.network.origins:
            for s in self.network.zones:
                if r.getDemand(s) > 0:
                    self.dem_cons[(r,s)] = self.rmp.add_constraint(sum(self.rmp.h[p] for p in self.paths[r][s]) == r.getDemand(s))
        
        
        
        # number of x bounds should be x_pieces+1
        
        
        for a in self.network.links:
            self.xbounds[a] = list()
            
            # last piece is at network total demand
            
            if totaldemand <= a.C*2:
                skip = totaldemand / (self.num_x_pieces)
                
                for i in range(0, self.num_x_pieces+1):
                    self.xbounds[a].append(i*skip)
            else:
                skip = a.C*2 / (self.num_x_pieces-1)

                for i in range(0, self.num_x_pieces):
                    self.xbounds[a].append(i*skip)

                self.xbounds[a].append(totaldemand)
            
            
            for i in range(0, self.num_x_pieces):
                
                self.rmp.add_constraint(self.rmp.x[a] - self.xbounds[a][i+1] <=  totaldemand * (1 - self.rmp.zeta[a][i]) + epsilon)
                self.rmp.add_constraint(self.xbounds[a][i] * self.rmp.zeta[a][i] - epsilon <= self.rmp.x[a])
            
            self.rmp.add_constraint(sum(self.rmp.zeta[a][i] for i in range(0, self.num_x_pieces)) == 1)
            
            
              
        for a in self.varlinks:
            self.ybounds[a] = list()
            
            skip = a.max_add_cap / self.num_y_pieces
            
            for i in range(0, self.num_y_pieces+1):
                self.ybounds[a].append(i*skip)
               
            for i in range(0, self.num_y_pieces): 
                self.rmp.add_constraint(self.rmp.y[a] - self.ybounds[a][i+1] <= a.max_add_cap * (1 - self.rmp.kappa[a][i]) + epsilon)
                self.rmp.add_constraint(self.ybounds[a][i] * self.rmp.kappa[a][i] - epsilon <= self.rmp.y[a])
        
        
            self.rmp.add_constraint(sum(self.rmp.kappa[a][i] for i in range(0, self.num_y_pieces)) == 1)
        
        
        
        
        
        # linear approximation of link travel times
        for a in self.network.links:
            
            
                
                
            if a in self.varlinks:
                for i in range(0, self.num_x_pieces):
                    for j in range(0, self.num_y_pieces):
                        min_tt = a.getTravelTimeC(self.xbounds[a][i], self.ybounds[a][j], "UE")
                        partial_tt = a.getTravelTimeC(self.xbounds[a][i+1], self.ybounds[a][j], "UE")
                        max_tt = a.getTravelTimeC(self.xbounds[a][i+1], self.ybounds[a][j+1], "UE")

                        
                        pa, pb, pc = self.equation_plane(self.xbounds[a][i], self.ybounds[a][j], min_tt, self.xbounds[a][i+1], self.ybounds[a][j], partial_tt, self.xbounds[a][i+1], self.ybounds[a][j+1], max_tt)
                        # Z = ax + by + c
                        
                        linear_approx = pa * self.rmp.x[a] + pb * self.rmp.y[a] + pc
                        L = - min_tt
                        U = pa*totaldemand + pc
                        L = -10000
                        U = 10000
                        self.rmp.add_constraint(self.rmp.tt[a] - linear_approx >= L*(2-self.rmp.zeta[a][i] - self.rmp.kappa[a][j]))
                        self.rmp.add_constraint(self.rmp.tt[a] - linear_approx <= U*(2-self.rmp.zeta[a][i] - self.rmp.kappa[a][j]))
            else:
                for i in range(0, self.num_x_pieces):
                    min_tt = a.getTravelTimeC(self.xbounds[a][i], 0, "UE")
                    max_tt = a.getTravelTimeC(self.xbounds[a][i+1], 0, "UE")

                    cons = min_tt
                    slope_x = (max_tt - min_tt) / (self.xbounds[a][i+1] - self.xbounds[a][i])
                    
                    linear_approx = cons + (self.rmp.x[a] - self.xbounds[a][i]) * slope_x
                    L = - min_tt
                    U = cons + slope_x * (1.1*totaldemand)
                    U = 10000
                    #print(a, L, U, totaldemand, cons + slope_x * (1.1*totaldemand))
                    self.rmp.add_constraint(self.rmp.tt[a] - linear_approx+ epsilon >= L*(1-self.rmp.zeta[a][i]))
                    self.rmp.add_constraint(self.rmp.tt[a] - linear_approx - epsilon <= U*(1-self.rmp.zeta[a][i]))
                  
                    
        
        self.rmp.minimize(sum(self.rmp.mu[(r,s)] * r.getDemand(s) for r in self.network.origins for s in self.paths[r].keys()) + sum(self.g[a] * self.rmp.y[a] for a in self.varlinks))
        #self.rmp.minimize(sum(self.rmp.mu[(r,s)] * r.getDemand(s) for r in self.network.origins for s in self.paths[r].keys()))
        
    
        
    def solveRMP(self):
    
        t_solve = time.time()
        self.rmp.solve(log_output=False)
        t_solve = time.time() - t_solve
        

        if self.rmp.solve_details.status == 'infeasible' or self.rmp.solve_details.status == 'integer infeasible':
            return 'infeasible',self.inf, dict(), dict()
        RMP_status = self.rmp.solve_details.status
        
        
        x = {a:self.rmp.x[a].solution_value for a in self.network.links}
        y = {a:self.rmp.y[a].solution_value for a in self.varlinks}
        #print(RMP_status)

        OFV = self.rmp.objective_value
        
        #print("ofv check", (OFV - sum(self.g[a] * self.rmp.y[a].solution_value for a in self.varlinks))/self.network.getTotalTrips())
        
        '''
        for a in self.network.links:
            print("\t", a, x[a], self.rmp.tt[a].solution_value)
            for i in range(0, self.num_x_pieces):
                
                print("\t\t", self.xbounds[a][i], self.xbounds[a][i+1], self.rmp.zeta[a][i].solution_value)
        
        
        for p in self.getPaths():
            print(p.links, self.rmp.sigma[p].solution_value, self.rmp.h[p].solution_value, self.rmp.c[p].solution_value, self.rmp.mu[(p.r, p.s)].solution_value)
        '''
        
        
        return RMP_status, OFV, x, y
        
   
        
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
        
    def equation_plane(self, x1, y1, z1, x2, y2, z2, x3, y3, z3): 
     
        a1 = x2 - x1
        b1 = y2 - y1
        c1 = z2 - z1
        a2 = x3 - x1
        b2 = y3 - y1
        c2 = z3 - z1
        a = b1 * c2 - b2 * c1
        b = a2 * c1 - a1 * c2
        c = a1 * b2 - b1 * a2
        d = (- a * x1 - b * y1 - c * z1)
        
        # ax+by+cz+d=0
        # z = -a/c x - b/c y - d/c
        
        return -a/c, -b/c, -d/c
 
 
    def addPath(self, nodelist):
        output = Path.Path()
        output.r = nodelist[0]
        output.s = nodelist[-1]
        
        for i in range(0, len(nodelist)-1):
            output.links.add(self.network.findLink(nodelist[i], nodelist[i+1]))
            
        #print(output.links)
        output.cost = output.getTravelTime("UE")
        self.paths[output.r][output.s].append(output)
 
    def printAllPathsUtil(self, u, d, visited, path):
 
        # Mark the current node as visited and store in path
        visited[u]= True
        path.append(u)
 
        # If current vertex is same as destination, then print
        # current path[]
        if u == d:
            #print (path)
            self.addPath(path)
        else:
            # If current vertex is not destination
            # Recur for all the vertices adjacent to this vertex
            for ui in u.outgoing:
                i = ui.end
                if visited[i]== False:
                    self.printAllPathsUtil(i, d, visited, path)
                     
        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[u]= False
  
  
    # Prints all paths from 's' to 'd'
    def printAllPaths(self, s, d):
 
        # Mark all the vertices as not visited
        visited ={i:False for i in self.network.nodes}
 
        # Create an array to store paths
        path = []
 
        # Call the recursive helper function to print all paths
        self.printAllPathsUtil(s, d, visited, path)
  

 