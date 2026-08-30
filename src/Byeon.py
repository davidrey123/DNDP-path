import time
from src import BB_node
from src import Params
from src import YDict
from docplex.mp.model import Model

class Byeon:

    #---instantiate a Byeon's algorithm
    def __init__(self, network):
        self.network = network
        self.INT_tol = 1e-4
        self.inf = 1e+9
        self.nit = 0
        self.LB = 0
        self.UBSP = self.inf
        self.UB = self.inf
        self.gap = self.inf
        self.yopt = None
        self.params = Params.Params()
        self.ydict = YDict.YDict()
        self.t0 = 0.0
        self.M = self.inf
        self.rootNodeLB = 0
        self.ytemp = None
        self.Dtemp = None
        
        self.nBD = 0
        self.nOA = 0
        self.nUE = 0
        self.rt = 0.0
        self.rt_TAP = 0.0
        self.rt_MP = 0.0
        self.rt_SP = 0.0
        self.rt_rootNode = 0.0        
        
        self.mp = None
        self.sp = None
        self.cntUnscaledInf = 0 
        self.nInitYvec = len(self.network.links2)
        
        self.BDcuts = []
        
        self.dict_cons = {}
      
    def getGap(self):           
        #return (self.UB - self.LB)/self.UB
        return (self.UBSP - self.LB)/self.UBSP
    
    def getOAcuts(self):
        
        #---determine OA of x_a^{p_a + 1} and add link-based OA cuts to LP

        for a in self.network.links:
            OAcut = {}
            
            if a.x > self.INT_tol:
                
                addOAcut = True
                for cut in a.OAcuts: 
                    
                    #---check if a.x is sufficiently different from existing OA cuts        
                    if abs(a.x - cut['x']) <= self.params.OAcut_tol*a.C:
                        addOAcut = False
                        break
                
                if addOAcut == True:
                    
                    OAcut['x'] = a.x
                    OAcut['a'] = (a.beta + 1) * pow(a.x, a.beta)
                    OAcut['b'] = - (a.beta) * pow(a.x, a.beta + 1)
                                                  
                    a.OAcuts.append(OAcut)
                    
                    self.nOA += 1
    
    def knapsack(self, yvec):
        
        knp = Model()
        knp.y = {a:knp.binary_var() for a in self.network.links2}
        knp.add_constraint(sum(knp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        #---No-good cuts
        for yv in yvec:
            knp.add_constraint(sum(knp.y[a] + yv[a] - 2*knp.y[a]*yv[a] for a in self.network.links2) >= 1)
        
        knp.maximize(sum(knp.y[a] * a.x for a in self.network.links2))
            
        knp.solve(log_output=False)
            
        if knp.solve_details.status == 'infeasible' or knp.solve_details.status == 'integer infeasible':
            print('infeasible instance?')
            return {}
        else:            
            yKNP = {}
            for a in self.network.links2:
                yKNP[a] = round(knp.y[a].solution_value)
                
            return yKNP
        
    def kBestKNP(self):
        
        tsttBest = self.inf
        yBest = None
        yvec = []
        for n in range(self.nInitYvec):
            yKNP = self.knapsack(yvec)
                       
            tstt = self.network.tapas('SO_OA_cuts',yKNP)
            self.ydict.insertSO(yKNP, tstt)
            self.getOAcuts()                            
            
            if tstt < tsttBest:
                tsttBest = tstt
                yBest = yKNP                
                
            #---interdict yKNP to get different y
            yvec.append(yKNP)
        
        return yBest
    
    def LocalSearchKNP(self, nLS):
        
        yvec = []
        for n in range(nLS):
        
            yKNP = self.knapsack(yvec) 
            yvec.append(yKNP)
            
            tstt = self.network.tapas('SO_OA_cuts',yKNP)
            self.ydict.insertSO(yKNP, tstt)
            self.getOAcuts()                            
            
            if n==0:
                yBest = yKNP
                tsttBest = tstt
            
            for a in self.network.links2:
                
                yLS = dict(yKNP)
                
                if yLS[a] == 1:
                    yLS[a] = 0
                    
                else:
                    yLS[a] = 1            
                    
                tstt = self.network.tapas('SO_OA_cuts',yLS)
                self.ydict.insertSO(yLS, tstt)
                self.getOAcuts()
                
                yCost = sum(b.cost*yLS[b] for b in yLS)
                
                if yCost <= self.network.B:
                
                    if tstt < tsttBest:
                        tsttBest = tstt
                        yBest = yLS
        
        return yBest
    
    def LocalSearchY1(self):
        
        yBest = self.knapsack([])                       
        tstt = self.network.tapas('SO_OA_cuts',yBest)
        self.ydict.insertSO(yBest, tstt)
        self.getOAcuts()  
        
        y1 = {a:1 for a in self.network.links2}
        
        for a in self.network.links2:
            
            yLS = dict(y1)            
            yLS[a] = 0
                
            tstt = self.network.tapas('SO_OA_cuts',yLS)
            self.ydict.insertSO(yLS, tstt)       
            self.getOAcuts() 
                    
        return yBest    
        
    def initOABDcuts(self, type):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('Initializing OA & BD cuts with',self.params.initOAheuristic)
            
        t0_OA = time.time()
        
        #---solve SO-TAP(y1) to initialize OA cuts and branching scores
        y1 = {a:1 for a in self.network.links2}
        
        self.network.tapas('SO_OA_cuts',y1)
        self.getOAcuts()
                
        #---select heuristic to explore y vectors and initialize OA cuts
        if type == 'kBestKNP':
            yBest = self.kBestKNP()
            
        elif type == 'LocalSearchKNP':
            yBest = self.LocalSearchKNP(1)
            
        elif type == 'LocalSearchY1':
            yBest = self.LocalSearchY1()            
            
        self.rt_OA = time.time() - t0_OA
                
        #---solve UE-TAP(yBest) to get UB
        t0_TAP = time.time()
        self.UB = self.network.tapas('UE',yBest)
        self.ydict.insertUE(yBest,self.UB)
        beck = self.network.getBeckmannOFV()         
        self.rt_TAP += time.time() - t0_TAP
        self.nUE += 1
        
        #print(yBest)
        #print('TAP TSTT %.1f\t\tBECK %.1f' % (self.UB,beck))

        #---solve SP to get BD cut and update MP
        self.createSP(yBest,beck)        
        status,OFV,UBSP = self.solveSP()        
        self.UBSP = UBSP
        self.updateMP(self.BDcuts[-1],yBest,beck)       
        #print('SP OFV %.1f' % OFV)
        self.yopt = yBest
    
    def createMP(self):
        
        mp = Model()
        mp.z = mp.continuous_var(lb=0)
        mp.y = {a:mp.binary_var() for a in self.network.links2}
        mp.add_constraint(sum(mp.y[a] * a.cost for a in self.network.links2) <= self.network.B)
        
        mp.minimize(mp.z)
        
        self.mp = mp
        
    def updateMP(self,duals,yhat,beck):      
        
        CTE = duals[0]
        psi = duals[1]
        w = duals[2]
        
        #print(CTE,psi,w)
        
        # note: using beck instead of computing follower obj based on SP solution because it should be the same        
        #self.mp.add_constraint(self.mp.z >= CTE + sum(psi[a] * (- self.mp.y[a] * self.network.TD) for a in self.network.links2) - w * (beck + (self.M - beck) * sum(self.mp.y[a] + yhat[a] - 2*self.mp.y[a]*yhat[a] for a in self.network.links2)))
        #self.mp.add_constraint(self.mp.z >= CTE + sum(psi[a] * (self.mp.y[a] * self.network.TD) for a in self.network.links2) - w * (beck + (self.M - beck) * sum(self.mp.y[a] + yhat[a] - 2*self.mp.y[a]*yhat[a] for a in self.network.links2)))        
        #self.mp.add_constraint(self.mp.z >= CTE + sum(psi[a] * (self.mp.y[a] * self.network.TD) for a in self.network.links2) + w * (beck + (self.M - beck) * sum(self.mp.y[a] + yhat[a] - 2*self.mp.y[a]*yhat[a] for a in self.network.links2)))        
        
        self.mp.add_constraint(self.mp.z >= CTE + sum(psi[a] * (-self.mp.y[a] * self.network.TD) for a in self.network.links2) + w * (beck + (self.M - beck) * sum(self.mp.y[a] + yhat[a] - 2*self.mp.y[a]*yhat[a] for a in self.network.links2)))                
        #self.mp.add_constraint(self.mp.z >= CTE + sum(psi[a] * (-self.mp.y[a] * self.network.TD) for a in self.network.links2) + w * (beck + sum((self.M - beck)*(1 - self.mp.y[a]) for a in self.network.links2 if a.x > self.INT_tol)))
           
    def solveMP(self):    

        t0_mp = time.time()  

        if self.params.PRINT_BB_INFO:
            print('MP nvars: %d, ncons: %d' % (self.mp.number_of_variables,self.mp.number_of_constraints))
        
        ###
        #for a in self.network.links2:
            #print(a.start.id,a.end.id,self.yopt[a])
            #self.mp.add_constraint(self.mp.y[a] == self.yopt[a])
        #    self.mp.add_constraint(self.mp.y[a] == 0)
        ###
        
        self.mp.solve(log_output=False)
            
        self.rt_MP += (time.time() - t0_mp)
            
        if self.mp.solve_details.status == 'infeasible' or self.mp.solve_details.status == 'integer infeasible':
            print('MP is infeasible')
            return self.inf,{}
        else:            
            yMP = {}
            for a in self.network.links2:
                yMP[a] = round(self.mp.y[a].solution_value)
                
            return self.mp.objective_value,yMP   
    
    def createSP(self, y, D):
                   
        sp = Model()
        
        sp.x = {a:sp.continuous_var(lb=0,ub=self.network.TD) for a in self.network.links}
        sp.xc = {(a,s):sp.continuous_var(lb=0) for a in self.network.links for s in self.network.zones}
        sp.mu = {a:sp.continuous_var(lb=0) for a in self.network.links}        
            
        for a in self.network.links:
            sp.add_constraint(sum(sp.xc[(a,s)] for s in self.network.zones) - sp.x[a] == 0, 'flow_%d' % (a.id))
           
        self.dict_cons = {}
        for i in self.network.nodes:          
            
            if len(i.outgoing)==0 and len(i.incoming)==0:
                continue
            
            for s in self.network.zones:            
                
                if i.id == s.id:
                    dem = - sum(r.getDemand(s) for r in self.network.zones)                
                elif isinstance(i, type(s)) == True:
                    dem = i.getDemand(s)
                else:
                    dem = 0
                
                self.dict_cons[(i,s)] = sp.add_constraint(sum(sp.xc[(a,s)] for a in i.outgoing) - sum(sp.xc[(a,s)] for a in i.incoming) == dem, 'dem_%d_%d' % (i.id,s.id))
               
        for a in self.network.links2:            
            #sp.add_constraint(sp.x[a] <= y[a] * self.network.TD, 'linking_%d' % (a.id))
            sp.add_constraint(-sp.x[a] >= -y[a] * self.network.TD, 'linking_%d' % (a.id))
            
        for a in self.network.links:
            ind = 0
            for OAcut in a.OAcuts:
                #sp.add_constraint(sp.mu[a] >= sp.x[a]*OAcut['a'] + OAcut['b'], 'oacut_%d_%d' % (a.id,ind))
                #sp.add_constraint(sp.x[a]*OAcut['a'] - sp.mu[a] <= -OAcut['b'], 'oacut_%d_%d' % (a.id,ind))
                sp.add_constraint(sp.mu[a] - sp.x[a]*OAcut['a'] >= OAcut['b'], 'oacut_%d_%d' % (a.id,ind))
                ind += 1
            
        sp.add_constraint(sum((a.t_ff * sp.x[a] + (a.gtilde/(a.beta + 1)) * sp.mu[a]) for a in self.network.links) <= D, 'SD')
        
        sp.minimize(sum((a.t_ff * sp.x[a] + a.gtilde * sp.mu[a]) for a in self.network.links))
        
        sp.parameters.threads = self.params.CPLEX_threads
        sp.parameters.timelimit = self.params.BB_timelimit
        sp.parameters.read.scale = 0
        #sp.context.cplex_parameters.preprocessing.presolve = 0
        
        if self.params.PRINT_BB_INFO:
            print('SP nvars: %d, ncons: %d' % (sp.number_of_variables,sp.number_of_constraints))
        
        self.sp = sp
        self.ytemp = y
        self.Dtemp = D
        
    def solveSP(self):
        
        t0_sp = time.time()
        
        self.sp.solve(log_output=False)

        if self.params.PRINT_BB_INFO:
            print('SP solve status:',self.sp.solve_details.status)     
            print(self.sp.solve_details)
            print(self.sp.objective_value)
        
        if self.sp.solve_details.status == 'optimal with unscaled infeasibilities':
            self.sp.parameters.read.scale = -1
            self.sp.solve(log_output=False)
            self.sp.parameters.read.scale = 0            
            self.cntUnscaledInf += 1
            
        if self.sp.solve_details.status == 'infeasible':
            return 'infeasible',self.inf
        
        else:
            OFV = self.sp.objective_value
            SP_status = self.sp.solve_details.status
            
            CTEE = 0
            dual_OFV = 0
            
            #---set link flows to generate OA cuts
            for a in self.network.links:
                a.x = round(self.sp.x[a].solution_value,self.params.rd)                                 
                
                """                
                if a.x > 0:    
                    diff = 100*(pow(a.x, a.beta + 1) - self.sp.mu[a].solution_value)/pow(a.x, a.beta + 1)
                else:
                    diff = 0.0
                print('a=%d\tx_a^{p_a + 1}=%.1f\tmu_a=%.1f\tdiff=%.1f%%' % (a.id,pow(a.x, a.beta + 1),self.sp.mu[a].solution_value,diff))
                """
                
            UBSP = self.network.getTSTT('UE')
            
            if self.params.PRINT_BB_INFO:
                OALHS = sum((a.t_ff * self.sp.x[a].solution_value + (a.gtilde/(a.beta + 1)) * self.sp.mu[a].solution_value) for a in self.network.links) 
                LHS = self.network.getBeckmannOFV()
                gapLHS = 100*(LHS - OALHS)/LHS
                print("OALHS %.1f \t LHS %.1f \t gapLHS %.2f%%" % (OALHS,LHS,gapLHS))   
                
            #---get duals values for BD optimality cut
            psi = {}
            for a in self.network.links2:
                psi[a] = self.sp.get_constraint_by_name('linking_%d' % (a.id)).dual_value
                dual_OFV += psi[a] * (-self.ytemp[a] * self.network.TD)
                CTEE += psi[a] * (-self.ytemp[a] * self.network.TD)
                #print(psi[a])
            
            w = self.sp.get_constraint_by_name('SD').dual_value
            dual_OFV += w * self.Dtemp            
            CTEE += w * self.Dtemp       
            #print(w)
            

            CTE = 0
            for a in self.network.links:
                ind = 0
                for OAcut in a.OAcuts:
                    CTE += OAcut['b'] * self.sp.get_constraint_by_name('oacut_%d_%d' % (a.id,ind)).dual_value                    
                    dual_OFV += OAcut['b'] * self.sp.get_constraint_by_name('oacut_%d_%d' % (a.id,ind)).dual_value                    
                    
                    #print(self.sp.get_constraint_by_name('oacut_%d_%d' % (a.id,ind)).dual_value)
                    
                    ind += 1   

            for i in self.network.nodes:      

                if len(i.outgoing)==0 and len(i.incoming)==0:
                    continue
                
                for s in self.network.zones:            
                    
                    if i.id == s.id:
                        dem = - sum(r.getDemand(s) for r in self.network.zones)                
                    elif isinstance(i, type(s)) == True:
                        dem = i.getDemand(s)
                    else:
                        dem = 0
                    
                    if type(self.dict_cons[(i,s)]) == None:
                        print(i.id,s.id,"None")
                    CTE += dem * self.sp.get_constraint_by_name('dem_%d_%d' % (i.id,s.id)).dual_value                    
                    dual_OFV += dem * self.sp.get_constraint_by_name('dem_%d_%d' % (i.id,s.id)).dual_value
                    
            CTEE = CTEE + CTE
            #print(OFV,dual_OFV,CTEE)
            
            self.BDcuts.append((CTE,psi,w))              

        self.rt_SP += (time.time() - t0_sp)
        
        if self.params.PRINT_BB_INFO:
            print('nvars: %d, ncons: %d, nOAcuts: %d, cplexTime: %.1f, lpTime: %.1f' % (self.sp.number_of_variables,self.sp.number_of_constraints,self.nOA,self.sp.solve_details.time,(time.time() - t0_sp)))
        
        return SP_status,OFV,UBSP

    def BD(self):
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+'---')
        
        self.network.resetTapas()
 
        self.t0 = time.time()
       
        #---create MP
        self.createMP()
        
        #---set bigM for BD cuts       
        y0 = {a:0 for a in self.network.links2}
        tstt = self.network.tapas('UE',y0)        
        self.nUE += 1         
        self.M = self.network.getBeckmannOFV()         
       
        #---initialize OA & BD cuts
        self.initOABDcuts(self.params.initOAheuristic)
         
        conv = False
        while conv == False:                    
             
            #---solve MP to get y 
            LB,y = self.solveMP()
            self.LB = max(self.LB,LB)
        
            if self.params.PRINT_BB_INFO:
                print('MP LB %.1f' % LB)
                #print('yMP',y)
            
            #---solve TAP to get UB and D            
            t0_TAP = time.time()                       
            UB = self.network.tapas('UE',y)   
            beck = self.network.getBeckmannOFV()   
            self.ydict.insertUE(y,UB)
            self.nUE += 1  
            self.rt_TAP += time.time() - t0_TAP
            
            if self.params.PRINT_BB_INFO:
                print('TAP TSTT %.1f \t BECK %.1f' % (UB,beck))               
            
            #---solve SP to get BD cut and update MP
            self.createSP(y,beck)
            SP_status,OFV,UBSP = self.solveSP()        
            self.updateMP(self.BDcuts[-1],y,beck)  
                                 
            if self.params.PRINT_BB_INFO:            
                print('OFVSP %.1f \t UBSP %.1f \t UB %.1f' % (OFV,UBSP,UB))
                        
            t0_OA = time.time()
            self.getOAcuts()
            self.rt_OA += time.time() - t0_OA
            
            if UBSP < self.UBSP:
                self.UBSP = UBSP
                if self.params.PRINT_BB_INFO:
                    print('--> update UBSP')
             
            if UB < self.UB:            
                self.UB = UB
                self.yopt = y
                if self.params.PRINT_BB_INFO:
                    print('--> update UB')              

            if self.nBD == 0:
                self.rootNodeLB = LB
                self.rt_rootNode = time.time() - self.t0
            
            self.gap = self.getGap()                 
            
            if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
                elapsed_time = time.time() - self.t0
                print('%d\t%d\t%d\t%.1f\t%.1f\t%.2f%%\t%.1f\t%.1f' % (self.nBD,self.nOA,self.nUE,self.LB,self.UBSP,100*self.gap,self.UB,elapsed_time))
            
            if self.gap <= self.params.BB_tol:
                conv = True
                if self.params.PRINT_BB_INFO:
                    print('--> convergence by optimality gap')
            
            if (time.time() - self.t0) >= self.params.BB_timelimit:
                if self.params.PRINT_BB_INFO:
                    print('--> time limit exceeded')
                break
            
            self.nBD += 1
 
        self.rt = time.time() - self.t0

        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('%s\t%.1f\t%d\t%d\t%d\t%.1f\t%.2f%%' % (conv,self.rt,self.nBD,self.nOA,self.nUE,self.UB,100*self.gap))
            print("Time TAP: %.2f" % self.rt_TAP)
            print("Time  OA: %.2f" % self.rt_OA)
            print("Time MP: %.2f" % self.rt_MP)
            print("Time SP: %.2f" % self.rt_SP)
            print(self.yopt)
            print(self.cntUnscaledInf)            
        
        if self.params.PRINT_BB_INFO or self.params.PRINT_BB_BASIC:
            print('---'+self.__class__.__name__+' end---')
        