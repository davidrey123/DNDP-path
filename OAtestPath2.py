#---modules
import time
from src import Params
from src import Network
import Bushify
from docplex.mp.model import Model

def getOAcut(network):
    OAcut = {'a':{}, 'b':{}}

    for a in network.links:
               
        OAcut['a'][a] = a.getTravelTime(a.x,'SO')
        OAcut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
        

    return OAcut
    
def pricing(network):
        
    t0_pricing = time.time()

    new = 0
    minrc = 1e15


    for r in network.origins:
        bush_rc, bush_new = r.bush.pricing(rmp)
        minrc = min(minrc, bush_rc)
        new += bush_new
        

    #if self.params.PRINT_BB_INFO:
    #    print('pricing',new,minrc)
      
    return minrc

def getPaths(network, paths):
        all_paths = []
        for r in network.origins:
            for s in network.zones:
                for p in paths[r][s]:
                    all_paths.append(p)
        return all_paths  

def rmp_path(rmp, network, y_ub):
        
    #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively

    t0_RMP = time.time()
    

    rmp.solve(log_output=False)

    #if self.params.PRINT_BB_INFO:
    #    print('nb of paths: %d, cplex time: %.1f, rmp time: %.1f' % (len(can.getPaths()),rmp.solve_details.time,(time.time() - t0_RMP)))


    if rmp.solve_details.status == 'infeasible' or rmp.solve_details.status == 'integer infeasible':
        return 'infeasible',1e15,{}

    else:
        OFV = rmp.objective_value
        RMP_status = rmp.solve_details.status

        yopt = {}
        for a in network.links2:
            yopt[a] = rmp.y[a].solution_value
            
            
        for a in network.links:
            a.x = rmp.x[a].solution_value


        for a in network.links:                    
            a.dual = max(rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0)                    

        for r in network.origins:
            for s in network.zones:
                if r.getDemand(s) > 0:
                    r.bush.dem_duals[s] = max(rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0)


        return RMP_status,OFV,yopt

def CG(rmp, network, y_ub):

    nCG = 0
    conv = False
    
    
    tot_t_solve = 0
    tot_t_price = 0
    while conv == False:        
        t1 = time.time()
        
        t0 = time.time()
        RMP_status,OFV,yRMP = rmp_path(rmp, network, y_ub)
        t_solve = time.time()-t0
        
        if RMP_status == 'infeasible':
            print("\tINFEASIBLE")
            CG_status = 'infeasible'
            break

        t0 = time.time()
        minrc = pricing(network)
        t_price = time.time()-t0
        
        
        nvars = 0
        for r in network.origins:
            for s in r.bush.paths[r]:
                nvars += len(r.bush.paths[s])
                
        tot_time = time.time() - t1
        print("\t", OFV, minrc, round(tot_time, 2), round(t_solve, 2), round(t_price, 2), nvars)
        
        tot_t_solve += t_solve
        tot_t_price += t_price

        #if self.params.PRINT_BB_INFO:
        #    npaths = len(self.getPaths())
        #    print('CG: %d\t%d\t%.1f\t%.2f' % (nCG,npaths,OFV,minrc))

        if minrc >= -0.0001:
            #if self.params.PRINT_BB_INFO:
            #    print('CG converged')
            conv = True

        nCG += 1

    CG_status = "solved"
    
    if conv == True:            


        #npaths = len(self.getPaths())
        #print('CG: %s\t%d\t%d\t%.1f\t%.2f' % (CG_status,nCG,npaths,OFV,minrc))        

        return CG_status,OFV,yRMP,tot_t_price,tot_t_solve

    else:
        return CG_status,1e15,yRMP,tot_t_price,tot_t_solve


  

    
    
net = 'SiouxFalls'
ins = 'SF_DNDP_20_1'

#net = 'EasternMassachusetts'
#ins = 'EM_DNDP_10_1'

tot_time = time.time()



b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}
network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])




bush_flows = {}

lastObj = 0

y_ub = 0

for a in network.links2:
    a.y = 0
    
    


            
for a in network.links2:
    a.y = y_ub
             
rmp = Model()

rmp.x = {a:rmp.continuous_var(lb=0,ub=network.TD) for a in network.links}
rmp.h = {}



for a in network.links:
    a.link_cons = rmp.add_constraint(rmp.x[a] + 0 >= 0, 'link_%d_%d' % (a.start.id,a.end.id)) # add +0 to make it linear expr

for r in network.origins:
    r.bush = Bushify.Bushify(r, network, rmp)
    
    
rmp.mu = {a:rmp.continuous_var(lb=0) for a in network.links}

rmp.y = {a:rmp.continuous_var(lb=0,ub=y_ub) for a in network.links2}


rmp.add_constraint(sum(rmp.y[a] * a.cost for a in network.links2) <= network.B)


dem_cons = {}




         
#for a in network.links:
#    rmp.add_constraint(rmp.mu[a] >= rmp.x[a]*a.t_ff)
       
rmp.minimize(sum(rmp.mu[a] for a in network.links))


for a in network.links2:
    rmp.add_constraint(rmp.x[a] <= rmp.y[a] * network.TD)


t_price = 0
t_solve = 0
for iter in range(0, 30):

    t0 = time.time()
    #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively
      
    
        
    CG_status, obj, y_sol, cg_t_price, cg_t_solve = CG(rmp, network, y_ub)
    t_price += cg_t_price
    t_solve += cg_t_solve
    
    print(iter, obj, round(time.time() - t0, 2), round(t_price, 2), round(t_solve, 2))
    
    
    
            
        
    if lastObj > 0 and (obj - lastObj) / lastObj < 0.0001:
        break
        
    lastObj = obj
        
    OAcut = getOAcut(network)
    
    for a in network.links:
        rmp.add_constraint(rmp.mu[a] >= rmp.x[a]*OAcut['a'][a] + OAcut['b'][a])
    

tot_time = time.time() - tot_time

print("total time", round(tot_time, 2))
'''
print("\n")

for a in network.links:
    if a.y > 0:
        print("\n")
        for r in network.origins:
            if a in r.bush:
                print(a, r, bush_flows[(a,r)])
'''