#---modules
import time
from src import Params
from src import Network
from docplex.mp.model import Model

def getOAcut(network):
    OAcut = {'a':{}, 'b':{}}

    for a in network.links:
               
        OAcut['a'][a] = a.getTravelTime(a.x,'SO')
        OAcut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
        

    return OAcut
    
def pricing(duals, network, paths):
        
    t0_pricing = time.time()

    new = 0
    minrc = 1e15
    
    for a in network.links:
        a.dual = duals['link'][a]

    for r in network.origins:
        network.dijkstras(r,'RC')

        for s in network.zones:

            if r.getDemand(s) > 0:

                rc = - duals['dem'][(r,s)] + s.cost                    

                if rc <= - 0.0001:
                    p = network.trace(r,s)
                    paths[r][s].append(p)
                    new += 1

                if rc < minrc:
                    minrc = rc

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

def rmp_path(network, paths, OAcuts, y_ub):
        
    #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively

    t0_RMP = time.time()
    rmp = Model()

    rmp.x = {a:rmp.continuous_var(lb=0,ub=network.TD) for a in network.links}
    rmp.h = {p:rmp.continuous_var(lb=0) for p in getPaths(network, paths)}
    rmp.mu = {a:rmp.continuous_var(lb=0) for a in network.links}

    rmp.y = {a:rmp.continuous_var(lb=0,ub=y_ub) for a in network.links2}

    rmp.add_constraint(sum(rmp.y[a] * a.cost for a in network.links2) <= network.B)

    for a in network.links2:
        rmp.add_constraint(rmp.x[a] <= rmp.y[a] * network.TD)

    for r in network.origins:
        for s in network.zones:
            if r.getDemand(s) > 0:
                rmp.add_constraint(sum(rmp.h[p] for p in paths[r][s]) >= r.getDemand(s), 'dem_%d_%d' % (r.id,s.id))                    

    for a in network.links:
        rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in getPaths(network, paths) if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))

    for OAcut in OAcuts:
        #---OA cuts
        for a in network.links:
            rmp.add_constraint(rmp.mu[a] >= rmp.x[a]*OAcut['a'][a] + OAcut['b'][a])

    
    rmp.minimize(sum(rmp.mu[a] for a in network.links))


    rmp.solve(log_output=False)

    #if self.params.PRINT_BB_INFO:
    #    print('nb of paths: %d, cplex time: %.1f, rmp time: %.1f' % (len(can.getPaths()),rmp.solve_details.time,(time.time() - t0_RMP)))


    if rmp.solve_details.status == 'infeasible' or rmp.solve_details.status == 'integer infeasible':
        return 'infeasible',1e15,{},None

    else:
        OFV = rmp.objective_value
        RMP_status = rmp.solve_details.status

        yopt = {}
        for a in network.links2:
            yopt[a] = rmp.y[a].solution_value
            
            
        for a in network.links:
            a.x = rmp.x[a].solution_value

        dual_link = {} 
        for a in network.links:                    
            dual_link[a] = max(rmp.get_constraint_by_name('link_%d_%d' % (a.start.id,a.end.id)).dual_value,0)                    

        dual_dem = {}
        for r in network.origins:
            for s in network.zones:
                if r.getDemand(s) > 0:
                    dual_dem[(r,s)] = max(rmp.get_constraint_by_name('dem_%d_%d' % (r.id,s.id)).dual_value,0)

        duals = {'link':dual_link,'dem':dual_dem}

        return RMP_status,OFV,yopt,duals

def CG(network, paths, OAcuts, y_ub):

    nCG = 0
    conv = False
    
    

    while conv == False:        

        RMP_status,OFV,yRMP,duals = rmp_path(network, paths, OAcuts, y_ub)

        if RMP_status == 'infeasible':
            CG_status = 'infeasible'
            break

        minrc = pricing(duals, network, paths)
        
        print("\t", OFV, minrc)

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

        return CG_status,OFV,yRMP

    else:
        return CG_status,1e15,yRMP


  

    
    
net = 'SiouxFalls'
ins = 'SF_DNDP_20_1'

net = 'EasternMassachusetts'
ins = 'EM_DNDP_10_1'

tot_time = time.time()



b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}
network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])




OAcuts = []

bush_flows = {}

lastObj = 0

y_ub = 0

for a in network.links2:
    a.y = 0
    
    
    
paths = {r:{s:[] for s in network.zones} for r in network.origins}
    
for r in network.origins:        
    network.dijkstras(r,'UE')

    for s in network.zones:

        if r.getDemand(s) > 0:

            p = network.trace(r,s)
            paths[r][s].append(p)

for iter in range(0, 30):

    t0 = time.time()
    #---to do: recode such that lp is setup once and only new paths and cuts are added iteratively
      
    
        
    CG_status, obj, y_sol = CG(network, paths, OAcuts, y_ub)
    
    print(iter, obj, round(time.time() - t0, 2))
    
    
    
            
        
    if lastObj > 0 and (obj - lastObj) / lastObj < 0.0001:
        break
        
    lastObj = obj
        
    OAcuts.append(getOAcut(network))
    

tot_time = time.time() - tot_time

print("total time", round(tot_time, 2))
'''
print("\n")

for a in network.links:
    if a.y > 0:
        print("\n")
        for r in network.origins:
            if a in bushes[r]:
                print(a, r, bush_flows[(a,r)])
'''