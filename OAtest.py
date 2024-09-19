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

net = 'SiouxFalls'
ins = 'SF_DNDP_20_1'

#net = 'EasternMassachusetts'
#ins = 'EM_DNDP_10_1'

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
    a.y = 1


cp = Model()
cp.x = {a:cp.continuous_var(lb=0,ub=network.TD) for a in network.links}
cp.xc = {(a,r):cp.continuous_var() for a in network.links for r in network.origins}

cp.y = {a:cp.continuous_var(lb=0, ub=y_ub) for a in network.links2}
cp.mu = {a:cp.continuous_var(lb=0,ub=1e10) for a in network.links}

cp.add_constraint(sum(cp.y[a] * a.cost for a in network.links2) <= network.B)

for a in network.links2:
    cp.add_constraint(cp.x[a] <= cp.y[a] * network.TD)
    
    
for a in network.links:
    cp.add_constraint(sum(cp.xc[(a,r)] for r in network.origins) == cp.x[a])
    
    
for i in network.nodes:                    
    for r in network.origins:            

        if i.id == r.id:
            dem = - sum(r.getDemand(s) for s in network.zones)                
        elif isinstance(i, type(r)) == True:
            dem = r.getDemand(i)
        else:
            dem = 0

        cp.add_constraint(sum(cp.xc[(a,r)] for a in i.incoming) - sum(cp.xc[(a,r)] for a in i.outgoing) == dem)


for a in network.links:
    cp.add_constraint(cp.mu[a] >= cp.x[a]*a.t_ff)
        
cp.minimize(sum(cp.mu[a] for a in network.links))


for iter in range(0, 30):

    t0 = time.time()

    cp.solve(log_output=False)
    
    if cp.solve_details.status == 'infeasible':
        print(iter, "INFEASIBLE")
        break
        
    obj = cp.objective_value
    
    print(iter, obj, round(time.time() - t0, 2))
    
    
    
    for a in network.links:
        a.x = cp.x[a].solution_value
        
        for r in network.origins:
            bush_flows[(a,r)] = cp.xc[(a,r)].solution_value
            
            
        
    if lastObj > 0 and (obj - lastObj) / lastObj < 0.0001:
        for a in network.links2:
            a.y = cp.y[a].solution_value
        break
        
    lastObj = obj
        
    OAcut = getOAcut(network)
    
    for a in network.links:
        cp.add_constraint(cp.mu[a] >= cp.x[a]*OAcut['a'][a] + OAcut['b'][a])
    

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