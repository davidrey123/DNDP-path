#---modules
import time
from src import Params
from src import Network
import CGbush
from docplex.mp.model import Model

def getOAcut(network):
    OAcut = {'a':{}, 'b':{}}

    for a in network.links:
               
        OAcut['a'][a] = a.getTravelTime(a.x,'SO')
        OAcut['b'][a] = - pow(a.x, 2) * a.getDerivativeTravelTime(a.x)
  
        #print(a, a.x, OAcut['a'][a], OAcut['b'][a])

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

lastObj = 0
obj = 0
lastCGobj = 0

bushes = {}


for a in network.links2:
    a.y = 0

for r in network.origins:
    bushes[r] = CGbush.CGbush(r, network)
    
    
            

for iter in range(0, 30): # OA loop

    t0 = time.time()
    
    for cg_iter in range(0, 30): # CG loop
    
        t1 = time.time()
        cp = Model()
        cp.x = {a:cp.continuous_var(lb=0,ub=network.TD) for a in network.links}
        cp.xc = {}

        cp.xc = {(a,r):cp.continuous_var() for r in network.origins for a in bushes[r].linkflows.keys()}

        cp.y = {a:cp.continuous_var(lb=0, ub=0) for a in network.links2}
        cp.mu = {a:cp.continuous_var(lb=0,ub=1E10) for a in network.links}

        mu_cons = cp.add_constraint(-sum(cp.y[a] * a.cost for a in network.links2) >= -network.B)

        phi_cons = {}

        for a in network.links2:
            phi_cons[a] = cp.add_constraint(-cp.x[a] + cp.y[a] * network.TD >= 0)

        gamma_plus_cons = {}
        gamma_minus_cons = {}

        for a in network.links:
            gamma_plus_cons[a] = cp.add_constraint(cp.x[a] - sum(cp.xc[(a,r)] for r in network.origins if bushes[r].hasLink(a)) >= 0)
            gamma_minus_cons[a] = cp.add_constraint(-cp.x[a] + sum(cp.xc[(a,r)] for r in network.origins if bushes[r].hasLink(a)) >= 0)


        eta_plus_cons = {}
        eta_minus_cons = {}

        for r in network.origins:  
            for i in network.nodes:                    


                if i.id == r.id:
                    dem = -sum(r.getDemand(s) for s in network.zones)                
                else:
                    dem = r.getDemand(i)


                eta_plus_cons[(i, r)] = cp.add_constraint(sum(cp.xc[(a,r)] for a in i.incoming if bushes[r].hasLink(a)) - sum(cp.xc[(a,r)] for a in i.outgoing if bushes[r].hasLink(a)) >= dem)
                eta_minus_cons[(i, r)] = cp.add_constraint(-sum(cp.xc[(a,r)] for a in i.incoming if bushes[r].hasLink(a)) + sum(cp.xc[(a,r)] for a in i.outgoing if bushes[r].hasLink(a)) >= -dem)

        lambda_cons = {}

        for OAcut_idx in range(0, len(OAcuts)):
            OAcut = OAcuts[OAcut_idx]
            for a in network.links:
                lambda_cons[(OAcut_idx, a)] = cp.add_constraint(cp.mu[a] - cp.x[a]*OAcut['a'][a] >=  OAcut['b'][a])


        cp.minimize(sum(cp.mu[a] for a in network.links))

        cp.solve(log_output=False)
        
        new_col = False

        if cp.solve_details.status == 'infeasible':
            print(cg_iter, "INFEASIBLE")
            break

        obj = cp.objective_value

        print("\t", cg_iter, obj, round(time.time() - t1, 2))
        
        
        for a in network.links:

            for r in network.origins:
                if bushes[r].hasLink(a):
                    bushes[r].linkflows[a] = cp.xc[(a,r)].solution_value

                    #if a.start.id==10 and a.end.id==16:
                    #    print(a, r, a.x, bush_flows[(a,r)])

        
        for r in network.zones:
            for a in network.links:
            
                
                price = -(-gamma_plus_cons[a].dual_value + gamma_minus_cons[a].dual_value)
                price += -(eta_plus_cons[(a.end, r)].dual_value - eta_plus_cons[(a.start, r)].dual_value)
                price += -(-eta_minus_cons[(a.end, r)].dual_value + eta_minus_cons[(a.start, r)].dual_value)

                if bushes[r].hasLink(a):
                    bushes[r].link_RC[a] = price
                elif price < -0.0001:
                    new_col = True
                    bushes[r].addLink(a)
                    
                    
            bushes[r].processNewLinks()
        
        if not new_col:
            break

        lastCGobj = obj

    print(iter, obj, round(time.time() - t0, 2))

    for a in network.links:
        a.x = cp.x[a].solution_value

        
    for a in network.links2:
        a.y = cp.y[a].solution_value

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