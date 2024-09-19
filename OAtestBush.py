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

net = 'EasternMassachusetts'
ins = 'EM_DNDP_10_1'



# test on BMC_DNDP_10_1 it has a difference between FS_NETS relaxed MILP and path-based CG


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
    
y_ub = 0

if y_ub == 1:
    for a in network.links2:
        a.y = 1
            
cp = Model()
cp.x = {a:cp.continuous_var(lb=0,ub=network.TD) for a in network.links}
cp.xc = {(a,r):cp.continuous_var(lb=0,ub=network.TD) for r in network.origins for a in bushes[r].linkflows.keys()}
cp.y = {a:cp.continuous_var(lb=0, ub=y_ub) for a in network.links2}
cp.mu = {a:cp.continuous_var(lb=0,ub=1E10) for a in network.links}

mu_cons = cp.add_constraint(-sum(cp.y[a] * a.cost for a in network.links2) >= -network.B)

phi_cons = {}

#for a in network.links:
#    cp.add_constraint(cp.mu[a] >= cp.x[a]*a.t_ff)
    
cp.minimize(sum(cp.mu[a] for a in network.links))

for a in network.links2:
    cp.add_constraint(-cp.x[a] + cp.y[a] * network.TD >= 0)
    
    



gamma_plus_cons = {}
gamma_minus_cons = {}

for a in network.links:
    gamma_plus_cons[a] = cp.add_constraint(cp.x[a] - sum(cp.xc[(a,r)] for r in network.origins if bushes[r].contains(a)) >= 0)
    gamma_minus_cons[a] = cp.add_constraint(-cp.x[a] + sum(cp.xc[(a,r)] for r in network.origins if bushes[r].contains(a)) >= 0)
    
    
eta_plus_cons = {}
eta_minus_cons = {}

OAcut_idx = 0

for r in network.origins:  
    for i in network.nodes:                    

        if i.id == r.id:
            dem = -sum(r.getDemand(s) for s in network.zones)                
        else:
            dem = r.getDemand(i)

        bush_inc = i.getBushIncoming(bushes[r])
        bush_out = i.getBushOutgoing(bushes[r])

        if len(bush_inc) > 0 or len(bush_out) > 0:
            eta_plus_cons[(i, r)] = cp.add_constraint(sum(cp.xc[(a,r)] for a in bush_inc) - sum(cp.xc[(a,r)] for a in bush_out) >= dem)
            eta_minus_cons[(i, r)] = cp.add_constraint(-sum(cp.xc[(a,r)] for a in bush_inc) + sum(cp.xc[(a,r)] for a in bush_out) >= -dem)


        
for iter in range(0, 30): # OA loop

    t0 = time.time()
    
    for cg_iter in range(0, 30): # CG loop
    
        t1 = time.time()

        cp.solve(log_output=False)
        
        new_col = False

        if cp.solve_details.status == 'infeasible':
            print(cg_iter, "INFEASIBLE")
            break

        obj = cp.objective_value

        
        
        
        for a in network.links:

            for r in network.origins:
                if bushes[r].contains(a):
                    bushes[r].linkflows[a] = cp.xc[(a,r)].solution_value

                    #if a.start.id==10 and a.end.id==16:
                    #    print(a, r, a.x, bush_flows[(a,r)])

        best_price = 0
        
        for r in network.origins:
            best_r_price = 0
            best_link = None
            
            added_vars = []
            
            for a in network.links:
            
                price = -(-gamma_plus_cons[a].dual_value + gamma_minus_cons[a].dual_value)

                if (a.end, r) in eta_plus_cons:
                    price += -(eta_plus_cons[(a.end, r)].dual_value) 
                    price += -(-eta_minus_cons[(a.end, r)].dual_value)
                #else:
                #    print("\t\tmissing", r, a.end)

                if (a.start, r) in eta_plus_cons:
                    price += -(- eta_plus_cons[(a.start, r)].dual_value)
                    price += -(eta_minus_cons[(a.start, r)].dual_value)
                #else:
                #    print("\t\tmissing", r, a.start)

                if bushes[r].contains(a):
                    bushes[r].link_RC[a] = price
                elif price < best_r_price - 0.0001:
                    new_col = True
                    '''
                    best_link = a
                    best_r_price = price
                    '''
                    bushes[r].addLink(a)
                    best_price = min(best_price, price)  
                    added_vars.append(a)
                    
            '''
            if best_link is not None:
                bushes[r].addLink(best_link)
                best_price = min(best_price, best_r_price)        
                bushes[r].processNewLinks()
            '''
                
            removed_vars = bushes[r].processNewLinks()      
            
            
            for a in added_vars:
                new_xc = None
                if (a,r) in cp.xc:
                    new_xc = cp.xc[(a,r)]
                else:
                    new_xc = cp.continuous_var(lb=0,ub=network.TD)
                    cp.xc[(a,r)] = new_xc
                
                #gamma_plus_cons[a] = cp.add_constraint(cp.x[a] - sum(cp.xc[(a,r)] for r in network.origins if bushes[r].contains(a)) >= 0)
                #gamma_minus_cons[a] = cp.add_constraint(-cp.x[a] + sum(cp.xc[(a,r)] for r in network.origins if bushes[r].contains(a)) >= 0)
                
                gamma_plus_cons[a].lhs.add_term(new_xc, -1)
                gamma_minus_cons[a].lhs.add_term(new_xc, 1)
            
                #eta_plus_cons[(i, r)] = cp.add_constraint(sum(cp.xc[(a,r)] for a in bush_inc) - sum(cp.xc[(a,r)] for a in bush_out) >= dem)
                #eta_minus_cons[(i, r)] = cp.add_constraint(-sum(cp.xc[(a,r)] for a in bush_inc) + sum(cp.xc[(a,r)] for a in bush_out) >= -dem)
                
                eta_plus_cons[(a.end, r)].lhs.add_term(new_xc, 1)
                eta_plus_cons[(a.start, r)].lhs.add_term(new_xc, -1)
                eta_minus_cons[(a.end, r)].lhs.add_term(new_xc, -1)
                eta_minus_cons[(a.start, r)].lhs.add_term(new_xc, 1)
            
            
            for a in removed_vars:
            
                old_xc = cp.xc[(a,r)]
                
                # remove from constraints
                gamma_plus_cons[a].lhs.remove_term(old_xc)
                gamma_minus_cons[a].lhs.remove_term(old_xc)
                
                eta_plus_cons[(a.end, r)].lhs.remove_term(old_xc)
                eta_plus_cons[(a.start, r)].lhs.remove_term(old_xc)
                eta_minus_cons[(a.end, r)].lhs.remove_term(old_xc)
                eta_minus_cons[(a.end, r)].lhs.remove_term(old_xc)
            
                # actually delete variable
                #cp.get_cplex().variables.delete(old_xc._index)
                cp.add_constraint(old_xc == 0)
                del cp.xc[(a,r)]
            
        
        print("\t", cg_iter, obj, best_price, round(time.time() - t1, 2))
        
        if not new_col or abs(lastCGobj - obj) < 0.0001:
            break

        lastCGobj = obj

    print(iter, obj, round(time.time() - t0, 2))

    for a in network.links:
        a.x = cp.x[a].solution_value

    if lastObj > 0 and (obj - lastObj) / lastObj < 0.0001:
        for a in network.links2:
            a.y = cp.y[a].solution_value
        break

    lastObj = obj

    OAcut = getOAcut(network)
    
    for a in network.links:
        cp.add_constraint(cp.mu[a] - cp.x[a]*OAcut['a'][a] >=  OAcut['b'][a])
        OAcut_idx += 1
    
    
tot_time = time.time() - tot_time

print("total time", round(tot_time, 2))

'''
print("\n")

for a in network.links:
    print("\n")

    print(a, r, a.y, a.x)
'''