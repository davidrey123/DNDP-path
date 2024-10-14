#---modules
import time
from src import Params
from src import Network
from docplex.mp.model import Model

rmp = Model()

c = [1, 2, 3, 4, 5]
g = [9, 8, 7]
d = [5, -1, 6, 3, 9]
e = [-1, 8, 3]

rmp.x = {i:rmp.continuous_var(lb=0,ub=1) for i in range(0, len(c))}
rmp.y = {j:rmp.continuous_var(lb=0,ub=1) for j in range(0, len(c))}


rmp.minimize(sum(rmp.x[i]*c[i] for i in range(0, len(c))) + sum(rmp.y[j]*g[j] for j in range(0, len(g))))


ub = 1e15
lb = -1e15

for iter in range(0, 10):
    rmp.solve(log_output=False)
    obj = rmp.objective_value
    
    lb = obj
    
    
    ll = Model()

    ll.x = {i:ll.continuous_var(lb=0,ub=1) for i in range(0, len(d))}
    ll.minimize(sum(ll.x[i]*d[i] for i in range(0, len(d))) + sum(rmp.y[j].solution_value*e[j] for j in range(0, len(e))))

    ll.solve(log_output=False)
    ll_obj = ll.objective_value
    
    obj_ll = sum(ll.x[i].solution_value * c[i] for i in range(0, len(c))) + sum(rmp.y[j].solution_value * e[j] for j in range(0, len(e)))
    ub = min(ub, obj_ll)
    
    # now take x, y point and make cut
    for i in range(0, len(d)):
        print("\t", i, ll.x[i].solution_value)
    rmp.add_constraint(sum(rmp.x[i] * d[i] for i in range(0, len(d))) + sum(rmp.y[j] * e[j] for j in range(0, len(e))) <= sum(ll.x[i].solution_value * d[i] for i in range(0, len(d))) + sum(rmp.y[j] * e[j] for j in range(0, len(e))))
    
    gap = ub-lb
    
    print(iter, lb, ub, gap, obj, ll_obj, obj_ll)
    
    if gap < 0.01:
        break