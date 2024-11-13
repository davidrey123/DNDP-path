#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC
from src import BC
from src import OA_CNDP
from src import OA_CNDP_CG
#from src import HY_CNDP
from src import DuGP_CNDP
from src import CNDP_MILP
#import polytope as pc

#import numpy as np
from src import OA_CNDP_CS



'''
A = np.array([[1.0, 0.0],
              [0.0, 1.0],
              [-1.0, -0.0],
              [-0.0, -1.0]])

b = np.array([2.0, 1.0, 0.0, 0.0])

p = pc.Polytope(A, b)

print( [0.5, 0.5] in p)
print(pc.extreme(p))
'''





net = 'Braess'
ins = 'Braess_CNDP_1'

net = 'SiouxFalls'
ins = 'SF_CNDP_1'

net = 'EasternMassachusetts'
ins = 'EM_CNDP_30_1'

net = 'HarkerFriesz'
ins = 'HF_CNDP_1'

#net = 'BerlinMitteCenter'
#ins = 'BMC_CNDP_30_1'

#net = 'Anaheim'
#ins = 'A_CNDP_30_1'

b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':0.25}
print(net,ins)

inflate_cost = 1


print(inflate_trips[net], inflate_cost, scal_flow[net])

network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])

#network.tapas( "UE", {a:0 for a in network.links})
#test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True)
#test = HY_CNDP.HY_CNDP(network)
test = CNDP_MILP.CNDP_MILP(network, 5, 5, 20, inflate_cost)
#obj, tot_time, tap_time, iter, = test.solve()
test.solve()
#print(len(test.varlinks), " & ", network.TD, "&", test.getAvgLinkCost(), "& x &", round(obj, 1), "&", round(100*test.gap, 3), "\% &", round(test.tstt, 1), "&", round(tot_time, 2), "s &", round(tap_time 2), "s &", iter)



