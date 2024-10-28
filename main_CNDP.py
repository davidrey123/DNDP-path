#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC
from src import BC
from src import OA_CNDP
from src import OA_CNDP_CG
from src import HY_CNDP
import polytope as pc

import numpy as np



A = np.array([[1.0, 0.0],
              [0.0, 1.0],
              [-1.0, -0.0],
              [-0.0, -1.0]])

b = np.array([2.0, 1.0, 0.0, 0.0])

p = pc.Polytope(A, b)

print( [0.5, 0.5] in p)
print(pc.extreme(p))




net = 'SiouxFalls'
ins = 'SF_CNDP_1'

#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_10_1'

b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}
network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
print(net,ins)

#test = OA_CNDP.OA_CNDP(network)
#test = OA_CNDP_CG.OA_CNDP_CG(network)
test = HY_CNDP.HY_CNDP(network)
test.solve()


