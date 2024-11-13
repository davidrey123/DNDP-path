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
import math
#import polytope as pc

#import numpy as np
from src import OA_CNDP_CS


CG = True

net = 'SiouxFalls'
ins = 'SF_CNDP_'


#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_'


b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':1}
print(net,ins)

#network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
#test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True)


inflate_cost = 1
scale_dem = 1


f = open("experiments"+net+".txt", "w")


for j in range (0, 2):
    yvars = 10+j*20
    
    
    for k in range (1, 3):
        
        if k == 1:
            inflate_cost = 1
            scale_dem = 1
        elif k == 2:
            inflate_cost = 1
            scale_dem = 2
        elif k == 3:
            inflate_cost = 5
            scale_dem = 1
            
        for i in range (1, 3):
            network = Network.Network(net,ins+str(j)+"_"+str(i),b_prop,1e-0,scal_flow[net],scale_dem * inflate_trips[net])
            test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True, useCG=CG)
            obj, tot_time, tap_time, iter, = test.solve()

            scientific_format = "{:.2e}".format(test.getAvgLinkCost())
            f.write(str(len(test.varlinks)) + " & "+ str(round(network.TD,1))+ " & " +  str(scientific_format) +  " & "+str(i)+" & " + str(round(obj, 1))+ " & "+ str(round(100*test.gap, 3))+ "\% & " + str(round(test.tstt, 1)) + " & " + str(round(tot_time, 2)) + "s & "+ str(round(tap_time, 2))+ "s &" + str(iter))
            f.write("\n")

            f.flush()
            
f.close()