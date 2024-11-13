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

nets = ['SiouxFalls', 'EasternMassachusetts', 'Anaheim', 'BerlinMitteCenter']
inss = ['SF_CNDP_','EM_CNDP_','A_CNDP_','BMC_CNDP_']

scale_dems = [
[1, 0.5],
[1, 2],
[1, 2],
[1, 2]
]

inflate_costs = [
[1,5],
[1,5],
[1,5],
[1,5],
]

yvarss = [
[10,30],
[10,30],
[30,60],
[30,60]
]

#net = 'EasternMassachusetts'
#ins = 'EM_CNDP_'


b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3, 'Braess':1, 'HarkerFriesz':1}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2, 'Braess':1, 'HarkerFriesz':1}

#network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
#test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True)


inflate_cost = 1
scale_dem = 1




for n in range (0, 1):
    net = nets[n]
    f = open("experiments_"+net+".txt", "w")
    
    for j in range (0, 2):
        yvars = yvarss[n][j]


        for k in range (1, 4):
            if k == 1:
                inflate_cost = inflate_costs[n][0]
                scale_dem = scale_dems[n][0]
            elif k == 2:
                inflate_cost = inflate_costs[n][0]
                scale_dem = scale_dems[n][1]
            elif k == 3:
                inflate_cost = inflate_costs[n][1]
                scale_dem = scale_dems[n][0]

            for i in range (1, 2):
                actual_ins = inss[n]+str(yvars)+"_"+str(i)
                print("\n\n\n", actual_ins, yvars, scale_dem, inflate_cost, CG, "\n\n\n")

                network = Network.Network(net,actual_ins,b_prop,1e-0,scal_flow[net],scale_dem * inflate_trips[net])
                test = OA_CNDP_CG.OA_CNDP_CG(network, inflate_cost, useLinkVF=True, useCG=CG)
                
                if net == "SiouxFalls":
                    test.params.min_gap = 1E-3
                elif net == "EasternMassachusetts":
                    test.params.min_gap = 4E-3
                    
                obj, tot_time, tap_time, iter, = test.solve()

                scientific_format = "{:.2e}".format(test.getAvgLinkCost())
                f.write(str(len(test.varlinks)) + " & "+ str(round(network.TD,1))+ " & " +  str(scientific_format) +  " & "+str(i)+" & " + str(round(obj, 1))+ " & "+ str(round(100*test.gap, 3))+ "\% & " + str(round(test.tstt, 1)) + " & " + str(round(tot_time, 2)) + "s & "+ str(round(tap_time, 2))+ "s &" + str(iter))
                f.write("\n")

                f.flush()
            
f.close()