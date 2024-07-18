#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BC
from src import BPC_singleTree
from src import BPC_singleTree_UEcuts
from src import BPC_twoPhase
from src import BPC_nestedTree

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

net = 'Anaheim'
ins = 'A_DNDP_10_1'

net = 'Barcelona'
ins = 'B_DNDP_10_1'

net = 'BerlinMitteCenter'
ins = 'BMC_DNDP_10_1'

inflate_trips = 1
network = Network.Network(net,ins,0.5,1e-0,1e-3,inflate_trips)
print(net,ins)

'''
import time
t0 = time.time()
ytemp = {}
for a in network.links2:
    ytemp[a] = 1
tstt = network.tapas('UE', ytemp)
runtime = time.time() - t0
nOD = 0
for r in network.origins:
    for s in network.zones:
        if r.getDemand(s) > 0:
            nOD += 1
print(tstt,runtime)
print(len(network.nodes),len(network.links),nOD)

sotstt = network.tapas('SO', ytemp)
ueso = 100*(tstt - sotstt)/tstt
print(sotstt,ueso)

'''
run = 'BPC_singleTree'

if run == 'BPC_singleTree':
    bpc_singleTree = BPC_singleTree.BPC_singleTree(network)
    bpc_singleTree.BB()
    
elif run == 'BPC_singleTree_UEcuts':
    bpc_singleTree_UEcuts = BPC_singleTree_UEcuts.BPC_singleTree_UEcuts(network)
    bpc_singleTree_UEcuts.BB()    
    
elif run == 'BPC_twoPhase':
    bpc_twoPhase = BPC_twoPhase.BPC_twoPhase(network)
    bpc_twoPhase.BB()

elif run == 'BPC_nestedTree':
    bpc_nestedTree = BPC_nestedTree.BPC_nestedTree(network)
    bpc_nestedTree.BB() 

elif run == 'BC':
    print('note: need to check that no UE solution is skipped')    
    bc = BC.BC(network)
    bc.BB()
    
elif run == 'FS_NETS':
    fs_nets = FS_NETS.FS_NETS(network)
    fs_nets.BB()

elif run == 'Leblanc':
    leblanc = Leblanc.Leblanc(network)
    leblanc.BB()


print('\n-------------------------------------------')
print('Notes:')
print('need to rules to control OA cut generation')
print('if UB > can.LB > UB-SO-DNDP: what can we do?')
print('interdiction cuts: check usefulness in RMP')
print('-------------------------------------------')
