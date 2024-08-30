#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import FS_NETS_piecewiselinear
from src import BC
from src import BPC_singleTree
from src import BPC_singleTree_link
from src import BPC_singleTree_UEcuts
from src import BPC_twoPhase
from src import BPC_nestedTree
from src import BPC_singleTree_piecewiselinear

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

'''net = 'Anaheim'
ins = 'A_DNDP_10_1'

net = 'Barcelona'
ins = 'B_DNDP_10_1'

net = 'BerlinMitteCenter'
ins = 'BMC_DNDP_10_1'
'''

b_prop = 0.25
inflate_trips = 1
network = Network.Network(net,ins,b_prop,1e-0,1e-3,inflate_trips)
print(net,ins)


#run = 'BPC_singleTree_link' 
run = 'FS_NETS' 
print(run)

if run == 'BPC_singleTree':
    bpc_singleTree = BPC_singleTree.BPC_singleTree(network)
    bpc_singleTree.BB()
    
elif run == 'BPC_singleTree_link':
    bpc_singleTree_link = BPC_singleTree_link.BPC_singleTree_link(network)
    bpc_singleTree_link.BB()
    
elif run == 'BPC_singleTree_piecewiselinear':
    
    bpc_singleTree_piecewiselinear = BPC_singleTree_piecewiselinear.BPC_singleTree_piecewiselinear(network)
    bpc_singleTree_piecewiselinear.BB()    
    
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
    
elif run == 'FS_NETS_piecewiselinear':
    fs_nets_piecewiselinear = FS_NETS_piecewiselinear.FS_NETS_piecewiselinear(network)
    fs_nets_piecewiselinear.BB()

elif run == 'Leblanc':
    leblanc = Leblanc.Leblanc(network)
    leblanc.BB()


print('\n-------------------------------------------')
print('Notes:')
print('need to rules to control OA cut generation')
print('if UB > can.LB > UB-SO-DNDP: what can we do?')
print('interdiction cuts: check usefulness in RMP')
print('-------------------------------------------')
