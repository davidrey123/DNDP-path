#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BC
from src import BPC_singleTree
from src import BPC_twoPhase
from src import BPC_nestedTree

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

#net = 'Anaheim'
#ins = 'A_DNDP_10_1'

network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)

run = 'BPC_nestedTree'

if run == 'BPC_singleTree':
    bpc_singleTree = BPC_singleTree.BPC_singleTree(network)
    bpc_singleTree.BB()
    
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
print('Try Beckmann OA cuts')
print('-------------------------------------------')