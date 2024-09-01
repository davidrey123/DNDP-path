#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC
from src import BPC_singleTree_link


net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

'''
net = 'Anaheim'
ins = 'A_DNDP_10_1'

net = 'Barcelona'
ins = 'B_DNDP_10_1'

net = 'BerlinMitteCenter'
ins = 'BMC_DNDP_10_1'
'''

b_prop = 0.5
inflate_trips = 1
network = Network.Network(net,ins,b_prop,1e-0,1e-3,inflate_trips)
print(net,ins)


run = 'BPC_singleTree_link' 
run = 'BPC'
#run = 'Leblanc' 
print(run)

if run == 'BPC':
    bpc = BPC.BPC(network)
    bpc.BB()
    
elif run == 'BPC_singleTree_link':
    bpc_singleTree_link = BPC_singleTree_link.BPC_singleTree_link(network)
    bpc_singleTree_link.BB()
    
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
