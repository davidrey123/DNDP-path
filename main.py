#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC

'''
net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

net = 'BerlinMitteCenter'
ins = 'BMC_DNDP_10_1'

net = 'Anaheim'
ins = 'A_DNDP_10_1'

net = 'Barcelona'
ins = 'B_DNDP_20_2'
'''

net = 'EasternMassachusetts'
ins = 'EM_DNDP_10_1'

b_prop = 0.5
scal_flow = 1e-1 #---default 1e-3
inflate_trips = 1
network = Network.Network(net,ins,b_prop,1e-0,scal_flow,inflate_trips)
print(net,ins)

run = 'BPC'
run = 'Leblanc' 
run = 'FS_NETS'
print(run)

if run == 'BPC':
    bpc = BPC.BPC(network)
    bpc.BB()
    
elif run == 'FS_NETS':
    fs_nets = FS_NETS.FS_NETS(network)
    fs_nets.BB()
    
elif run == 'Leblanc':
    leblanc = Leblanc.Leblanc(network)
    leblanc.BB()


print('\n-------------------------------------------')
print('Notes:')
print('need rules to control OA cut generation: stop once enough?')
print('if UB > can.LB > UB-SO-DNDP: what can we do?')
print('-------------------------------------------')
