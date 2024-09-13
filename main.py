#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC
from src import BC

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

'''
net = 'BerlinMitteCenter'
ins = 'BMC_DNDP_10_1'

net = 'Anaheim'
ins = 'A_DNDP_10_1'

net = 'Barcelona'
ins = 'B_DNDP_20_2'

net = 'EasternMassachusetts'
ins = 'EM_DNDP_10_1'
'''

b_prop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}
network = Network.Network(net,ins,b_prop,1e-0,scal_flow[net],inflate_trips[net])
print(net,ins)

run = 'BPC'
run = 'BC'
#run = 'Leblanc' 
#run = 'FS_NETS'
print(run)

if run == 'BPC':
    bpc = BPC.BPC(network)
    bpc.BB()
    
elif run == 'BC':
    bc = BC.BC(network)
    bc.BB()    
    
elif run == 'FS_NETS':
    fs_nets = FS_NETS.FS_NETS(network)
    fs_nets.BB()
    
elif run == 'Leblanc':
    leblanc = Leblanc.Leblanc(network)
    leblanc.BB()
