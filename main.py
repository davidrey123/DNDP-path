#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC


net = 'SiouxFalls'
ins = 'SF_DNDP_20_1'

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

'''
#---start debug
y_fix = {(19, 22): 0, (22, 19): 1, (11, 15): 0, (15, 11): 1, (9, 11): 1, (11, 9): 1, (13, 14): 1, (14, 13): 0, (3, 11): 0, (11, 3): 1, (4, 10): 1, (10, 4): 0, (2, 13): 0, (13, 2): 0, (1, 18): 0, (18, 1): 1, (13, 18): 1, (18, 13): 1, (2, 12): 0, (12, 2): 0}
y = {}

for a in network.links2:
    y[a] = 0
    
    for yp in y_fix:
        if a.start.id==yp[0] and a.end.id==yp[1]:
            y[a] = y_fix[yp]
        
print(y)
tstt = network.tapas('UE',y)
print(tstt)
#---end debug
'''

run = 'BPC'
#run = 'Leblanc' 
#run = 'FS_NETS'
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
