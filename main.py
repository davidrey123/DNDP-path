#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BC
from src import BPC

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)

run = 'BPC'

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
