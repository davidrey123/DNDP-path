#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BC
from src import BPC

net = 'SiouxFalls'
ins = 'SF_DNDP_10_4'


network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)
bc = BC.BC(network)
bc.BB()
print()

'''
network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)
bpc = BPC.BPC(network)
bpc.BB()
print()
'''
'''
network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)
leblanc = Leblanc.Leblanc(network)
leblanc.BB()
print()
'''

network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)
fs_nets = FS_NETS.FS_NETS(network)
fs_nets.BB()
print()
