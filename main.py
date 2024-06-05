# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:46:55 2024

@author: david
"""

#---modules
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'


network = Network.Network(net,ins,0.5,1e-0,1e-3)
print(net,ins)
bpc = BPC.BPC(network)
bpc.BB()
print()

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



'''
network = Network.Network(net,ins,0.5,1e-0,1e-3)
yy = {(7, 16): 0, (16, 7): 0, (19, 22): 1, (22, 19): 1, (11, 15): 1, (15, 11): 1, (9, 11): 0, (11, 9): 0, (13, 14): 0, (14, 13): 1}
y = {}
for a in network.links2:
    y[a] = yy[(a.start.id,a.end.id)]
print(y)
print(network.tapas('UE', y))
'''