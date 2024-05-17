# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 18:46:55 2024

@author: david
"""

#---modules
from src import Network
from src import Leblanc

net = 'SiouxFalls'
ins = 'SF_DNDP_10_1'

network = Network.Network(net,ins,0.5,1e-0,1e-3)
#network = Network.Network("grid3",1,1,1,1,1)
print(net,ins)
print()

'''
y1 = {}

for ij in network.links2:
    y1[ij] = 1


y2 = dict(y1)
y3 = dict(y2)
y4 = dict(y3)

id = 0

for ij in y2:
    id = id + 1
    if id % 2 == 0:
        y2[ij] = 0
        
id = 0

for ij in y3:
    id = id + 1
    if id % 2 == 0:
        y3[ij] = 0
        
id = 0

for ij in y4:
    id = id + 1
    if id % 4 == 0:
        y4[ij] = 0

network.tapas('SO', y1)

network.tapas('SO', y2)

network.tapas('SO', y3)

network.tapas('SO', y4)


'''

timelimit = 600

leblanc = Leblanc.Leblanc(network, timelimit)

leblanc.BB()


'''
yy = {(11, 15): 1, (15, 11): 1, (12, 14): 1, (14, 12): 1, (10, 19): 1, (19, 10): 1, (2, 13): 0, (13, 2): 0, (2, 12): 0, (12, 2): 0}
y = {}
for a in network.links2:
    y[a] = yy[(a.start.id,a.end.id)]
print(y)
print(network.tapas('UE', y))
'''
