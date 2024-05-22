# Created on : May 22, 2024, 11:13:03 AM
# Author     : michaellevin

from src import Network
from src import YDict

def test():
    net = 'SiouxFalls'
    ins = 'SF_DNDP_10_6'
    network = Network.Network(net,ins,0.5,1e-0,1e-3)
    
    y1 = {}

    for a in network.links2:
        y1[a] = 1

    y2 = {}
    count = 0
    for a in network.links2:

        if (count % 2 == 0):
            y2[a] = 1
        else:
            y2[a] = 0
        count += 1

    y3 = {}
    count = 0
    for a in network.links2:

        if (count % 3 == 0):
            y3[a] = 1
        else:
            y3[a] = 0
        count += 1

    print(y1)

    test = YDict.YDict()

    test.insertSO(y1, 100)
    test.insertSO(y2, 200)

    print("y1 is", test.getSO(y1))
    print("y2 is", test.getSO(y2))
    print("y3 is", test.getSO(y3))
    test.insertSO(y3, 300)

    print("y3 is", test.getSO(y3))
