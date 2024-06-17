from src import TNTP
from numpy.random import choice

net = 'SiouxFalls'
net = 'Anaheim'

tntp = TNTP.TNTP(net)
print(net)

print(len(tntp.zones))
print(len(tntp.nodes))
print(len(tntp.links))
print(len(tntp.origins))

y = {}
tstt = tntp.tapas('UE',y)
print(tstt)

TF = sum(a.x for a in tntp.links)

flows = []
for a in tntp.links:
    flows.append(a.x/TF)
    
newLinks = 10

A2 = choice(tntp.links, newLinks, replace=False, p=flows)


for a in A2:
    print(a.id,a.start.id,a.end.id)
    tntp.resetY()
    a.y = 0
    newtstt = tntp.tapas('UE',y)
    loss = 100*(newtstt - tstt)/tstt
    print(newtstt,loss)

