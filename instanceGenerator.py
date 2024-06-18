from src import TNTP
from src import Params
from numpy.random import choice
from numpy.random import rand

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
print('tstt',tstt)

TF = sum(a.x for a in tntp.links)

flows = []
for a in tntp.links:
    flows.append(a.x/TF)
    
newLinks = 20
nbInstances = 10

for n in range(nbInstances):
    
    #---reset link costs
    for a in tntp.links:
        a.cost = 0
       
    allConnected = 0
    nbNotEnough = 0
    while allConnected < newLinks:
    
        #---randomly choose new links 
        A2 = choice(tntp.links, newLinks, replace=False, p=flows)
        
        for a in A2:            
            
            tntp.resetY()
            a.y = 0
            
            #---check connectivity if only link a is closed
            connected = True
            for r in tntp.origins:
                tntp.dijkstras(r,'UE')
                for s in tntp.zones:
                    if r.getDemand(s) > 0:
                        if s.cost >= Params.INFTY - 1:
                            connected = False
                            
            if connected == True: 
                
                #---run UE-TAP to get TSTT and check loss wrt to base TSTT
                newtstt = tntp.tapas('UE',y)
                loss = 100*(newtstt - tstt)/tstt
                if loss <= 0.1:
                    print('not enough loss %.3f%%' % loss)
                    nbNotEnough += 1

                allConnected += 1
                #print('%d\t%d\t%d\t%.1f%%' % (a.id,a.start.id,a.end.id,loss))
            
            else:
                print('not connected')
                allConnected = 0
                
                #---set link flow to 0 to avoid this link being chosen again
                for b in tntp.links:
                    if a.id == b.id:
                        idx = tntp.links.index(b)                        
                        TF2 = TF - flows[idx]*TF
                        flows[idx] = 0
                        #print('flow set to 0 on link',b.id)
                        break
                
                #---update total flow to ensure that sum of flows = 1 and forms a probabilty distribution
                flows = [f*TF/TF2 for f in flows]
                TF = TF2
                break
            
        if allConnected == newLinks and nbNotEnough <= newLinks/3:
            
            #---test instance with all A2 links closed to check connectivity
            tntp.resetY()
            for a in A2:
                a.y = 0                
                
            connected = True
            for r in tntp.origins:
                tntp.dijkstras(r,'UE')
                for s in tntp.zones:
                    if r.getDemand(s) > 0:
                        if s.cost >= Params.INFTY - 1:
                            connected = False
            
            if connected == True:
            
                if net == 'Anaheim':
                    #---coefs for generating A2 link costs: idea is to generate costs correlated to fftt and C of the order of magniture 1e3
                    beta_t_ff = 1e3
                    beta_C = 1/1e1                    
                    
                    for a in A2:
                        scal = 0.8 + (rand(1)[0] * (1.2 - 0.8))
                        a.cost = round(scal*(beta_t_ff*a.t_ff + beta_C*a.C))
                        #print(a.id,a.cost)
                    
                    insName = 'A_DNDP_'+str(newLinks)+'_'+str(n+1)+'.txt'
                
                else:
                    print('unknown network')
                
                #---generate instance
                print(insName)
                tntp.writeInstance(insName,A2)
                
            else:
                print('not all connected')
                allConnected = 0
                
        else:
            #---either the network is not connected if all A2 links are closed 
            #   or not enough links generate significant TSTT loss when closed
            print('allConnected',allConnected)
            print('nbNotEnough',nbNotEnough)
            allConnected = 0
            nbNotEnough = 0