#---modules
import time
from src import Params
from src import Network
from src import BPC

t0_exp = time.time()

filename = 'sensitivity.txt'
f = open(filename, "w")

#nets = ['SiouxFalls','BerlinMitteCenter','Anaheim','Barcelona']
nets = ['SiouxFalls','EasternMassachusetts','BerlinMitteCenter']

configs = [{'min_gap_SO_OA_cuts':1e-1,'OAcut_tol':0.01, 'nInitKNP':1},
           {'min_gap_SO_OA_cuts':1e-2,'OAcut_tol':0.01, 'nInitKNP':1},
           {'min_gap_SO_OA_cuts':1e-3,'OAcut_tol':0.01, 'nInitKNP':1},
           {'min_gap_SO_OA_cuts':1e-1,'OAcut_tol':0.005, 'nInitKNP':1},
           {'min_gap_SO_OA_cuts':1e-2,'OAcut_tol':0.005, 'nInitKNP':1},
           {'min_gap_SO_OA_cuts':1e-3,'OAcut_tol':0.005, 'nInitKNP':1}]

bprop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-1,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':1,'BerlinMitteCenter':2,'Anaheim':2,'Barcelona':2}

headers = {'BPC':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\'}

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.2f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().CPLEX_threads))
f.write("Params: TAP_tol %.1f, SO_OA_cuts_tol %.1f, OAcuts_tol %.1f, \n" % (Params.Params().min_gap,Params.Params().min_gap_SO_OA_cuts,Params.Params().OAcut_tol))
f.write("Budget/Total cost: %.2f\n" % (bprop))
f.write("Flow scaling: %s\n" % (scal_flow))
f.write("Trips inflation: %s\n" % (inflate_trips))

header = 'Instance '+headers['BPC']+'\n'
f.write(header)
f.flush()

for config in configs:
    print('-->',config)
    t0_config = time.time()
    
    f.write("config: %s\n" % config)

    for net in nets:
        
        if net == 'SiouxFalls':
            ins0 = 'SF'
            
        elif net == 'Anaheim':
            ins0 = 'A'
            
        elif net == 'Barcelona':
            ins0 = 'B'
            
        elif net == 'BerlinMitteCenter':
            ins0 = 'BMC'  

        elif net == 'EasternMassachusetts':
            ins0 = 'EM'              
            
        for nA2 in ['10','20']:
            
            for ID in range(1,4):
                
                ins = ins0+'_DNDP_'+nA2+'_'+str(ID)
                insshort = ins0+'\\_'+nA2+'\\_'+str(ID)
                print('running',ins)
             
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,scal_flow[net],inflate_trips[net])

                bpc = BPC.BPC(network)
                
                bpc.params.min_gap_SO_OA_cuts = config['min_gap_SO_OA_cuts']
                bpc.params.OAcut_tol = config['OAcut_tol']
                bpc.nInitKNP = round(config['nInitKNP']*len(network.links2))
                
                bpc.BB()
                f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                f.flush()
            
    total_config_hours = (time.time() - t0_config)/3600
    print('config completed in %.2f hours' % total_config_hours)
            
f.close()

total_hours = (time.time() - t0_exp)/3600
print('total runtime of %.2f hours' % total_hours)
