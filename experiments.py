#---modules
import time
from src import Params
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC

t0_exp = time.time()

filename = 'results.txt'
f = open(filename, "w")

#nets = ['SiouxFalls','BerlinMitteCenter','Anaheim','Barcelona']
nets = ['SiouxFalls','EasternMassachusetts','BerlinMitteCenter']
algs = ['BPC','FS_NETS','Leblanc']

bprop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-1,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':1,'BerlinMitteCenter':2,'Anaheim':2,'Barcelona':2}

headers = {'BPC':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\',
           'FS_NETS':'& UB & Gap (\\%) & Time (s) & MILP (s) & TAP (s) & SO & UE \\\\',
           'Leblanc':'& UB & Gap (\\%) & Time (s) & TAP (s) & SO & UE \\\\'}

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.2f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().CPLEX_threads))
f.write("Params: TAP_tol %.1f, SO_OA_cuts_tol %.1f, OAcuts_tol %.1f, \n" % (Params.Params().min_gap,Params.Params().min_gap_SO_OA_cuts,Params.Params().OAcut_tol))
f.write("Budget/Total cost: %.2f\n" % (bprop))
f.write("Flow scaling: %s\n" % (scal_flow))
f.write("Trips inflation: %s\n" % (inflate_trips))
f.write('%s\n' % algs)

for alg in algs:
    print('-->',alg)
    t0_alg = time.time()
    
    header = 'Instance '+headers[alg]+'\n'
    f.write(header)
    f.flush()

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
            
            for ID in range(1,2):
                
                ins = ins0+'_DNDP_'+nA2+'_'+str(ID)
                insshort = ins0+'\\_'+nA2+'\\_'+str(ID)
                print('running',ins)
             
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,scal_flow[net],inflate_trips[net])
                
                if alg == 'BPC':
                    bpc = BPC.BPC(network)
                    bpc.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    f.flush()
                    
                elif alg == 'FS_NETS':
                    fs_nets = FS_NETS.FS_NETS(network)
                    fs_nets.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,fs_nets.UB,100*fs_nets.gap,fs_nets.rt,fs_nets.rt_MILP,fs_nets.rt_TAP,fs_nets.nSO,fs_nets.nUE))
                    f.flush()
                
                elif alg == 'Leblanc':
                    leblanc = Leblanc.Leblanc(network)
                    leblanc.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,leblanc.UB,100*leblanc.gap,leblanc.rt,leblanc.rt_TAP,leblanc.nSO,leblanc.nUE))
                    f.flush()
            
    total_alg_hours = (time.time() - t0_alg)/3600
    print('alg completed in %.2f hours' % total_alg_hours)
            
f.close()

total_hours = (time.time() - t0_exp)/3600
print('total runtime of %.2f hours' % total_hours)
