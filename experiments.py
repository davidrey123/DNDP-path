#---modules
import time
from src import Params
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC_singleTree
from src import BPC_twoPhase
from src import BPC_nestedTree
from src import BPC_singleTree_link

t0_exp = time.time()

filename = 'results.txt'
f = open(filename, "w")

nets = ['SiouxFalls','BerlinMitteCenter','Anaheim','Barcelona']
#algs = ['BPC_singleTree_link','FS_NETS','Leblanc']
algs = ['Leblanc']

bprop = 0.25
inflate_trips = {'SiouxFalls':1,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}

headers = {'BPC_singleTree_link':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\',
           'BPC_singleTree':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\',
           'BPC_nestedTree':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\',
           'BPC_twoPhase':'& UB & Gap (\\%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE \\\\',
           'FS_NETS':'& UB & Gap (\\%) & Time (s) & MILP (s) & TAP (s) & SO & UE \\\\',
           'Leblanc':'& UB & Gap (\\%) & Time (s) & TAP (s) & SO & UE \\\\'}

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.2f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().CPLEX_threads))
f.write("Params: TAP_tol %.1f, SO_OA_cuts_tol %.1f, \n" % (Params.Params().min_gap,Params.Params().min_gap_SO_OA_cuts))
f.write("Budget/Total cost: %.2f\n" % (bprop))
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
            
        for nA2 in ['10','20']:
            
            for ID in range(1,2):
                
                ins = ins0+'_DNDP_'+nA2+'_'+str(ID)
                insshort = ins0+'\\_'+nA2+'\\_'+str(ID)
                print('running',ins)
             
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,1e-3,inflate_trips[net])
                
                if alg == 'BPC_singleTree_link':
                    bpc = BPC_singleTree_link.BPC_singleTree_link(network)
                    bpc.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    f.flush()
                    
                elif alg == 'BPC_singleTree':
                    bpc = BPC_singleTree.BPC_singleTree(network)
                    bpc.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    f.flush()                    
                    
                elif alg == 'BPC_twoPhase':
                    bpc = BPC_twoPhase.BPC_twoPhase(network)
                    bpc.BB()
                    f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    f.flush()
                
                elif alg == 'BPC_nestedTree':
                    bpc = BPC_nestedTree.BPC_nestedTree(network)
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
