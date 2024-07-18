#---modules
import time
from src import Params
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BPC_singleTree
from src import BPC_twoPhase
from src import BPC_nestedTree

t0_exp = time.time()

filename = 'results.txt'
f = open(filename, "w")

nets = ['SiouxFalls','Anaheim','Barcelona']
#nets = ['Anaheim','Barcelona']
#runs = ['BPC_singleTree','BPC_nestedTree','BPC_twoPhase','FS_NETS','Leblanc']
runs = ['BPC_singleTree','BPC_nestedTree','BPC_twoPhase']
bprop = 0.5
inflate_trips = {'SiouxFalls':1,'Anaheim':4,'Barcelona':2}

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.2f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().CPLEX_threads))
f.write("Budget/Total cost: %.2f\n" % (bprop))
f.write("Trips inflation: %s\n" % (inflate_trips))
f.write('%s\n' % runs)

header = 'Instance '
header += '& UB & Gap (%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE '
header += '& UB & Gap (%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE '
header += '& UB & Gap (%) & Time (s) & RMP (s) & Prc (s) & TAP (s) & SO & UE '
#header += '& UB & Gap (%) & Time (s) & MILP (s) & TAP (s) & SO & UE '
#header += '& UB & Gap (%) & Time (s) & TAP (s) & SO & UE '
header += '\n'
f.write(header)
f.flush()

for net in nets:
    print(net)
    
    if net == 'SiouxFalls':
        ins0 = 'SF_DNDP'
        
    elif net == 'Anaheim':
        ins0 = 'A_DNDP'
        
    elif net == 'Barcelona':
        ins0 = 'B_DNDP'
        
    for nA2 in ['10','20']:
        
        for ID in range(1,3):
            
            ins = ins0+'_'+nA2+'_'+str(ID)
            print('starting to run',ins)
            t0_ins = time.time()
            
            f.write('%s ' % (ins))
            
            for run in runs:
                print(run)
                
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,1e-3,inflate_trips[net])
                #network = Network.Network(net,ins,bprop,1e-0,1e-3)
                
                if run == 'BPC_singleTree':
                    bpc = BPC_singleTree.BPC_singleTree(network)
                    bpc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d' % (bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    
                elif run == 'BPC_twoPhase':
                    bpc = BPC_twoPhase.BPC_twoPhase(network)
                    bpc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d' % (bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                
                elif run == 'BPC_nestedTree':
                    bpc = BPC_nestedTree.BPC_nestedTree(network)
                    bpc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d' % (bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                    
                elif run == 'FS_NETS':
                    fs_nets = FS_NETS.FS_NETS(network)
                    fs_nets.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %d & %d' % (fs_nets.UB,100*fs_nets.gap,fs_nets.rt,fs_nets.rt_MILP,fs_nets.rt_TAP,fs_nets.nSO,fs_nets.nUE))
                
                elif run == 'Leblanc':
                    leblanc = Leblanc.Leblanc(network)
                    leblanc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %d & %d' % (leblanc.UB,100*leblanc.gap,leblanc.rt,leblanc.rt_TAP,leblanc.nSO,leblanc.nUE))
            
            f.write('\\\\\n')
            f.flush()
            
            total_ins_hours = (time.time() - t0_ins)/3600
            print('instance runs completed in %.2f hours' % total_ins_hours)
            
f.close()

total_hours = (time.time() - t0_exp)/3600
print('total runtime of %.2f hours' % total_hours)
