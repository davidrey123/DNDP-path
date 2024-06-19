#---modules
import time
from src import Params
from src import Network
from src import Leblanc
from src import FS_NETS
from src import BC
from src import BPC

t0_exp = time.time()

filename = 'results.txt'
f = open(filename, "w")

nets = ['SiouxFalls','Anaheim']
runs = ['BPC','FS_NETS','Leblanc']
#runs = ['BPC','BC','FS_NETS','Leblanc']
bprop = 0.5

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.2f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().CPLEX_threads))
f.write("Budget/Total cost: %.2f\n" % (bprop))
f.write('%s\n' % runs)

if len(runs) == 4:
    f.write('Instance & UB & Gap (%%) & Time (s) & RMP (s) & Pricing (s) & TAP (s) & SO & UE & UB & Gap (%%) & Time (s) & LP (s) & TAP (s) & SO & UE & UB & Gap (%%) & Time (s) & MILP (s) & TAP (s) & SO & UE & UB & Gap (%%) & Time (s) & TAP (s) & SO & UE\n')
elif len(runs) == 3:
    f.write('Instance & UB & Gap (%%) & Time (s) & RMP (s) & Pricing (s) & TAP (s) & SO & UE & UB & Gap (%%) & Time (s) & MILP (s) & TAP (s) & SO & UE & UB & Gap (%%) & Time (s) & TAP (s) & SO & UE\n')

for net in nets:
    print(net)
    
    if net == 'SiouxFalls':
        ins0 = 'SF_DNDP'
        
    elif net == 'Anaheim':
        ins0 = 'A_DNDP'
        
    for nA2 in ['10','20']:
        
        for ID in range(1,2):
            
            ins = ins0+'_'+nA2+'_'+str(ID)
            print('starting to run',ins)
            t0_ins = time.time()
            
            f.write('%s ' % (ins))
            
            for run in runs:
                print(run)
                
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,1e-3)
                
                if run == 'BPC':
                    bpc = BPC.BPC(network)
                    bpc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %d & %d' % (bpc.UB,100*bpc.gap,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_TAP,bpc.nSO,bpc.nUE))
                
                elif run == 'BC':
                    bc = BC.BC(network)
                    bc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %d & %d' % (bc.UB,100*bc.gap,bc.rt,bc.rt_LP,bc.rt_TAP,bc.nSO,bc.nUE))
                    
                elif run == 'FS_NETS':
                    fs_nets = FS_NETS.FS_NETS(network)
                    fs_nets.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %.1f & %d & %d' % (fs_nets.UB,100*fs_nets.gap,fs_nets.rt,fs_nets.rt_MILP,fs_nets.rt_TAP,fs_nets.nSO,fs_nets.nUE))
                
                elif run == 'Leblanc':
                    leblanc = Leblanc.Leblanc(network)
                    leblanc.BB()
                    f.write('& %.1f & %.2f & %.1f & %.1f & %d & %d' % (leblanc.UB,100*leblanc.gap,leblanc.rt,leblanc.rt_TAP,leblanc.nSO,leblanc.nUE))
            
            f.write('\\\\\n')
            
            total_ins_hours = (time.time() - t0_ins)/3600
            print('instance runs completed in %.2f hours' % total_ins_hours)
            
f.close()

total_hours = (time.time() - t0_exp)/3600
print('total runtime of %.2f hours' % total_hours)
