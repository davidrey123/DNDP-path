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

configs = [{'initOAheuristic':'kBestKNP','OAcut_tol':0.05,'solveSO':False,'useValueFunctionCuts':False},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.05,'solveSO':False,'useValueFunctionCuts':False},
           {'initOAheuristic':'LocalSearchY1','OAcut_tol':0.05,'solveSO':False,'useValueFunctionCuts':False},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.05,'solveSO':True,'useValueFunctionCuts':False},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.05,'solveSO':False,'useValueFunctionCuts':True},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.05,'solveSO':True,'useValueFunctionCuts':True},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.1,'solveSO':False,'useValueFunctionCuts':False},
           {'initOAheuristic':'LocalSearchKNP','OAcut_tol':0.01,'solveSO':False,'useValueFunctionCuts':False}]

bprop = 0.5
scal_flow = {'SiouxFalls':1e-3,'EasternMassachusetts':1e-3,'BerlinMitteCenter':1e-3,'Anaheim':1e-3,'Barcelona':1e-3}
inflate_trips = {'SiouxFalls':1,'EasternMassachusetts':4,'BerlinMitteCenter':2,'Anaheim':4,'Barcelona':2}

headers = {'BPC':'& UB & Gap (\\%) & rooNodeLB & Time (s) & RMP (s) & Prc (s) & OA (s) & TAP (s) & nBB & nUE \\\\'}

f.write("Params: BB_timelimit (s): %.1f, BB_tol: %.3f, TAP_tol %.3f, CPLEX_threads: %d\n" % (Params.Params().BB_timelimit,Params.Params().BB_tol,Params.Params().min_gap,Params.Params().CPLEX_threads))
f.write("Budget/Total cost: %.2f\n" % (bprop))
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
            
            for ID in range(1,11):
                
                ins = ins0+'_DNDP_'+nA2+'_'+str(ID)
                insshort = ins0+'\\_'+nA2+'\\_'+str(ID)
                print('running',ins)
             
                #---reset network object
                network = Network.Network(net,ins,bprop,1e-0,scal_flow[net],inflate_trips[net])

                bpc = BPC.BPC(network)
                
                bpc.params.initOAheuristic = config['initOAheuristic']
                bpc.params.OAcut_tol = config['OAcut_tol']
                bpc.params.solveSO = config['solveSO']
                bpc.params.useValueFunctionCuts = config['useValueFunctionCuts']
  
                bpc.BB()
                f.write('%s & %.1f & %.2f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %d & %d \\\\\n' % (insshort,bpc.UB,100*bpc.gap,bpc.rootNodeLB,bpc.rt,bpc.rt_RMP,bpc.rt_pricing,bpc.rt_OA,bpc.rt_TAP,bpc.nBB,bpc.nUE))
                f.flush()
            
    total_config_hours = (time.time() - t0_config)/3600
    print('config completed in %.2f hours' % total_config_hours)
            
f.close()

total_hours = (time.time() - t0_exp)/3600
print('total runtime of %.2f hours' % total_hours)
