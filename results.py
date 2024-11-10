from matplotlib import pyplot as plt

filename = 'experiments25.txt'
f = open(filename, "r")
lines25 = f.readlines()
f.close()

filename = 'experiments50.txt'
f = open(filename, "r")
lines50 = f.readlines()
f.close()

filename = 'experiments75.txt'
f = open(filename, "r")
lines75 = f.readlines()
f.close()

nIns = 60

headers = {'BPC':'& UB & Gap (\\%) & rootNodeLB & Time (s) & RMP (s) & Prc (s) & OA (s) & TAP (s) & rootNode (s) & nBB & nUE \\\\',
           'BC':'& UB & Gap (\\%) & rootNodeLB & Time (s) & LP (s) & OA (s) & TAP (s) & rootNode (s) & nBB & nUE \\\\',
           'FS_NETS':'& UB & Gap (\\%) & Time (s) & MILP (s) & TAP (s) & rootNode (s) & nBB & nSO & nUE \\\\',
           'Leblanc':'& UB & Gap (\\%) & Time (s) & TAP (s) & rootNode (s) & nBB & nSO & nUE \\\\'}

algs = ['BPC','BC','FS_NETS','Leblanc']
colors = {'BPC':'cyan','BC':'blue','FS_NETS':'orange','Leblanc':'red'}
times = {alg:[] for alg in algs}
times_per_network = {'SF':{alg:[] for alg in algs},'EM':{alg:[] for alg in algs},'BMC':{alg:[] for alg in algs}}

cnt=0
alg_index=0
for line in lines25:
    
    if cnt > 5:
        
        if len(times[algs[alg_index]]) == nIns:
            alg_index += 1
            continue
    
        ls = line.rstrip('\\\\\n').split('&')        
        net = ls[0].split('\_')[0]
        
        if algs[alg_index] == 'BPC' or algs[alg_index] == 'BC':
            times[algs[alg_index]].append(float(ls[4]))            
            times_per_network[net][algs[alg_index]].append(float(ls[4]))
            
        else:            
            times[algs[alg_index]].append(float(ls[3]))
            times_per_network[net][algs[alg_index]].append(float(ls[3]))
            
    cnt+=1
print('processed lines25')

cnt=0
alg_index=0
for line in lines50:
    
    if cnt > 5:
        
        if len(times[algs[alg_index]]) == 2*nIns:
            alg_index += 1
            continue
    
        ls = line.rstrip('\\\\\n').split('&')        
        net = ls[0].split('\_')[0]
        
        if algs[alg_index] == 'BPC' or algs[alg_index] == 'BC':
            times[algs[alg_index]].append(float(ls[4]))            
            times_per_network[net][algs[alg_index]].append(float(ls[4]))
            
        else:            
            times[algs[alg_index]].append(float(ls[3]))
            times_per_network[net][algs[alg_index]].append(float(ls[3]))
            
    cnt+=1
print('processed lines50')

cnt=0
alg_index=0
for line in lines75:
    
    if cnt > 5:
        
        if len(times[algs[alg_index]]) == 3*nIns:
            alg_index += 1
            continue
    
        ls = line.rstrip('\\\\\n').split('&')        
        net = ls[0].split('\_')[0]
        
        if algs[alg_index] == 'BPC' or algs[alg_index] == 'BC':
            times[algs[alg_index]].append(float(ls[4]))            
            times_per_network[net][algs[alg_index]].append(float(ls[4]))
            
        else:            
            times[algs[alg_index]].append(float(ls[3]))
            times_per_network[net][algs[alg_index]].append(float(ls[3]))
            
    cnt+=1    
print('processed lines75')

#---compute performance profiles
timelimit = 3600
timesteps = [t for t in range(timelimit)]
scores = {alg:[0 for t in range(timelimit)] for alg in algs}
scores_per_network = {'SF':{alg:[0 for t in range(timelimit)] for alg in algs},'EM':{alg:[0 for t in range(timelimit)] for alg in algs},'BMC':{alg:[0 for t in range(timelimit)] for alg in algs}}

step = 1
for alg in algs:
    for t in range(timelimit):
        for time in times[alg]:
            if time <= t:
                scores[alg][t] += 1
                
        for net in times_per_network:
            for time in times_per_network[net][alg]:
                if time <= t:
                    scores_per_network[net][alg][t] += 1
                
for alg in algs:
    for t in range(timelimit):     
        scores[alg][t] = 100*scores[alg][t]/(3*nIns)
        
        for net in scores_per_network:
            scores_per_network[net][alg][t] = 100*scores_per_network[net][alg][t]/60

print(len(scores_per_network['SF']['BPC']))

fig = plt.figure(figsize=(10,8))

for alg in algs:
    x = timesteps
    y = scores[alg]
    plt.plot(x, y, color=colors[alg], lw=3, alpha=0.6, label=alg)
plt.xlabel('Runtime (s)')
plt.ylabel('Percentage of instances solved out of 180 instances')
plt.title('Performance profiles of DNDP algorithms')
plt.xlim(0, timelimit)
plt.ylim(0, 101)
plt.legend()
plt.savefig('performanceProfiles.pdf')
plt.show()

for net in scores_per_network:
        
    fig = plt.figure(figsize=(10,8))
    
    for alg in algs:
        x = timesteps
        y = scores_per_network[net][alg]
        plt.plot(x, y, color=colors[alg], lw=3, alpha=0.6, label=alg)
                
    plt.xlabel('Runtime (s)')
    plt.ylabel('Percentage of instances solved out of 60 instances')
    plt.title('Performance profiles of DNDP algorithms on '+net+' network')
    plt.xlim(0, timelimit)
    plt.ylim(0, 101)
    plt.legend()
    plt.savefig('performanceProfiles_'+net+'.pdf')
    plt.show()