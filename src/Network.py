from src import Node
from src import Link
from src import Path
from src import Zone
from src import Bush
from src import Params
from src import PASList
from src import Heap
from src import HeapMP
import time
import threading
from multiprocessing import Process, Manager, Array, Value




                
def dijkstras_mp(origin, num_nodes, link_costs, node_out, SP_tol):



    node_costs = {n: Params.INFTY for n in range(1, num_nodes+1)}
    node_pred = {n: None for n in range(1, num_nodes+1)}


    node_costs[origin] = 0.0

    Q = HeapMP.HeapMP(node_costs = node_costs)
    Q.insert(origin)

    #if type == 'RC':
    #    print('origin',origin)

    while Q.size() > 0:

        u = Q.removeMin()

        #if type == 'RC':
        #    print('u',u)

        
        for i in range(0, len(node_out[u])):
            v = node_out[u][i]
            u_cost = node_costs[u]
            v_cost = node_costs[v]
            tt = link_costs[(u,v)]

            #if u.cost + tt < v.cost:
            if u_cost + tt < v_cost and v_cost - u_cost - tt >= SP_tol:
                node_costs[v] = u_cost + tt
                node_pred[v] = u

                #if type == 'RC':
                #    print('v',v,v.pred,v.cost)

                #if v.isThruNode():
                Q.insert(v)



    return node_costs, node_pred
                        
def trace_mp(node_pred, r, s):
    curr = s

    output = []


    while curr != r and curr is not None:
        i = node_pred[curr]

        if i is not None:
            output.append(i)
            curr = node_pred[curr]

    #print('trace',r,s,output)

    return output
    
    
def apsp_threaded_mp_help(trips, num_nodes, link_costs, node_out, SP_tol, get_paths, start_idx, end_idx, output_costs, output_paths):
    
    
    for r in range(start_idx, end_idx):
        node_costs, node_pred = dijkstras_mp(r, num_nodes, link_costs, node_out, SP_tol)

        output_costs_r = dict()
        output_paths_r = dict()

        for s in range(1, num_nodes+1):
            if trips[(r,s)] > 0:
                output_costs_r[s] = node_costs[s]

                if get_paths:
                    path_rs = trace_mp(node_pred, r, s)
                    output_paths_r[s] = path_rs

                    #path_rs.cost = node_costs[s]


        output_costs[r] = output_costs_r

        if get_paths:
            output_paths[r] = output_paths_r
            
            
        

class Network:

    # construct this Network with the name; read files associated with network name
    def __init__(self,name,ins,B_prop,scal_time,scal_flow,inflate_trips):
        self.nodes = [] 
        self.links = []
        self.zones = []
        self.origins = []
        
        self.links2 = []
        self.type = 'UE'
        self.TD = 0
        self.TC = 0 # total cost
        self.params = Params.Params()
        
        self.ins = ins        
        
        self.dijkstras_time = 0
        self.allPAS = PASList.PASList()        
        
        self.inf = 1e+9
        self.tol = 1e-2
        
        self.node_out = None
        self.trips = None
        
        self.manager = Manager()
        
        if len(ins) == 0:
            ins = "net"
            
        self.readNetwork("data/"+name+"/"+ins+".txt",scal_time,scal_flow)
        self.readTrips("data/"+name+"/trips.txt",scal_time,scal_flow,inflate_trips)
        
        
        self.links1 = []
        
        self.linkid_map = dict()
        self.nodeid_map = dict()
        
        for a in self.links:
            if a not in self.links2:
                self.links1.append(a)
        
        self.B = self.TC * B_prop # budget        
        
        #print('Total scaled demand %.1f' % self.TD)
        #print('Total cost %.1f - Budget %.1f' % (self.TC, self.B))
        

            
    def setType(self, type):
        self.type = type
        
    # read file "/net.txt"
    def readNetwork(self,netFile,scal_time,scal_flow):
        
        firstThruNode = 1
        numZones = 0
        numNodes = 0
        numLinks = 0
        newLinks = 0
        file = open(netFile, "r")

        line = ""
        
        while line.strip() != "<END OF METADATA>":
            line = file.readline()
            if "<NUMBER OF ZONES>" in line:            
                numZones = int(line[line.index('>') + 1:].strip())            
            elif "<NUMBER OF NODES>" in line:
                numNodes = int(line[line.index('>') + 1:].strip())
            elif "<NUMBER OF LINKS>" in line:
                numLinks = int(line[line.index('>') + 1:].strip())
            elif "<NUMBER OF NEW LINKS>" in line:
                newLinks = int(line[line.index('>') + 1:].strip())
            elif "<FIRST THRU NODE>" in line:
                firstThruNode = int(line[line.index('>') + 1:].strip())

        for i in range(0, numZones):
            self.zones.append(Zone.Zone(i + 1))

        for i in range(0, numNodes):
            if i < numZones:
                self.nodes.append(self.zones[i])
                
                if i + 1 < firstThruNode:
                    self.zones[i].setThruNode(False)

            else:
                self.nodes.append(Node.Node(i + 1))

        line = ""
        id = 0
        while len(line) == 0:
            line = file.readline().strip()

        for i in range(0, numLinks + newLinks):
            line = file.readline().split()
            if len(line) == 0:
                continue
            start = self.nodes[int(line[0]) - 1]
            end = self.nodes[int(line[1]) - 1]
            C = float(line[2]) * scal_flow             

            t_ff = float(line[4]) * scal_time
            alpha = float(line[5])
            beta = float(line[6])
            
            cost = float(line[10])
            
            self.TC += cost
            
            link = Link.Link(id, start ,end, t_ff, C, alpha, beta, cost)
            id = id +1
            #print(start,end)
            self.links.append(link)
            #print(self.links)
            
            if i >= numLinks:
                self.links2.append(link)
            
        file.close()
        
        self.nodeid_map = {n.id:n for n in self.nodes}
        self.linkid_map = {(ij.start.id, ij.end.id):ij for ij in self.links}
        
        
        


    def readTrips(self,tripsFile,scal_time,scal_flow,inflate_trips):
        
        file = open(tripsFile, "r")
        
        lines = file.readlines()
        
        line_idx = 0
        
        while lines[line_idx].strip() != "<END OF METADATA>":
            line_idx += 1
            
        line_idx += 1
        
        while lines[line_idx].strip() == "":
            line_idx += 1
            
        r = None
        
        idx = 0
        
        splitted = lines[line_idx].split()
        #print(splitted)
        
        while len(lines) < line_idx or idx < len(splitted):

            next = splitted[idx]

            if next == "Origin":
                
                idx += 1
                r = self.zones[int(splitted[idx]) - 1]

            else:
                s = self.zones[int(splitted[idx]) - 1]

                #print(s)
                idx += 2
                next = splitted[idx]
                d = float(next[0:len(next) - 1]) * scal_flow
                d = d * inflate_trips
                
                r.addDemand(s, d)
                self.TD += d

            idx += 1

            if idx >= len(splitted):
                line_idx += 1
                while line_idx < len(lines) and lines[line_idx].strip() == "":
                    line_idx += 1
                    
                if line_idx < len(lines):
                    line = lines[line_idx].strip()
                    splitted = line.split()
                    idx = 0
            
        file.close()
        
        for r in self.zones:
            if r.getProductions() > self.params.flow_epsilon:
                self.origins.append(r)

    def getLinks(self):
        return self.links
    
    def getNodes(self):
        return self.nodes
    
    def getZones(self):
        return self.zones

    # find the node with the given id
    def findNode(self, id):
        if id <= 0 or id > len(self.nodes):
            return None
        return self.nodes[id - 1]

    # find the link with the given start and end nodes
    def findLink(self, i, j):
        if i is None or j is None:
            return None

        for link in i.outgoing:
            if link.end == j:
                return link

        return None
    
    def memoizeLinkCosts(self, type):
        for a in self.links:
            a.saved_tt = a.getTravelTime(a.x, type)
            
            
    
            
    def apsp(self, type, get_paths=True):
    
        
        if self.params.num_threads == 0:
            t1 = time.time()
            self.memoizeLinkCosts(type)


            output_costs = {r: dict() for r in self.origins}
        
            output_paths = None

            if get_paths:
                output_paths = {r: dict() for r in self.origins}


            for r in self.origins:
                node_costs, node_pred = self.dijkstras_sep(r, 'memoized')


                for s in self.zones:
                    if r.getDemand(s) > 0:
                        output_costs[r][s] = node_costs[s]

                        if get_paths:
                            path_rs = self.trace_sep(node_pred, r, s)
                            output_paths[r][s] = path_rs

                            path_rs.cost = node_costs[s]

            t1 = time.time() - t1
        
            self.dijkstras_time += t1    

            return output_costs, output_paths
        else:
            #thread_test()
            return self.apsp_threaded(type, get_paths)
            
        
        
        
    def apsp_threaded(self, type, get_paths=True):
    
        
        t1 = time.time()
    
        self.memoizeLinkCosts(type)
        
        output_costs = {r: dict() for r in self.origins}
        
        output_paths = None
        
        if get_paths:
            output_paths = {r: dict() for r in self.origins}
            
            
        break_points = []
        
        for i in range(0, self.params.num_threads+1):
            break_points.append(round(i*len(self.origins) / self.params.num_threads))
            
        threads = []
            
        for i in range(0, self.params.num_threads):
            thread = threading.Thread(target=self.apsp_threaded_help, args=('memoized', get_paths, break_points[i], break_points[i+1], output_costs, output_paths))
            threads.append(thread)
            
        for i in range(0, len(threads)):
            threads[i].start()

            
        for i in range(0, len(threads)):
            threads[i].join()
            
        
        
        t1 = time.time() - t1
        
        self.dijkstras_time += t1    
        
        
        return output_costs, output_paths
        
    def apsp_threaded_help(self, type, get_paths, start_idx, end_idx, output_costs, output_paths):
    
        
        for idx in range(start_idx, end_idx):
            r = self.origins[idx]
            node_costs, node_pred = self.dijkstras_sep(r, 'memoized')
            
            output_costs_r = dict()
            output_paths_r = dict()
            
            for s in self.zones:
                if r.getDemand(s) > 0:
                    output_costs_r[s] = node_costs[s]
                
                    if get_paths:
                        path_rs = self.trace_sep(node_pred, r, s)
                        output_paths_r[s] = path_rs

                        #path_rs.cost = node_costs[s]
            
            
            output_costs[r] = output_costs_r
            
            if get_paths:
                output_paths[r] = output_paths_r
            
    
    

    def apsp_threaded_mp(self, type, get_paths=True):
    
        self.trips = self.manager.dict()
        
        if self.trips is None:
            for r in self.origins:
                for s in self.nodes:
                    self.trips[(r.id, s.id)] = r.getDemand(s)
                

                
                    
        if len(self.nodeid_map) == 0:
            self.nodeid_map = {n.id:n for n in self.nodes}
            self.linkid_map = {(ij.start.id, ij.end.id):ij for ij in self.links}
        
        
        t1 = time.time()

        self.memoizeLinkCosts(type)

        link_costs = self.manager.dict()


        for a in self.links:
            link_costs[(a.start.id, a.end.id)] = a.getTravelTime(a.x, type)

        self.node_out = self.manager.dict()
        for i in self.nodes:
            self.node_out[i.id] = self.manager.list()
            for ij in i.outgoing:
                self.node_out[i.id].append(ij.end.id)

        

        for r in self.origins:
            for s in self.nodes:
                self.trips[(r.id, s.id)] = r.getDemand(s)

        break_points = []

        for i in range(0, self.params.num_threads+1):
            break_points.append(1+round(i*len(self.origins) / self.params.num_threads))

        threads = []

        temp_costs = dict()
        temp_paths = dict()

        for i in range(0, self.params.num_threads):
            temp_costs[i] = self.manager.dict()
            temp_paths[i] = self.manager.dict()

            thread = Process(target=apsp_threaded_mp_help, args=(self.trips, len(self.nodes), link_costs, self.node_out, get_paths, self.params.SP_tol, break_points[i], break_points[i+1], temp_costs[i], temp_paths[i]))
            threads.append(thread)

        for i in range(0, len(threads)):
            threads[i].start()


        for i in range(0, len(threads)):
            threads[i].join()



        t1 = time.time() - t1

        self.dijkstras_time += t1    

        output_costs = {r: dict() for r in self.origins}
        output_paths = dict()

        #print(self.nodeid_map)


        for i in range(0, len(threads)):
            for r in temp_costs[i].keys():
                for s in temp_costs[i][r].keys():

                    output_costs[self.nodeid_map[r]][self.nodeid_map[s]] = temp_costs[i][r][s]

        output_paths = None

        return output_costs, output_paths
                    
    
        
    def dijkstras_sep(self, origin, type):
        
        
        
        node_costs = {n: Params.INFTY for n in self.nodes}
        node_pred = {n: None for n in self.nodes}
        

        node_costs[origin] = 0.0

        Q = Heap.Heap(node_costs = node_costs, ext_compare = True)
        Q.insert(origin)

        #if type == 'RC':
        #    print('origin',origin)

        while Q.size() > 0:

            u = Q.removeMin()

            #if type == 'RC':
            #    print('u',u)

            for uv in u.outgoing:

                v = uv.end
                u_cost = node_costs[u]
                v_cost = node_costs[v]
                tt = uv.getTravelTime(uv.x, type)

                #if u.cost + tt < v.cost:
                if u_cost + tt < v_cost and v_cost - u_cost - tt >= self.params.SP_tol:
                    node_costs[v] = u_cost + tt
                    node_pred[v] = uv

                    #if type == 'RC':
                    #    print('v',v,v.pred,v.cost)

                    if v.isThruNode():
                        Q.insert(v)
                        
        
                        
        return node_costs, node_pred
                        
    def trace_sep(self, node_pred, r, s):
        curr = s

        output = Path.Path()
        output.r = r
        output.s = s
        
        
        while curr != r and curr is not None:
            ij = node_pred[curr]

            if ij is not None:
                output.add(ij)
                curr = node_pred[curr].start
              
        #print('trace',r,s,output)
              
        return output
            
                    
    def dijkstras(self, origin, type):
        
        #t1 = time.time()
        
        for n in self.nodes:
            n.cost = Params.INFTY
            n.pred = None

        origin.cost = 0.0

        Q = Heap.Heap()
        Q.insert(origin)

        #if type == 'RC':
        #    print('origin',origin)

        while Q.size() > 0:

            u = Q.removeMin()

            #if type == 'RC':
            #    print('u',u)

            for uv in u.outgoing:

                v = uv.end
                tt = uv.getTravelTime(uv.x, type)

                #if u.cost + tt < v.cost:
                if u.cost + tt < v.cost and v.cost - u.cost - tt >= self.params.SP_tol:
                    v.cost = u.cost + tt
                    v.pred = uv

                    #if type == 'RC':
                    #    print('v',v,v.pred,v.cost)

                    if v.isThruNode():
                        Q.insert(v)
        
        #t1 = time.time() - t1
        
        #self.dijkstras_time += t1

    def trace(self, r, s):
        curr = s

        output = Path.Path()
        output.r = r
        output.s = s
        
        while curr != r and curr is not None:
            ij = curr.pred

            if ij is not None:
                output.add(ij)
                curr = curr.pred.start
              
        #print('trace',r,s,output)
              
        return output
        
    def traceTree(self, tree, r, s):
            curr = s

            output = []

            while curr != r and curr is not None:
                ij = tree[curr]

                if ij is not None:
                    output.append(ij)
                    curr = ij.start
            
            return output

    def getSPTree(self, r):
        return self.dijkstras_sep(r, self.type)   
        
    
    # returns the total system travel time
    def getTSTT(self, type):
        output = 0.0
        for ij in self.links:
            #if ij.y == 1 or ij.x > len(self.origins) * self.params.flow_epsilon:
            if ij.y == 1:
                tt = ij.getTravelTime(ij.x, type)
                output += ij.x * tt
                
        return output
    
    def validateLinkFlows(self):
        output = True
        for ij in self.links:
            totbushflow = 0
            
            for r in self.origins:
                totbushflow += r.bush.getFlow(ij)
            
            if abs(ij.x - totbushflow) > self.params.flow_epsilon:
                #print(ij, ij.x, totbushflow, ij.x-totbushflow, ij.getTravelTime(ij.x, self.type))
                output = False
        return output
            
    
    # returns the total system travel time if all demand is on the shortest path
    def getSPTT(self, type):
        output = 0.0

        rs_costs, rs_paths = self.apsp(type, get_paths=False)
        
        for r in self.origins:
            for s in self.zones:
                if r.getDemand(s) > 0:
                    output += rs_costs[r][s] * r.getDemand(s)

        return output

    # returns the total number of trips in the network
    def getTotalTrips(self):
        output = 0.0

        for r in self.origins:
            output += r.getProductions()

        return output

    # returns the average excess cost
    def getAEC(self):
        return (self.getTSTT() - self.getSPTT()) / self.getTotalTrips()
    
    # returns the UE TAP objective function value
    def getBeckmannOFV(self):
        output = 0.0
        for a in self.links:
            if a.y == 1:
                output += a.getPrimitiveTravelTime(a.x)                
                
        return output

    # find the step size for the given iteration number
    def calculateStepsize(self, iteration):
        return 1.0 / iteration
        #print(1.0 / iteration)


    # calculate the new X for all links based on the given step size
    def calculateNewX(self, stepsize):
        for ij in self.links:
            ij.calculateNewX(stepsize)


    # calculate the all-or-nothing assignment
    def calculateAON(self):
        for r in self.origins:
            self.dijkstras(r, self.type)

            for s in self.zones:
                if r.getDemand(s) > 0:
                    pi_star = self.trace(r,s)
                    pi_star.addHstar(r.getDemand(s))
                    
    def setAON(self, type, y):
        self.setY(y)
        
        for r in self.origins:
            self.dijkstras(r, type)

            for s in self.zones:
                if r.getDemand(s) > 0:
                    pi_star = self.trace(r,s)
                    pi_star.addHstar(r.getDemand(s))
                    
        for ij in self.links:
            ij.calculateNewX(1)
            
        return self.getTSTT('UE')
    
    def getFlowMap(self):
        output = {}
        
        for a in self.links:
            output[a] = a.x
        
        return output
        
    def setFlows(self, flowmap):
        for a in self.links:
            a.x = flowmap[a]
        
    def setY(self, y):
    
        if y is not None:
            newlinks = []
            removedlinks = []

            for ij in self.links2:
                if y[ij] != ij.y:
                    if y[ij] == 0:
                        removedlinks.append(ij)
                    else:
                        newlinks.append(ij)
                ij.y = y[ij]

            # delete PAS using removedlinks        
            for r in self.origins:
                if r.bush != None:
                    r.bush.addLinks(newlinks)
                    r.bush.removeLinks(removedlinks)
                


    def msa(self, type, y):
        self.setY(y)
        self.setType(type)        
        
        max_iteration = self.params.tapas_max_iter
        min_gap = self.params.min_gap
        
       
        if self.params.PRINT_TAP_ITER:
            print("Iteration\tTSTT\tSPTT\tgap\tAEC")
        
        
        for iteration in range(1, max_iteration + 1):
            self.calculateAON()
            stepsize = self.calculateStepsize(iteration)
            
            self.calculateNewX(stepsize)
            
            tstt = self.getTSTT(self.type)
            sptt = self.getSPTT(self.type)
            gap = (tstt - sptt)/tstt
            aec = (tstt - sptt)/self.TD
            
            if self.params.PRINT_TAP_ITER:
                print(str(iteration)+"\t"+str(tstt)+"\t"+str(sptt)+"\t"+str(gap)+"\t"+str(aec))
                
            if gap < min_gap:
                break
        
        return self.getTSTT('UE')
        
    def resetTapas(self):
        for r in self.origins:
            r.bush = None
            
        for a in self.links:
            a.x = 0
                
        self.params = Params.Params()
        
        self.allPAS = PASList.PASList()
    
    def tapas(self, type, y):
        
        if not self.params.warmstart:
            self.resetTapas()
            
        self.setY(y)
        self.setType(type)
        
        iter = 1
        max_iter = self.params.tapas_max_iter
        
        if type == 'UE' or type == 'SO':
            min_gap = self.params.min_gap
        elif type == 'SO_OA_cuts':
            min_gap = self.params.min_gap_SO_OA_cuts
        
        #self.params.line_search_gap = pow(10, math.floor(math.log10(self.TD) - 6))
        
        if self.params.PRINT_TAP_ITER:
            print("Iteration\tTSTT\tSPTT\tgap\tAEC")
            
        last_iter_gap = 1
        
        for r in self.origins:
            if r.bush == None:
                r.bush = Bush.Bush(self, r)
        
        while iter <= max_iter:
                        
            #custom_x = {link: link.x for link in self.links}
            #print(custom_x)
            
            self.params.good_pas_flow_mu = 0
            self.params.good_pas_cost_mu = 0
            self.params.good_bush_gap = 0
            self.params.good_pas_cost_epsilon = 0
            
            self.checkBushes()
                
            for r in self.origins:        
                r.bush.branchShifts()
            
                printed = False                              
                for a in r.bush.relevantPAS.forward:
                    for p in r.bush.relevantPAS.forward[a]:
                        if self.params.PRINT_TAPAS_INFO and not printed:
                            print("initial flow shifts", r)
                            printed = True
                        p.flowShift(self.type, self.params)
                        
                        # for every active PAS
             
            if self.params.PRINT_TAPAS_INFO:
                print("general flow shifts")
                               
            modified = False
            for shiftIter in range(0, self.params.tapas_equilibrate_iter):
                # check if it should be eliminated
                self.removePAS(iter)
                # perform flow shift to equilibrate costs
                modified = self.equilibratePAS(iter)
                # redistribute flows between origins by the proportionality condition
                            
                # in the case that no flow shifting occurred, do not try to equilibrate more
                if not modified:
                    break
                            
            tstt = self.getTSTT(type)
            sptt = self.getSPTT(type)
            gap = (tstt - sptt)/tstt
            aec = (tstt - sptt)/self.TD

            #print(iter, sptt)
                            
            if self.params.PRINT_TAP_ITER:
                print(str(iter)+"\t"+str(tstt)+"\t"+str(sptt)+"\t"+str(gap)+"\t"+str(aec))
                
                #printLinkFlows();

            if gap < min_gap:
                break
                
            # there's an issue where PAS are labeled as not cost effective because the difference in cost is small, less than 5% of the reduced cost
            # for low network gaps, this is causing PAS to not flow shift
            # when the gap is low, increase the flow shift sensitivity
            '''
            if (last_iter_gap - gap) / gap < 0.01:
                print('TAPAS check 1')
                
                #for r in self.origins:
                #    r.bush.algBShift()
                
                self.params.bush_gap = max(self.params.bush_gap/10, 1e-6)
                self.params.pas_cost_mu = max(self.params.pas_cost_mu/10, 1e-6)
                self.params.pas_flow_mu = max(self.params.pas_flow_mu/10, 1e-6)
                self.params.flow_epsilon = max(self.params.flow_epsilon/10, 1e-6)
               
                self.params.line_search_gap = max(self.params.line_search_gap/10, self.params.min_line_search_gap)
                    
                if self.params.PRINT_TAPAS_INFO:
                    print("Adjusting parameters due to small gap "+str(self.params.pas_cost_mu)+" "+str(self.params.line_search_gap))
            '''             
            
            if (last_iter_gap - gap) < min_gap:
                
                self.params.line_search_gap = max(self.params.line_search_gap/10, self.params.min_line_search_gap)
                
                if self.params.good_pas_cost_mu <= 10 and self.params.pas_cost_mu > 5e-5:
                    self.params.pas_cost_mu = max(self.params.pas_cost_mu/10, 1e-6)   
                    
                elif self.params.good_pas_flow_mu <= 10 and self.params.pas_flow_mu > 5e-5:
                    self.params.pas_flow_mu = max(self.params.pas_flow_mu/10, 1e-6) 
                   
                elif self.params.good_pas_cost_epsilon == 0 and self.params.pas_cost_epsilon > 5e-6:
                    self.params.pas_cost_epsilon = max(self.params.pas_cost_epsilon/10, 1e-6) 
                    
                elif self.params.bush_gap > 5e-6:
                    self.params.bush_gap = max(self.params.bush_gap/10, 1e-6)
                    self.params.resetPAS()
                
                if self.params.PRINT_TAP_ITER:
                    print('TAPAS gap check', self.params.bush_gap, self.params.pas_cost_mu, self.params.pas_flow_mu, self.params.pas_cost_epsilon)
                    print("\t", self.params.good_bush_gap, self.params.good_pas_cost_mu, self.params.good_pas_flow_mu, self.params.good_pas_cost_epsilon)
                
            last_iter_gap = gap
            iter += 1
        
        #---rounding link flows for numerical stability
        for a in self.links:
            a.x = round(a.x,self.params.rd)
            #print(a.start.id,a.end.id,a.x)
            
        return self.getTSTT('UE')
        
    def findPAS(self, ij, bush):
        
        if not self.allPAS.containsKey(ij):
            return None
        
        #best = None
        #max = self.params.bush_gap
        
        if ij in self.allPAS.backward:
            for p in self.allPAS.backward[ij]:
                '''
                bwdcost = p.getBackwardCost(self.type)
                fwdcost = p.getForwardCost(self.type)
                
                if fwdcost > bwdcost * (1 + self.params.pas_cost_mu):
                    bwdflow = p.maxBackwardBushFlowShift(bush)
                    fwdflow = p.maxForwardBushFlowShift(bush)
                    if fwdflow > max:
                        best = p
                elif bwdcost > fwdcost * (1 + self.params.pas_cost_mu):
                    bwdflow = p.maxBackwardBushFlowShift(bush)
                    fwdflow = p.maxForwardBushFlowShift(bush)
                    if bwdflow > max:
                        best = p
                
                '''
                if p.isEffective(self.type, bush, self.params.pas_cost_mu, bush.origin.getProductions()*self.params.pas_flow_mu, self.params):
                    return p
                        
        if ij in self.allPAS.forward:
            for p in self.allPAS.forward[ij]:
                if p.isEffective(self.type, bush, self.params.pas_cost_mu, bush.origin.getProductions()*self.params.pas_flow_mu, self.params):
                    return p
          
        return None        
        
    def equilibratePAS(self, iter):
        output = False
        
        for a in self.allPAS.forward:
            for p in self.allPAS.forward[a]:
                if p.flowShift(self.type, self.params):
                    output = True
                    p.lastIterFlowShift = iter

        return output        
        
    def removeAPAS(self, p):
        self.allPAS.remove(p)
            
        for r in p.relevant:
            r.bush.relevantPAS.remove(p)
    
    def removePAS(self, iter):
        removed = []
        
        for a in self.allPAS.forward:
            for p in self.allPAS.forward[a]:
                if p.lastIterFlowShift < iter-2:
                    removed.append(p)
        
        for p in removed:
            self.removeAPAS(p)


    def checkBushes(self):
        minPathTrees = dict()
        minPathCosts = dict()

        for r in self.origins:
            minPathCosts[r], minPathTrees[r] = self.getSPTree(r)


        # for every origin
        for r in self.origins:

            minPathTree = minPathTrees[r]
            minPathCost = minPathCosts[r]

            #minPathCost, minPathTree = self.getSPTree(r)

            # remove all cyclic flows and topological sort
            if self.params.PRINT_TAPAS_INFO:
                print("removing cycles", r)

            r.bush.removeCycles()
            # find tree of least cost routes



            if self.params.PRINT_TAPAS_INFO:
                print("checking for PAS", r)

            r.bush.checkPAS(minPathTree, minPathCost)

            if self.params.PRINT_PAS_INFO:
                print("num PAS", r, r.bush.relevantPAS.size())
            # for every link used by the origin which is not part of the tree
                # if there is an existing effective PAS
                    # make sure the origin is listed as relevant
                # else
                    # construct a new PAS    

            # choose a random subset of active PASs
            # shift flow within each chosen PAS