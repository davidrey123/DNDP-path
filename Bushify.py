import llist
from llist import sllist
from src import Node
from src import Path
from collections import deque
import heapq

class Bushify:

    def __init__(self, r, network, rmp):
        self.linkflows = {}
        self.origin = r
        self.network = network
        self.paths = {s:[] for s in network.zones}
        self.dem_cons = {}
        self.sorted = []
        self.dem_duals = {}
      
        self.network.dijkstras(self.origin,'UE')

        for s in self.network.zones:

            if self.origin.getDemand(s) > 0:

                p = self.network.trace(self.origin, s)
                
                self.addPath(s, p, rmp)
                
                self.dem_cons[s] = rmp.add_constraint(sum(rmp.h[p] for p in self.paths[s]) >= self.origin.getDemand(s), 'dem_%d_%d' % (r.id,s.id)) 
                
                
        self.topologicalSort()
    
    def addPath(self, s, path, rmp):
        rmp.h[path] = rmp.continuous_var(lb=0,ub=self.origin.getDemand(s))
                
        self.paths[s].append(path)
        
        if s in self.dem_cons: # remove corner case of initialization
            self.dem_cons[s].lhs.add_term(rmp.h[path], 1)

        for a in path.links:
            if a not in self.linkflows:
                self.linkflows[a] = 0
            a.link_cons.lhs.add_term(rmp.h[path], -1)
            
    def addPathWithCycle(self, s, path, rmp):
        
        # add TD to each path

        added = set()
        
        if self.origin.id == 1:
            print("newpath", s, path, self.origin.totaldemand)
            
        # look for cycles created by each link
        hasCycle = False
        for a in path.links:
            #print("trying to add link ", a)
            
            #if a not in added:
            flowremoved = self.removeCycle(a.end, a.start)
            if flowremoved > 0:
                hasCycle = True
                #added.add(a)
                
            
            if a in self.linkflows:
                self.linkflows[a] += self.origin.totaldemand - flowremoved
            else:
                self.linkflows[a] = self.origin.totaldemand - flowremoved
            
            print("\t\t", a, flowremoved, self.linkflows[a])
            self.topologicalSort()
        
        self.addPath(s, path, rmp)
        
        return hasCycle
                    
    def regeneratePaths(self, rmp, rem_dem, newpaths):

        #if self.origin.id == 1:
        #    print("\tstarting regenerate", self.origin, sum(self.origin.getDemand(s) for s in self.network.zones))
        new = 0
        self.removeAllPaths(rmp)
    
        origin_out = self.origin.getBushOutgoing(self)
        
        outflow = sum(self.linkflows[a] for a in origin_out)
        
        
        for s in newpaths.keys():
            new += 1
            
            maxflow = 1e15
            for a in newpaths[s].links:
                maxflow = min(self.linkflows[a], maxflow)

            for a in newpaths[s].links:
                self.linkflows[a] -= maxflow
                self.addPath(s, newpaths[s], rmp)

        while outflow > 0.000001:
            removed_flow = 0
            
            if self.origin.id == 1:
                print("\tstart regenerate loop", outflow, self.origin, sum(rem_dem[s] for s in rem_dem.keys()), sum(self.origin.getDemand(s) for s in self.network.zones))
                for a in self.linkflows.keys():
                    print("\t\t", a, self.linkflows[a])

            #    print("rem dem", sum(rem_dem[s] for s in rem_dem.keys()))
            #    for s in rem_dem.keys():
            #        print("\t", s, rem_dem[s])
                
            # max used tree
            self.maxUsedTree()
            
            if self.origin.id == 1:
                print("maxtree check", self.sorted)
                for n in self.sorted:
                    print("\t", n, n.cost, n.pred)
            
            # or s in destination with positive remaining inflow, trace path to origin
            
            for s in self.network.zones:
                if self.origin != s and self.origin.getDemand(s) > 0:
                
                    cutoff = 0.000001
                    
                    if s.cost > 0:
                        maxflow = min(rem_dem[s], s.cost)
                        curr = s
                        trace = Path.Path()
                        
                        if maxflow > 0.00001:
                            while curr != self.origin:
                                pred = curr.pred
                                if pred is None:
                                    maxflow = 0
                                    break
                                if self.linkflows[pred] > cutoff:
                                    maxflow = min(maxflow, self.linkflows[pred])
                                    trace.links.add(pred)
                                    curr = pred.start
                                else:
                                    maxflow = 0
                                    break
                        
                        if maxflow > cutoff:
                            new += 1
                            for a in trace.links:
                                self.linkflows[a] -= maxflow
                            self.addPath(s, trace, rmp)
                            
                            removed_flow += maxflow
                        
                            if self.origin.id == 1:
                                
                                print("\tremoved", maxflow, s, self.origin.getDemand(s), rem_dem[s], trace.links,  sum(rem_dem[s] for s in rem_dem.keys()), sum(self.linkflows[a] for a in origin_out))
                                for a in self.linkflows.keys():
                                    print("\t\t", a, self.linkflows[a])
                                    
                            rem_dem[s] -= maxflow
                        #elif self.origin.id == 1:
                        #    print("\tcannot remove", maxflow, s, rem_dem[s], trace.links)
                            
                            
                        
            outflow = sum(self.linkflows[a] for a in origin_out)
            
            if outflow > 0.001 and removed_flow < 0.00001:
                print("regenerate failed", self.origin.id)

                exit()
            # remove flow and store path if flow is large enough. Maybe 10% of dest inflow?
        
        if sum(rem_dem[s] for s in rem_dem.keys()) > 0.001:
            print("end regenerate", self.origin.id, sum(rem_dem[s] for s in rem_dem.keys()))
            for a in self.linkflows.keys():
                print("\t\t", a, round(self.linkflows[a], 3))
                
            exit()
            
        return new


    def maxUsedTree(self):
        for u in self.sorted:
            u.cost = 0
            u.pred = None
            
        self.origin.cost = 1e15
            
        for u in self.sorted:
            for uv in u.getBushOutgoing(self):
                temp = min(self.linkflows[uv], u.cost)
                if temp > uv.end.cost:
                    uv.end.cost = temp
                    uv.end.pred = uv
        
        
        

        
        
        
        
    def removeAllPaths(self, rmp):
        for s in self.network.zones:
            for p in self.paths[s]:
                self.removePath(s, p, rmp)
            self.paths[s] = []
    
    def removePath(self, s, path, rmp):
        # deleting a variable is hard in docplex, so set the value to 0 to remove it in presolve
        rmp.add_constraint(rmp.h[path] == 0)
        #self.dem_cons[s].lhs.remove_term(rmp.h[path])
        #for a in path.links:
        #    a.link_cons.rhs.remove_term(rmp.h[path])
        
               
    def pricing(self, rmp):
    
        self.linkflows = {}
        

        # rebuild link flows
        for s in self.network.zones:
        
            #if self.origin.id==2:
            #    print("check dem", s, self.origin.getDemand(s), sum(rmp.h[p].solution_value for p in self.paths[s]))
                
            for p in self.paths[s]:
                pathflow = rmp.h[p].solution_value
                
                #if self.origin.id == 2:
                #    print("\t", pathflow, p)
                for a in p.links:
                    if a in self.linkflows:
                        self.linkflows[a] += pathflow
                    else:
                        self.linkflows[a] = pathflow
        
        
        '''
        # reduced cost pricing on bush
        minrc_b, new_b = self.minRCPath(rmp)
        
        if new_b > 0:
            return minrc_b, new_b
        '''
        
        new = 0
        minrc = 1e15
        
        rem_dem = {}
        
        for s in self.network.zones:
            rem_dem[s] = self.origin.getDemand(s)
        
        # RC dijkstras
        self.network.dijkstras(self.origin,'RC')

        
        # attempt to add paths. If cycle is detected, add links from new paths with max flow and remove cycles for added links. Then regenerate paths from link flows
        cycles = False
        
        newpaths = {}
        
        for s in self.network.zones:

            if self.origin.getDemand(s) > 0:

                rc = - self.dem_duals[s] + s.cost                    

                if rc <= - 0.0001:
                    p = self.network.trace(self.origin, s)
                    
                    
                    newpaths[s] = p
                        
                    new += 1
                    
                if rc < minrc:
                    minrc = rc
        
        #if self.origin.id == 2:
        #    for s in newpaths.keys():
        #        print("foundpath", s, newpaths[s])
                
        for s in newpaths.keys():          
            if self.addPathWithCycle(s, newpaths[s], rmp):
                cycles = True
            rem_dem[s] += self.origin.totaldemand
            
        
                  
        
        # I use topological order to do shortest path on bush. If I change the bush, I need to sort again
        self.topologicalSort()
        
        # if any cycles were removed, I need to regenerate paths because some of the paths have cycles
        #if cycles:
        #    new = self.regeneratePaths(rmp, rem_dem)
        
        
        new = self.regeneratePaths(rmp, rem_dem, newpaths)
        
        #print("\t", self.origin.id)
        
        #for s in self.network.zones:
        #    if self.origin.getDemand(s) > 0:
        #        print("\t\t", s, len(self.paths[s]), self.origin.getDemand(s))

        return minrc, new
        
    
    def minRCPath(self, rmp):
    
        minrc = 1e15
        for n in self.network.nodes:
            n.pred = None
            n.cost = 1e15
        
        self.origin.cost = 0
        
        new = 0
        
        for u in self.sorted:
            for uv in u.getBushOutgoing(self):
                temp = uv.dual + u.cost
                if temp < uv.end.cost:
                    uv.end.cost = temp
                    uv.end.pred = uv
    
        for s in self.network.zones:
            if self.origin.getDemand(s) > 0:
                s_price = s.cost - self.dem_duals[s]
                minrc = min(minrc, s_price)
                
                if s_price < 0:
                    self.addPath(s, self.tracePath(self.origin, s), rmp)
                    new += 1
                    
        return minrc, new
        
        
         
    def normal_pricing(self, rmp):
        self.network.dijkstras(self.origin,'RC')

        new = 0
        minrc = 1e15
        
        for s in self.network.zones:

            if self.origin.getDemand(s) > 0:

                rc = - self.dem_duals[s] + s.cost                    

                if rc <= - 0.0001:
                    p = self.network.trace(self.origin, s)
                    self.paths[s].append(p)
                    new += 1
                    
                    self.addPath(s, p, rmp)
                    

                if rc < minrc:
                    minrc = rc
        return minrc, new
        
    def removeCycle(self, start, end): # remove the connection(s) from start to end
        # find path from start to end, if 1 exists
        # this may not be the most efficient, but let's try it for now.
        # BFS or DFS from end towards start. Probably DFS to find path faster
        # then remove links to break the connection and repeat while a connection exists
        
        #print(start, self.top_order[start], end, self.top_order[end])
        
        
        #if start.top_order < end.top_order:
        #    return False
        
        foundStart = True
        
        #removed = []
        
        total_flow = 0
        
        while foundStart:
            for n in self.network.nodes:
                n.visited = False
                n.pred = None

            end.visited = True



            unvisited = deque()
            unvisited.append(end)

            foundStart = False

            while len(unvisited) > 0:
                u = unvisited.pop()

                for a in u.getBushIncoming(self):
                    if a.start.visited == False:
                        a.start.visited = True
                        a.start.pred = a

                        if a.start == start:
                            foundStart = True
                            break

                        unvisited.append(a.start)

                if foundStart:
                    break

            if foundStart:

                #print("cycle search", a.start, a.end)
                #for n in self.network.nodes:
                #    print("\t", n, n.pred)
                #print("--")
                #for a in self.linkflows.keys():
                #    print("\t", a, self.linkflows[a])
                    
                    
                # remove link from trace
                # prioritize links with 0 flow and high reduced cost
                
                cycleflow = self.linkflows[a]
                
                trace = []
                
                curr = start
                while curr != end:
                    #print("\t", curr, start, end, curr.pred)
                    cycleflow = min(cycleflow, self.linkflows[curr.pred])
                    trace.append(curr.pred)
                    curr = curr.pred.end
                    
                #print("cycle found", trace, cycleflow, start.pred)
                    
                
                total_flow += cycleflow
                
                curr = start
                while curr != end:
                    a = curr.pred
                    flow = self.linkflows[a]
                    flow -= cycleflow
                    if flow > 0.000001:
                        self.linkflows[a] = flow
                    else:
                        del self.linkflows[a]
                        cycleFound = True
                    curr = a.end
                    
                
            else:
                break
            
        return total_flow
        
    def topologicalSort(self):
        
        for n in self.network.nodes:
            n.in_degree = len(n.getBushIncoming(self))
            n.top_order = -1

        # Queue to store vertices with indegree 0
        q = deque()
        for n in self.network.nodes:
            if n.in_degree == 0:
                q.append(n)
        self.sorted = []
        
        idx = 0
        while q:
            i = q.popleft()
            self.sorted.append(i)
            i.top_order = idx
            idx += 1
            # Decrease indegree of adjacent vertices as the current node is in topological order
            for ij in i.getBushOutgoing(self):
                j = ij.end
                j.in_degree -= 1
                # If indegree becomes 0, push it to the queue
                if j.in_degree == 0:
                    q.append(j)
                    
        if len(self.sorted) != len(self.network.nodes):
            print("Graph contains cycle!", self.origin.id, self.sorted)
            for n in self.network.nodes:
                if n.top_order < 0:
                    print(n)
            print("--")
            for a in self.linkflows.keys():
                print("\t", a, self.linkflows[a])
            exit()
        
    def testTopologicalSort(self):
    
        #self.topologicalSort()
        
        for l in self.linkflows:
            if self.contains(l) and l.start.top_order > l.end.top_order:
                exit()
                return False

        return True
        
    def contains(self, a):
        return a in self.linkflows
        
    def tracePath(self, i, j):
        output = []
        
        curr = j
        #print("tracing", i, j)
        
        while curr != i:
            if curr.pred == None:
                return None
            output.append(curr.pred)
            #print(curr, curr.pred, i, j, self.origin, output)
            curr = curr.pred.start
        
        return output
    
    def tracePath2(self, i, j):
        output = []
        
        curr = j
        
        while curr != i:
            if curr.pred2 == None:
                return None
            #print(str(curr)+" "+str(curr.pred2)+" "+str(i)+" "+str(j))
            output.append(curr.pred2)
            curr = curr.pred2.start
        
        return output