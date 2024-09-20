import llist
from llist import sllist
from src import Node
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
        self.duals = {}
      
        self.network.dijkstras(self.origin,'UE')

        for s in self.network.zones:

            if self.origin.getDemand(s) > 0:

                p = self.network.trace(self.origin, s)
                
                self.addPath(r, s, p, rmp)
                
                
        self.topologicalSort()
    
    def addPath(self, r, s, path, rmp):
        rmp.h[path] = rmp.continuous_var(lb=0,ub=self.origin.getDemand(s))
                
        self.paths[s].append(path)

        for a in path.links:
            self.linkflows[a] = 0
            a.link_cons.lhs.add_term(rmp.h[path], -1)
                    
                    
    def pricing(self, rmp):
        self.network.dijkstras(self.origin,'RC')

        new = 0
        minrc = 1e15
        
        for s in self.network.zones:

            if self.origin.getDemand(s) > 0:

                rc = - self.duals[s] + s.cost                    

                if rc <= - 0.0001:
                    p = self.network.trace(self.origin, s)
                    self.paths[s].append(p)
                    new += 1
                    
                    new_hp = rmp.continuous_var(lb=0)
                    rmp.h[p] = new_hp
                    self.dem_cons[s].lhs.add_term(new_hp, 1)
                    
                    #link_cons[a] = rmp.add_constraint(rmp.x[a] - sum(rmp.h[p] for p in getPaths(network, paths) if a in p.links) >= 0, 'link_%d_%d' % (a.start.id,a.end.id))
                    for a in p.links:
                        a.link_cons.lhs.add_term(new_hp, -1)

                if rc < minrc:
                    minrc = rc
        return minrc, new
        
    def removeCycle(self, start, end): # remove the connection(s) from start to end
        # this may not be the most efficient, but let's try it for now.
        # BFS or DFS from end towards start. Probably BFS to find shortest path
        # then remove one of the links to break the connection and repeat while cycle exists
        # remove link based on combination of link flow and reduced cost
        
        #print(start, self.top_order[start], end, self.top_order[end])
        
        foundStart = True
        
        removed = []
        
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

                    
                # remove link from trace
                # prioritize links with 0 flow and high reduced cost
                
                cycleflow = 1e15
                
                curr = start
                while curr != end:
                    cycleflow = min(cycleflow, self.linkflows[curr.pred])
                    curr = curr.pred.end
                    
                # find link with highest reduced cost and 0 flow to remove
                
                best_rc = -1e15
                best = None
                rem = False
                
                curr = start
                while curr != end:
                    a = curr.pred
                    flow = self.linkflows[a]
                    flow -= cycleflow
                    if flow > 0.0001:
                        self.linkflows[a] = flow
                    elif a.dual > 0.0001:
                        rem = True
                        del self.linkflows[a]
                    elif not rem and a.dual > best_rc:
                        best_rc = a.dual
                        best = a
                    curr = a.end
                    
                if not rem:
                    del self.linkflows[best]
                
            else:
                break
            
        return removed
        
    def topologicalSort(self):
        
        for n in self.network.nodes:
            n.in_degree = len(n.getBushIncoming(self))

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
        
    def testTopologicalSort(self):
    
        #self.topologicalSort()
        
        for l in self.linkflows:
            if self.contains(l) and l.start.top_order > l.end.top_order:
                exit()
                return False

        return True
        
    def contains(self, a):
        return a in self.linkflows