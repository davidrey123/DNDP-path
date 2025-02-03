import llist
from llist import sllist
from src import Node
from collections import deque
import heapq

class CGbush:

    def __init__(self, r, network):
        self.linkflows = {}
        self.origin = r
        self.network = network
        
        
        sptree = network.getSPTree(self.origin)
        for n in network.nodes:
            if n in sptree and n != self.origin:
                self.linkflows[sptree[n]] = 0
                
                
        #self.topologicalSort()
    

        
    def addLinks(self, newlinks):
    
        removed = []
        added = set()
        
        for a in newlinks:
            
            # check if adding it creates cycle
            self.linkflows[a] = self.origin.totaldemand
            
            if not a in added:
                removed.extend(self.removeCycle(a.end, a.start))
                added.add(a)
                
            # add link to bush
            #print("\t\tadding", self.origin, a)
            
        
        
        #self.topologicalSort()

        return removed
        
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
                        removed.append(a)
                        del self.linkflows[a]
                    elif not rem and a.dual > best_rc:
                        best_rc = a.dual
                        best = a
                    curr = a.end
                    
                if not rem:
                    removed.append(best)
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
