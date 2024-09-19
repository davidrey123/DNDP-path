import llist
from llist import sllist
from src import Node
from collections import deque
import heapq

class CGbush:

    def __init__(self, r, network):
        self.linkflows = {}
        self.link_RC = {}
        self.origin = r
        self.network = network
        
        self.newlinks = []
        
        sptree = network.getSPTree(self.origin)
        for n in network.nodes:
            if n != self.origin:
                self.linkflows[sptree[n]] = 0
                self.link_RC[sptree[n]] = 0
                
                
        self.topologicalSort()
    
    def addLink(self, a):
    
        self.newlinks.append(a)
        
        
        
        
    def processNewLinks(self):
    
        removed = []
        
        for a in self.newlinks:
            
            # check if adding it creates cycle
            removed.extend(self.removeCycle(a.end, a.start))
                
            # add link to bush
            #print("\t\tadding", self.origin, a)
            self.linkflows[a] = 0
            self.link_RC[a] = 0
        
        if len(removed) > 0:
            self.topologicalSort()
            
        self.newlinks = []
        
        
        
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
                # trace
                trace = []
                curr = start
                while curr != end:
                    trace.append(curr.pred)
                    curr = curr.pred.end
                    
                # remove link from trace
                # prioritize links with 0 flow and high reduced cost
                
                best = None
                best_val = -1
                best_zero = None
                
                
                for a in trace:
                    
                    if self.linkflows[a] == 0:
                        val = a.dual
                        
                        if val > best_val:
                            best_val = val
                            best_zero = a
                            
                    elif best_zero is None:
                        val = a.dual
                        
                        if val > best_val:
                            best_val = val
                            best = a
                    '''
                    val = self.link_RC[a]
                        
                    if val > best_val:
                        best_val = val
                        best = a
                    '''
                
                if best_zero is not None:
                    #print("\t\tdeleting", self.origin, best_zero, self.linkflows[best_zero], self.link_RC[best_zero])
                    del self.linkflows[best_zero]
                    del self.link_RC[best_zero]
                    removed.append(best_zero)
                else:
                    #print("\t\tdeleting", self.origin, best, self.linkflows[best], self.link_RC[best])
                    del self.linkflows[best]
                    del self.link_RC[best]
                    removed.append(best)
                
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