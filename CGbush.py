import llist
from llist import sllist
from src import Node
from collections import deque
import heapq

class CGbush:

    def __init__(self, r, network):
        self.linkflows = {}
        self.link_RC = {}
        self.top_order = {}
        self.origin = r
        self.network = network
        
        self.newlinks = []
        
        sptree = network.getSPTree(self.origin)
        for n in network.nodes:
            if n != self.origin:
                self.linkflows[sptree[n]] = 0
                self.link_RC[sptree[n]] = 0
                
                
        #self.topologicalSort()
    
    def addLink(self, a):
    
        self.newlinks.append(a)
        
        
        
        
    def processNewLinks(self):
    
        for a in self.newlinks:
            
            # check if adding it creates cycle
            self.removeCycle(a.end, a.start)
                
            # add link to bush
            #print("\t\tadding", self.origin, a)
            self.linkflows[a] = 0
            self.link_RC[a] = 0
            
        self.newlinks = []
        
        #self.topologicalSort()
        
    def removeCycle(self, start, end): # remove the connection(s) from start to end
        # this may not be the most efficient, but let's try it for now.
        # BFS or DFS from end towards start. Probably BFS to find shortest path
        # then remove one of the links to break the connection and repeat while cycle exists
        # remove link based on combination of link flow and reduced cost
        
        #print(start, self.top_order[start], end, self.top_order[end])
        
        foundStart = True
        
        
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
                        val = self.link_RC[a]
                        
                        if val > best_val:
                            best_val = val
                            best_zero = a
                            
                    elif best_zero is None:
                        val = self.link_RC[a]
                        
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
                else:
                    #print("\t\tdeleting", self.origin, best, self.linkflows[best], self.link_RC[best])
                    del self.linkflows[best]
                    del self.link_RC[best]
                
            else:
                break
            
        
    def topologicalSort(self):
        
        # Initialize in-degrees and visited flags
        for n in self.network.nodes:
            n.in_degree = len(n.getBushIncoming(self))
            #print(n.in_degree)

            n.visited = False
            n.top_order = -1
        
        # Use a list as a priority queue
        queue = []
        heapq.heappush(queue, self.origin.id)  # Assume each node has a unique node_id
        self.origin.visited = True
        
        self.sorted = []
        idx = 0
        
        while queue:
            #print(queue)
            # Use heappop for consistent smallest element first
            vertex_id = heapq.heappop(queue)
            #print(vertex_id)
            vertex = self.network.findNode(vertex_id)  # You need to be able to fetch nodes by ID
            #print(vertex)
            self.sorted.append(vertex)
            vertex.top_order = idx
            self.top_order[vertex] = vertex.top_order
            idx += 1
            
            # Process outgoing edges
            #with open('result2.txt', 'a') as file, contextlib.redirect_stdout(file):
            for ij in vertex.getBushOutgoing(self):
                #print(f"This is ij{ij}")
                j = ij.end
                #print(f"This is j{j}")
                if not j.visited:
                    j.in_degree -= 1
                    if j.in_degree == 0:
                        heapq.heappush(queue, j.id)
                        j.visited = True   
        
        
    def contains(self, a):
        return a in self.linkflows