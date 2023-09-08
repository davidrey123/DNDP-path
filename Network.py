import Node
import Link
import Path
import Zone


class Network:

    # construct this Network with the name; read files associated with network name
    def __init__(self, name):
        self.nodes = [] 
        self.links = []
        self.zones = []
        self.readNetwork("data/"+name+"/net.txt")
        self.readTrips("data/"+name+"/trips.txt")
        
    # read file "/net.txt"
    def readNetwork(self, netFile):
        # **********
        # Exercise 5(b)
        # ********** 
        
        # **********
        # Exercise 5(c)
        # ********** 
        
        # **********
        # Exercise 5(d)
        # ********** 
        
        firstThruNode = 1
        numZones = 0
        numNodes = 0
        numLinks = 0
        
        file = open(netFile, "r")

        line = ""
        
        while line.strip() != "<END OF METADATA>":
            line = file.readline()
            
            if "<NUMBER OF ZONES>" in line:
            
                numZones = int(line[line.index('>') + 1:].strip());
            
            elif "<NUMBER OF NODES>" in line:
                numNodes = int(line[line.index('>') + 1:].strip())
            elif "<NUMBER OF LINKS>" in line:
                numLinks = int(line[line.index('>') + 1:].strip())
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
        while len(line) == 0:
            line = file.readline().strip()
        
        for i in range(0, numLinks):
            line = file.readline().split()
            
            start = self.nodes[int(line[0]) - 1]
            end = self.nodes[int(line[1]) - 1]
            C = float(line[2])

            t_ff = float(line[4])
            alpha = float(line[5])
            beta = float(line[6])
            
            self.links.append(Link.Link(start, end, t_ff, C, alpha, beta))
            
        file.close()
        
    def readTrips(self, tripsFile):

        # **********
        # Exercise 5(d)
        # ********** 
        
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
        
        while len(lines) < line_idx or idx < len(splitted):

            next = splitted[idx]
            if next == "Origin":
                idx += 1
                r = self.zones[int(splitted[idx]) - 1]
            else:
                s = self.zones[int(splitted[idx]) - 1]
                idx += 2
                next = splitted[idx]
                d = float(next[0:len(next) - 1])
                
                #print("add demand " + str(d) + " " + str(r) + " " + str(s))
                r.addDemand(s, d)
                
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

    def getLinks(self):
        return self.links
    
    def getNodes(self):
        return self.nodes
    
    def getZones(self):
        return self.zones

    # **********
    # Exercise 5(e)
    # ********** 
    # find the node with the given id
    def findNode(self, id):
        if id <= 0 or id > len(self.nodes):
            return None
        return self.nodes[id - 1]

    # find the link with the given start and end nodes
    def findLink(self, i, j):
        if i is None or j is None:
            return None

        for link in i.getOutgoing():
            if link.getEnd() == j:
                return link

        return None

    # find the node with the smallest cost (node.cost)
    def argmin(self, set):
        best = None
        min = float("inf")

        for n in set:
            if n.cost < min:
                min = n.cost
                best = n

        return best

    def dijkstras(self, origin):
        # **********
        # Exercise 6(b)
        # ********** 

        for n in self.nodes:
            n.cost = float("inf")
            n.predecessor = None

        origin.cost = 0.0

        Q = {origin}

        # **********
        # Exercise 6(c)
        # ********** 
        
        while len(Q) > 0:

            u = self.argmin(Q)
            Q.remove(u)

            for uv in u.getOutgoing():
                v = uv.getEnd()

                if u.cost + uv.getTravelTime() < v.cost:
                    v.cost = u.cost + uv.getTravelTime()
                    v.predecessor = u

                    if v.isThruNode():
                        Q.add(v)
            
  
    # **********
    # Exercise 6(d)
    # ********** 

    def trace(self, r, s):
        curr = s

        output = Path.Path()
        
        while curr.getId() != r.getId() and curr is not None:
            ij = self.findLink(curr.predecessor, curr)

            if ij is not None:
                output.addFront(ij)
                curr = curr.predecessor
              
        return output

    # **********
    # Exercise 7
    # ********** 
    # returns the total system travel time
    def getTSTT(self):
        output = 0.0
        for ij in self.links:
            output += ij.getFlow() * ij.getTravelTime()
        return output

    # returns the total system travel time if all demand is on the shortest path
    def getSPTT(self):
        output = 0.0

        for r in self.zones:
            if r.getProductions() > 0:
                self.dijkstras(r)

                for s in self.zones:
                    if r.getDemand(s) > 0:
                        output += r.getDemand(s) * s.cost

        return output

    # returns the total number of trips in the network
    def getTotalTrips(self):
        output = 0.0

        for r in self.zones:
            output += r.getProductions()

        return output

    # returns the average excess cost
    def getAEC(self):
        return (self.getTSTT() - self.getSPTT()) / self.getTotalTrips()

    # **********
    # Exercise 8(a)
    # ********** 
    # find the step size for the given iteration number
    def calculateStepsize(self, iteration):
        return 1.0 / iteration

    # **********
    # Exercise 8(b)
    # ********** 
    # calculate the new X for all links based on the given step size
    def calculateNewX(self, stepsize):
        for ij in self.links:
            ij.calculateNewX(stepsize)

    # **********
    # Exercise 8(c)
    # ********** 
    # calculate the all-or-nothing assignment
    def calculateAON(self):
        for r in self.zones:
            if r.getProductions() > 0:
                self.dijkstras(r)

                for s in self.zones:
                    if r.getDemand(s) > 0:
                        pi_star = self.trace(r, s)
                        pi_star.addHstar(r.getDemand(s))

    # **********
    # Exercise 8(d)
    # ********** 
    def msa(self, max_iteration):

        output = "Iteration\tAEC\n"

        for iteration in range(1, max_iteration + 1):
            self.calculateAON()
            stepsize = self.calculateStepsize(iteration)
            
            self.calculateNewX(stepsize)
            
            output += str(iteration) + "\t" + str(self.getAEC()) + "\n"
        
        return output

