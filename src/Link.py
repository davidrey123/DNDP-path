from src import Params

class Link:

    # construct this Link with the given parameters
    def __init__(self, id, start, end, t_ff, C, alpha, beta, cost):
        self.start = start
        self.id = id
        self.end = end
        self.t_ff = t_ff
        self.C = C
        self.alpha = alpha
        self.beta = beta
        self.x = 0
        self.y = 1
        self.cost = cost # for DNDP
        
        if start is not None:
            start.addOutgoingLink(self)
            
        if end is not None:
            end.addIncomingLink(self)
            
        self.xstar = 0
        self.lbdcost = 0
        
    def setlbdCost(self, lbdcost):
        self.lbdcost = lbdcost

    # updates the flow to the given value
    def setFlow(self, x):
        self.x = x
    
    def __repr__(self):
        return str(self)
        
    


    def getTravelTime(self, x, type):
        if self.y == 0:
            return Params.INFTY
            
        
        output = self.t_ff * (1 + self.alpha * pow(x / self.C, self.beta))
        #print(type)
        
        #t = tff (1 + a * (x/c)^b)
        #t' = tff * a * b * (x/c)^b-1 / c
        
        if type == 'SO':
            output += x * self.t_ff * self.alpha * self.beta * pow(x / self.C, self.beta-1) / self.C
        
        #if type != 'TT':
        #    output += self.lbdcost
        return output
        

    def getCapacity(self):
        return self.C
    
    def getFlow(self):
        return self.x
        
    
        

    def __str__(self):
        return "(" + str(self.start.getId()) + ", " + str(self.end.getId()) + ")"
        

    def addXstar(self, flow):
        self.xstar += flow
        #print(xstar)
        
    
    def calculateNewX(self, stepsize):
        #print(str(self.x)+"\t"+ str(self.xstar))
        
        self.x = (1 - stepsize) * self.x + stepsize * self.xstar
        self.xstar = 0
        #print(f"After recalculating, new x = {self.x}, reset xstar = {self.xstar}")
        #print(self.xstar)
        
    def hasHighReducedCost(self, type, percent):
        reducedCost = self.end.cost - self.start.cost
        tt = self.getTravelTime(self.x, type)
        
        return tt - reducedCost > tt*percent
 