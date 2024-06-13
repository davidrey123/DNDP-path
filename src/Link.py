from src import Params

class Link:

    # construct this Link with the given parameters
    def __init__(self, id, start, end, t_ff, C, alpha, beta, cost):
        self.id = id
        self.start = start
        self.end = end
        self.t_ff = t_ff
        self.C = C
        self.alpha = alpha
        self.beta = beta
        self.x = 0
        self.y = 1
        self.cost = cost # for DNDP
        
        self.visit_order = -1
        
        if start is not None:
            start.addOutgoingLink(self)
            
        if end is not None:
            end.addIncomingLink(self)
            
        self.xstar = 0
        self.dual = 0 # for CG

    def setFlow(self, x):
        self.x = x
    
    def __repr__(self):
        return str(self)

    def getTravelTime(self, x, type):
        if self.y == 0 and type != 'RC':
            return Params.INFTY
            
        if type == 'UE':
            output = self.t_ff * (1 + self.alpha * pow(x / self.C, self.beta))
        
        elif type == 'SO':
            output = self.t_ff * (1 + self.alpha * pow(x / self.C, self.beta))
            output += x * self.t_ff * self.alpha * self.beta * pow(x / self.C, self.beta-1) / self.C
            
        elif type == 'RC':
            output = -self.dual
            
        else:
            raise Exception("wrong type "+str(type))

        return output

    def getDerivativeTravelTime(self, x):
        
        if self.y == 0:
            return Params.INFTY
            
        return self.t_ff * self.alpha * self.beta * pow(x / self.C, self.beta-1) / self.C   

    def getCapacity(self):
        return self.C
    
    def getFlow(self):
        return self.x
        
    def __str__(self):
        return "(" + str(self.start.getId()) + ", " + str(self.end.getId()) + ")"
        
    def addXstar(self, flow):
        self.xstar += flow   
    
    def calculateNewX(self, stepsize):        
        self.x = (1 - stepsize) * self.x + stepsize * self.xstar
        self.xstar = 0
        
    def hasHighReducedCost(self, type, percent):
        reducedCost = self.end.cost - self.start.cost
        tt = self.getTravelTime(self.x, type)        
        return tt - reducedCost > tt*percent
 
    def getReducedCost(self, type):
        reducedCost = self.end.cost - self.start.cost
        tt = self.getTravelTime(self.x, type)
        return tt - reducedCost
