class Path:
    # construct this Path; it contains a set of links representing the links in this path
    next_id = 0
    
    def __init__(self):
        self.links = set()
        self.id = self.next_id
        self.next_id += 1
        self.cost = 0
        
        
    def equals(self, rhs):
    	if len(self.links) == len(rhs.links):
    		for a in self.links:
    			if a not in rhs.links:
    				return False
    		return True
    	else:
    		return False
        
    def __hash__(self):
        output = 0
        
        counter = 1
        
        for a in self.links:
        	output += counter * a.__hash__()
        	counter += 1
        
        return output
        
    def __str__(self):
        return str(self.links)     
        
       
        
    def size(self):
        return len(self.links)
   
    def add(self, a):
        self.links.add(a)
    
    def addHstar(self, h):
        for a in self.links:
            a.addXstar(h)
            
    def getTravelTime(self, type):
        output = 0
        for a in self.links:
            output += a.getTravelTime(a.x, type)
        return output
        
    def getTT(self, x, y):
        output = 0
        
        for a in self.links:
            y_e = 0
            
            if a in y:
                y_e = y[a]
            output += a.getTravelTimeC(a.x, y_e, "UE")
        return output
        
    