class Path:
    # construct this Path; it contains a set of links representing the links in this path
    def __init__(self):
        self.links = set()
        
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