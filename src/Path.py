class Path:
    # construct this Path; it contains a set of links representing the links in this path
    next_id = 0
    
    def __init__(self):
        self.links = set()
        self.id = self.next_id
        self.next_id += 1
        
    def __hash__(self):
        return hash(self.id)
        
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