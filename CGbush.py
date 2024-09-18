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
    
    def addLink(self, a):
    
        self.newlinks.append(a)
        
        self.linkflows[a] = 0
        self.link_RC[a] = 0
        
        
    def processNewLinks(self):
        self.newlinks = []
        
    def hasLink(self, a):
        return a in self.linkflows