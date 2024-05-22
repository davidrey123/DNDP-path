# Created on : May 22, 2024, 9:54:48 AM
# Author     : michaellevin

class YDict:
    def __init__(self):
        self.tstt_map = {}
        
    def insertSO(self, yvec, tstt):
        self.insertTSTT(yvec, tstt, 1)
            
    def insertUE(self, yvec, tstt):
        self.insertTSTT(yvec, tstt, 2)
    
    def hasSO(self, yvec):
        return self.hasTSTT(yvec, 1)
        
    def hasUE(self, yvec):
        return self.hasTSTT(yvec, 2)
        
    def getSO(self, yvec):
        return self.getTSTT(yvec, 1)
    
    def getUE(self, yvec):
        return self.getTSTT(yvec, 2)
        
        
    
    
    def insertTSTT(self, yvec, tstt, idx):
        hash = self.hashcode(yvec)
        
        ylist = None
        
        if hash in self.tstt_map:
            ylist = self.tstt_map[hash]
        else:
            ylist = []
            self.tstt_map[hash] = ylist
            
        for ytuple in ylist:
            if self.equals(yvec, ytuple[0]):
                ytuple[idx] = tstt
                return
        
        newtuple = [yvec, None, None]
        newtuple[idx] = tstt
        ylist.append(newtuple)
                
    def getTSTT(self, yvec, idx):
        hash = self.hashcode(yvec)
        
        if hash in self.tstt_map:
            ylist = self.tstt_map[hash]

            for ytuple in ylist:
                if self.equals(yvec, ytuple[0]):
                    return ytuple[idx]
        return None
        
    def hasTSTT(self, yvec, idx):
        return self.getTSTT(yvec, idx) is not None
        
    def hashcode(self, yvec):
        output = 0
        for a in yvec:
            output += yvec[a] * (a.start.id * 1000 + a.end.id)
        
        return output
        
    def equals(self, yvec1, yvec2):
        for a in yvec1:
            if yvec2[a] != yvec1[a]:
                return False
        return True
    
     