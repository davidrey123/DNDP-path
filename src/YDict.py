# Created on : May 22, 2024, 9:54:48 AM
# Author     : michaellevin

class YDict:
    def __init__(self):
        self.tstt_map = {}
        
    def insertSO(self, y, tstt):
        self.insertTSTT(y, tstt, 1)
            
    def insertUE(self, y, tstt):
        self.insertTSTT(y, tstt, 2)
    
    def hasSO(self, y):
        return self.hasTSTT(y, 1)
        
    def hasUE(self, y):
        return self.hasTSTT(y, 2)
        
    def getSO(self, y):
        return self.getTSTT(y, 1)
    
    def getUE(self, y):
        return self.getTSTT(y, 2)
        
    
    def insertTSTT(self, y, tstt, idx):
        hash = self.hashcode(y)
        
        ylist = None
        
        if hash in self.tstt_map:
            ylist = self.tstt_map[hash]
        else:
            ylist = []
            self.tstt_map[hash] = ylist
            
        for ytuple in ylist:
            if self.equals(y, ytuple[0]):
                ytuple[idx] = tstt
                return
        
        newtuple = [y, None, None]
        newtuple[idx] = tstt
        ylist.append(newtuple)
                
    def getTSTT(self, y, idx):
        hash = self.hashcode(y)
        
        if hash in self.tstt_map:
            ylist = self.tstt_map[hash]

            for ytuple in ylist:
                if self.equals(y, ytuple[0]):
                    return ytuple[idx]
        return None
        
    def hasTSTT(self, y, idx):
        return self.getTSTT(y, idx) is not None
        
    def hashcode(self, y):
        output = 0
        for a in y:
            output += y[a] * (a.start.id * 1000 + a.end.id)
        
        return output
        
    def equals(self, y1, y2):
        for a in y1:
            if y2[a] != y1[a]:
                return False
        return True
    
     