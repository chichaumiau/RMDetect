import time

class TimeEstimate:
    def __init__(self, total):
        self.total = total
        
        self.t0 = 0.0
        self.telapsed = 0.0
        self.tmean = 0.0
        self.texpected = 0.0
    
    def start(self):
        self.t0 = time.time()
    
    def update(self, current):
        if( current > 0 ):
            self.telapsed = (time.time() - self.t0)
    
            self.tmean = self.telapsed / float(current)
            self.texpected = int(self.tmean * float(self.total - current))

    def elapsed(self):
        return( self.telapsed )
        
    def expected(self):
        return( self.texpected )