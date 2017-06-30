class BPProbs:
    def __init__( self ):
        self.bases = {}
        self.matrix = {}
    
    def add(self, b1, b2, p):
        if( p != 0.0 ):
            self.matrix[b1] = self.matrix.get( b1, {} )
            self.matrix[b1][b2] = p
            
            self.matrix[b2] = self.matrix.get( b2, {} )
            self.matrix[b2][b1] = p
            
            self.bases[b1] = self.bases.get( b1, 0.0 ) + p
            self.bases[b2] = self.bases.get( b2, 0.0 ) + p

    def prob_pair(self, b1, b2 ):
        return( self.matrix.get( b1, {} ).get( b2, 0.0 ) )

    def prob_base(self, b ):
        return( self.bases.get( b, 0.0 ) )

    def prob_local(self, b1s, b2s ):
        ptotal = 0.0
        count = 0
        for b1 in b1s:
            for b2 in b2s:
                p = self.matrix.get( b1, {} ).get( b2, 0.0 )
                
                if( p != 0.0 ):
                    count += 1
                    ptotal += p
            
        return( ptotal, count )
    
    def prob_total(self, bs ):
        ptotal = 0.0
        count = 0
        for b in bs:
            p = self.bases.get( b, 0.0 )
                
            if( p != 0.0 ):
                count += 1
                ptotal += p
            
        return( ptotal, count )