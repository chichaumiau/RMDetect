class Fold:
    def __init__( self, struct, mfe=0.0 ):
        self.struct = struct
        self.mfe = mfe
        
        self.tpairs = []    # list of all tuples (open base, close base) 
        self.bpairs = []    # list of all bases forming a pair
        self.dpairs = {}    # dict indexed base base number stores the pair number
        self.lpairs = {}    # dict indexed base base number stores the pairing base 

        self.parse()
        
    # computes how many bases of the current structure are correctly paired/unpaired in the other one
    # this comparison only makes sense if both structures correspond to the same underlying sequence
    def compare(self, other ):
        ok = 0.0

        for i in xrange( len(self.struct) ):
            if( self.lpairs[i] == other.lpairs[i] ):
                ok += 1.0
        
        return( ok / float(len(self.struct)) )
    
    def is_pair(self, i, j ):
        return( (min(i, j), max(i, j)) in self.tpairs )

    def is_single(self, i ):
        return( not (i in self.bpairs) )

    def parse( self ):
        self.bpairs = []
        self.dpairs = {}
        self.lpairs = {}

        count = 0
        stack = []
        
        self.lpairs = [None] * len(self.struct) 

        for (i, char) in enumerate( self.struct ):
            if( char == "(" ):
                self.bpairs.append( i )
                stack.append( i )
                
            elif( char == ")" ):
                self.bpairs.append( i )
                if( len(stack) == 0 ):
                    print "Fold.parse() > ERROR: Unexpected ')' in position '%d'" %i
                    quit()
                    
                j = stack.pop()
                
                self.tpairs.append( (min(i, j), max(i, j)) )

                count += 1
                self.dpairs[i] = count
                self.dpairs[j] = count
                
                self.lpairs[i] = j
                self.lpairs[j] = i
            elif( char == "." ):
                pass
            
            else:
                print "Fold.parse() > ERROR: Bad character '%s' in position '%d'" %(char, i)
                quit()
    
        if( len(stack) > 0 ):
            print "ERROR: Unclosed parenthesis in positions '%s'" %(", ".join( stack ))
            quit()
