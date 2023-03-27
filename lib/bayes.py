class Model:
    def __init__(self, data):
        self.obss = []
        self.total = 0.0
       
        self.load( data )

    # Exclude from the model all observations for which at least one of the variables 'vars' is empty 
    def exclude(self, vars):
        for obs in self.obss:
            obs.exclude = False
            for var in vars:
                if( obs.data[var] == "" ):
                    obs.exclude = True
                    break 

    # Include all observations (resets the model to the state before any call to 'exclude')  
    def include(self):
        for obs in self.obss:
            obs.exclude = False

    # Creates a new variable as the combination of previous ones
    def new_var(self, name, vars, pos_from=None, pos_to=None, keep_gaps=False ):
        for obs in self.obss:
            data = ""
            for var in vars:
                data += obs.data[var]
            if( not keep_gaps ):
                data = data.replace( ".", "" ) 

            data = data[pos_from:pos_to]
            
            if( data == "" ):
                data = "."
            
            obs.data[name] = data
        
    def load( self, data ):
        self.obss = []
        self.total = 0.0
        
        for entry in data:
            obs = Observation( entry )
            
            self.total += entry['weight']
            self.obss.append( obs )
                
    def get_vars(self):
        result = []
        if( len(self.obss) > 0 ):
            result = filter(lambda k: k != 'weight', self.obss[0].data.keys())
            
        result.sort()
        return( result )
        
class Observation:
    def __init__(self, entry):
        self.prob = 0.0
        self.exclude = False
        
        self.weight = entry['weight']

        self.data = {}

        for (k, v) in entry.items():
            self.data[k] = v
            
class JointProb:
    def __init__(self, model, vars ):
        self.model = model
        self.vars = vars
        
        self.probs = {}
        
        self.__compute()
        
    def show(self, limit=1E-3, sort_by_value=True):
        print("P( %s ):" %(", ".join( self.vars )))
        
        residues = 0.0
        
        # sort data
        data = []
        
        if( sort_by_value ):
            for k, v in self.probs.items():
                data.append( (k, v) )
                
            data.sort( cmp=lambda x, y: int( x[1] * 100000.0 - y[1] * 100000.0 ), reverse=True )
        else:
            keys = self.probs.keys()
            keys.sort()
            
            for key in keys:
                data.append( (key, self.probs[key]) )  
            
        # show data
        for (k, v) in data:
            if( v > limit ):
                print(",".join( k ), "== %.4f" %v)
            else:
                residues += v
        
        print("residues (< %f): %.6f\n" %(limit, residues))

    def __compute(self):
        self.probs = {}
        
        total = 0.0
        for obs in self.model.obss:
            data = []
            for var in self.vars:
                dt = obs.data[var]
                
                if( (dt == ".") and obs.exclude ):
                    break
                else:
                    data.append( dt )

            tdata = tuple( data )

            self.probs[tdata] = self.probs.get( tdata, 0.0 ) + obs.weight
            total += obs.weight
        
        # discount all parameters that are too small  
        total = self.__discount( total )

        # compute final parameters
        for (k, v) in self.probs.items():
            self.probs[k] = (v / total)

    def __is_gap(self, k):
        return ("".join(k).replace( ".", "" ) == "")
    
    def __discount(self, total):
        new_total = total
        
        for (k, v) in self.probs.items():
            p = (v / total)
            
            gap = self.__is_gap( k )
            if( (not gap and (p < 0.001)) or (gap and (p < 0.010))):
                new_total -= v
                del self.probs[k]
                
        return( new_total )

class MultProb:
    def __init__(self, model, vars ):
        self.model = model
        self.vars = vars
        
        self.probs = {}
        
        self.__compute()
        
    def show(self, limit=1E-3, sort_by_value=True):
        print("P( %s ):" %(", ".join( self.vars )))
        
        residues = 0.0
        
        # sort data
        data = []
        
        if( sort_by_value ):
            for k, v in self.probs.items():
                data.append( (k, v) )
                
            data.sort( cmp=lambda x, y: int( x[1] * 100000.0 - y[1] * 100000.0 ), reverse=True )
        else:
            keys = self.probs.keys()
            keys.sort()
            
            for key in keys:
                data.append( (key, self.probs[key]) )  
            
        # show data
        for (k, v) in data:
            if( v > limit ):
                print(k, "== %.4f" %v)
            else:
                residues += v
        
        print("residues (< %f): %.6f\n" %(limit, residues))

    def __compute(self):
        self.probs = {}

        # get the individual probabilities
        jps = []
        jkeys = []
        pointers = []
        
        for var in self.vars:
            jp = JointProb( self.model, [var] )
            
            jps.append( jp )
            keys = jp.probs.keys()
            jkeys.append( keys )
            
            if( len(keys) > 0 ):
                pointers.append( 0 )
            else:
                pointers.append( None )

        ok = True
        while ok:
            key = []
            prob = 1.0
            for (i, jkey) in enumerate(jkeys):
                if( not pointers[i] is None ):
                    k = jkey[pointers[i]]
                    p = jps[i].probs[k]
                    
                    key.append( k )
                    prob *= p
            
            self.probs[tuple(key)] = prob
            
            # increment pointers
            ok = False
            for i in xrange(len(pointers)):
                if( not pointers[i] is None ):
                    if( pointers[i] + 1 < len(jkeys[i]) ):
                        pointers[i] += 1
                        ok = True
                        break
                    else:
                        pointers[i] = 0
    
class CondProb:
    def __init__(self, model, vars_a, vars_b ):
        self.model = model
        self.vars_a = vars_a
        self.vars_b = vars_b
        
        self.probs = {}
        
        self.__compute()
        
    def show(self, limit=1E-3, sort_by_value=True):
        print("P( %s | %s ):" %(", ".join( self.vars_a ), ", ".join( self.vars_b )))
        
        residues = 0.0
        
        # sort data
        data = []
        
        if( sort_by_value ):
            for k, v in self.probs.items():
                data.append( (k, v) )
                
            data.sort( cmp=lambda x, y: int( x[1] * 100000.0 - y[1] * 100000.0 ), reverse=True )
        else:
            keys = self.probs.keys()
            keys.sort()
            
            for key in keys:
                data.append( (key, self.probs[key]) )  
            
        # show data
        for (k, v) in data:
            if( v > limit ):
                print(",".join( k[0] ), "|", ",".join( k[1] ), "== %.4f" %v)
            else:
                residues += v
        
        print("residues (< %f): %.6f\n" %(limit, residues))
        
    def __compute(self):
        jp_ab = JointProb( self.model, self.vars_a + self.vars_b )
        jp_b = JointProb( self.model, self.vars_b )
        
        self.probs = {}
        
        for (k, v) in jp_ab.probs.items():
            k1 = k[:len(self.vars_a)]
            k2 = k[len(self.vars_a):]
            
            kc = (k1, k2)
            
            # conditional probability formula
            if( jp_b.probs.has_key( k2 ) ):
                self.probs[kc] = (v / jp_b.probs[k2])

            #print(k, "/",  k1, "/", k2, "/", kc, "/", v, "/", jp_b.probs[k2], "/", self.probs[kc])
