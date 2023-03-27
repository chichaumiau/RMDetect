import math

class GC:
    def compute( sequence ):
        result = {}
        
        for c in sequence:
            result[c] = result.get( c, 0.0 ) + 1.0
        
        for (nt, count) in result.items():
            result[nt] = math.log( count / len(sequence) )
        
        return( result )
    
    def neutral():
        p = math.log( 0.25 )
        return( {"A":p, "C":p, "G":p, "U":p} )

    def specific(gc):
        pgc = math.log( gc/2.0 )
        pau = math.log( (1.0-gc)/2.0 )
        
        return( {"A":pau, "C":pgc, "G":pgc, "U":pau} )
    
    def get_classes(resolution=10, min_p=10, max_p=90):
        result = []
        
        r = resolution
        
        # why 4 loops? because we don't want any priviledge variable
        
        a = min_p
        while a <= max_p:
            c = min_p
            while c <= max_p:
                g = min_p
                while g <= max_p:
                    u = min_p
                    while u <= max_p:
                        if( (a+c+g+u) == 100 ):
                            result.append( {"A":float(a)/100, "C":float(c)/100, "G":float(g)/100, "U":float(u)/100 } )
                        u += r
                    g += r
                c += r
            a += r
                        
        
        return( result )
                    
            
        
    compute = staticmethod( compute )
    neutral = staticmethod( neutral )
    specific = staticmethod( specific )
    get_classes = staticmethod( get_classes )