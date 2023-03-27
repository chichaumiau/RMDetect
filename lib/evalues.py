import math
import os

class EValues:
    def __init__(self):
        self.gc_classes = []
        self.gc_values = []
        
        self.selected = None
        
    # selects the closest gc content from the list and stores the index 
    def select_gc(self, gc):
        min_d = None
        for (i, gc_class) in enumerate(self.gc_classes):
            x = 0.0
            for (k,v) in gc_class.items():
                x += (v - math.e ** gc.get(k, 0))**2  

            d = math.sqrt( x )
            
            if( (min_d is None) or (d < min_d) ):
                min_d = d
                self.selected = i
                
    def compute(self, score):
        result = 1.0
        
        if( not self.selected is None ):
            key = "%.1f" %round(score, 1)
            
            if( key == "-0.0" ):
                key = "0.0"
            
            result = self.gc_values[self.selected].get( key, 0.0 )
        
        return( result )

class EValuesParser:
    def parse(self, fname):
        result = EValues()
        
        READ_NONE = 0
        READ_CLASSES = 1
        READ_EVALUES = 2
        
        state = READ_NONE
        
        if( os.path.isfile( fname ) ):
            fi = open(fname)
            
            for line in fi:
                line = line.strip()
                
                if( (line == "") or line.startswith( "#" ) ):
                    continue
                
                if( "GC_CLASSES" in line ):
                    state = READ_CLASSES
                elif( "E_VALUES" in line ):
                    state = READ_EVALUES
                elif( state == READ_CLASSES ):
                    data = line.split()
                    
                    gc_class = {}
                    
                    for i in range(1, len(data), 2):
                        gc_class[data[i]] = float(data[i+1])
                    
                    result.gc_classes.append( gc_class )
                    result.gc_values.append( {} )
                    
                elif( state == READ_EVALUES ):
                    data = line.split()
                    
                    for i in range(1, len(data)):
                        result.gc_values[i-1][data[0]] = float(data[i])
                        
            # TODO: check if the scores are increasing 0.1
             
        return( result )
