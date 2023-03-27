import re
import os

class MCAnnotate:
    def __init__(self):
        self.mc_file = ""
        self.interactions = []
    
    def load(self, pdb_file, mc_dir):
        # defines the annotation file
        if( mc_dir == "" ):
            mc_dir = "."
            
        self.mc_file = "%s/%s.mcout" %( (mc_dir != "" and mc_dir or ".") , os.path.basename(pdb_file))
        
        # check if the annotation file exists
        if not os.path.isfile( self.mc_file ):
            # create a new annotation file
            cmd = "MC-Annotate %s > %s" %(pdb_file, self.mc_file)
            os.system( cmd )
        
        # parse the annotation file
        self.parse()
    
    def parse(self):
        STATE_OUT = 0
        STATE_PAIR = 1
        STATE_STACK = 2
        
        pattern_pair = "^([A-Z]|\'[0-9]\'|)(\d+)-([A-Z]|\'[0-9]\'|)(\d+) : (\w+)-(\w+) ([\w\']+)/([\w\']+)(?:.*)pairing( (parallel|antiparallel) (cis|trans))"
        pattern_stack = "^([A-Z]|\'[0-9]\'|)(\d+)-([A-Z]|\'[0-9]\'|)(\d+) :.*(inward|upward|downward|outward).*"
              
        # opens and parses the annotation file
        f = open( self.mc_file, "r" )
            
        model_count = 0
        state = STATE_OUT 
        for line in f:
            line = line.strip()
    
            # for now we only care about the first model
            if( line.startswith( "Residue conformations" ) ):
                if( model_count == 0 ):
                    model_count += 1
                    continue
                else:
                    break
            
            if( line.startswith( "Base-pairs" ) ):
                state = STATE_PAIR
                continue
    
            if( line.startswith( "Adjacent stackings" ) or line.startswith( "Non-Adjacent stackings" ) ):
                state = STATE_STACK
                continue
    
            if( line.endswith( "----------" ) ):
                state = STATE_OUT
                continue
            
            interaction = None
            
            if( state == STATE_PAIR ):
                match = re.match( pattern_pair, line )
                
                if( match != None ):
                    g = match.groups()
                    interaction = self.convert_pair( g )
            
            if( state == STATE_STACK ):
                match = re.match( pattern_stack, line )
                
                if( match != None ):
                    g = match.groups()
                    interaction = self.convert_stack( g )
    
            if( interaction != None ):
                self.interactions.append( interaction )
            
        f.close()
    
    def convert_pair( self, match ):
        int_a = match[6][0].upper()
        int_b = match[7][0].upper()
        
        result = None
        
        if( (int_a in ["W", "H", "S"]) and (int_b in ["W", "H", "S"]) ):
            chain_a = match[0].replace( "'", "" )
            pos_a = int(match[1])
            nt_a = match[4]
            
            chain_b = match[2].replace( "'", "" )
            pos_b = int(match[3])
            nt_b = match[5]
            
            int_type = "%s%s" %(int_a, int_b)
            int_orientation = match[10].lower()
            
            # check if the smallest 'pos' is always the first 
            if( ((chain_a == chain_b) and (pos_a < pos_b)) or (chain_a < chain_b) ):
                result = ("PAIR", chain_a, pos_a, nt_a, chain_b, pos_b, nt_b, int_type, int_orientation, "")
            else:
                result = ("PAIR", chain_b, pos_b, nt_b, chain_a, pos_a, nt_a, int_type, int_orientation, "")
        
        return( result )
    
    def convert_stack( self, match ):
        chain_a = match[0].replace( "'", "" )
        pos_a = int(match[1])

        chain_b = match[2].replace( "'", "" )
        pos_b = int(match[3])
        
        int_type = match[4]
            
        return( "STACK", chain_a, pos_a, "", chain_b, pos_b, "", int_type, "", "" )
    

if __name__ == '__main__':
    mca = MCAnnotate()
    
    mca.load( "/media/Elements/3-PHD/resources/structures/mcannotate/1A3M.pdb", "/media/Elements/3-PHD/resources/structures/mcannotate" )
    
    for entry in mca.interactions:
        print(entry)
