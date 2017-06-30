class InclusionList:
    def read( ifile ):
        result = []
        fi = open( ifile )
        
        for line in fi:
            line = line.strip()
            
            if( line.startswith( "#" ) or (line == "") ):
                continue
            
            result.append( line )
        
        fi.close()
        
        return( result )
    
    read = staticmethod( read )

class Target:
    def __init__(self, model, cols):
        self.model = model
        self.cols = cols
    
    def check(self):
        return( not self.cols is None )

class TestSet:
    FTYPE_STK = "stk"
    FTYPE_FASTA = "fasta"
    FTYPE_INCLUSION = "in"
    FTYPE_CANDS = "cands"
    FTYPE_DISREG = "disreg"
    FTYPE_RESULTS = "results"
        
    def __init__(self, fname):
        self.fname = fname
        
        self.parse()

    def parse(self):
        fi = open( self.fname )
    
        self.ftype = fi.readline().strip()
        self.dir = fi.readline().strip()
        self.prefix = fi.readline().strip()
        self.models = map( lambda s: s.strip(), fi.readline().split( "," ) )

        self.targets = []
        for line in fi:
            data = line.strip().split()
            
            cols = None
            if( data[1] != "-" ):
                cols = map( int, data[1:] )
            
            self.targets.append( Target(data[0], cols) )
        
        fi.close()

    def get_fname(self, ftype, n=None):
        if( n is None ):
            nstr = ""
        else:
            nstr = ".%02d" %n
            
        return( "%s/%s%s.%s" %(self.dir, self.prefix, nstr, ftype) )

class Cands:
    def read( ifile ):
        fi = open( ifile )
        
        comments = ""
        result = {}
        for line in fi:
            line = line.strip()
            
            if( line == "" ):
                continue
            elif( line.startswith( "#" ) ):
                comments += line + "\n"
            else:
                dt = line.split()
                
                if( not result.has_key( dt[0] ) ):
                    result[dt[0]] = []
                
                result[dt[0]].append( Cand( dt ) )
        
        fi.close()
        
        return( comments, result )
    
    def write( dcands, ofile=None, fo=None ):
        if( not ofile is None ):
            fo = open( ofile, "w" )

        for (key, cands) in dcands.items():
            Cands.write_list( cands, fo )
        
        if( not ofile is None ):
            fo.close()
        
    
    def write_list( cands, fo ):
        for cand in cands:
            fo.write( "%s\t%s\t%s\t" %(cand.seq_key.ljust(70), cand.model_name.ljust(10), cand.model_version.ljust(4)) )
            
            if( cand.empty ):
                fo.write( "EMPTY\n" )
            else:
                fo.write( "%6.3f\t%4d\t" %(cand.score, cand.bpp) )
                fo.write( "%s\t" %"\t".join( map( lambda col: "%4d" %col, cand.cols ) ) )
                fo.write( "%s\n" %cand.seqs )

    read = staticmethod( read )
    write = staticmethod( write )
    write_list = staticmethod( write_list )

class Cand:
    def __init__(self, dt):
        self.seq_key = dt[0]
        self.model_name = dt[1]
        self.model_version = dt[2]

        self.empty = (dt[3] == "EMPTY")

        if( not self.empty ):
            self.score = float(dt[3])
            self.bpp = float(dt[4])
            self.value = float(dt[5])
            self.cols = dt[6].split( ";" )
            self.seqs = dt[7].split( ";" )

class Filter:
    def __init__(self):
        self.SCORE_LIMIT = 1.0
        self.BOLTZ_LIMIT = 1
        
    def rules(self):
        txt = "#\n"
        txt += "# Filtered all candidates that:"
        txt += "# 1. Have a score less than %d.\n" %self.SCORE_LIMIT
        txt += "# 2. Have a bpp score less than %f.\n" %self.BPP_LIMIT
        txt += "#\n"
        
        return (txt)
    
    def filter(self, dcands):
        dcands = self.filter_absolute( dcands )
        dcands = self.filter_relative( dcands )
        
        return( dcands )
    
    def filter_absolute(self, dcands):
        dcands_new = {}
        
        for (key, cands) in dcands.items():
            cands_new = []

            # filter sequence by absolute rules 
            for (i, cand) in enumerate(cands):
                # RULE 0: filter empty candidates
                if( cand.empty ):
                    continue
                
                # RULE 1: filter by score
                if( cand.score < self.SCORE_LIMIT ):
                    continue
                    
                # RULE 2: filter by boltz
                if( (cand.boltz < self.BOLTZ_LIMIT) and (not cand.auto) ):
                    continue
                
                # If reached this point we keep it
                cands_new.append( cand )
    
            dcands_new[key] = cands_new
        
        return( dcands_new )
    
    def filter_relative(self, dcands):
        dcands_new = {}
        
        # filter sequence by relative rules
        for (key, cands) in dcands.items():

            discard = []
            for i in xrange(len(cands)):
                if( not i in discard ):
                    for j in xrange(i+1, len(cands)):
                        # RULE 1: if both candidates start in the same position keep only the highest 'boltz'
                        for (ci, cj) in zip(cands[i].cols, cands[j].cols):
                            if( ci == cj ):
                                if( cands[i].boltz < cands[j].boltz ):
                                    discard.append( i )
                                else:
                                    discard.append( j )
                                    
                                break
            
            # last pass
            cands_new = []
            for (i, cand) in enumerate(cands):
                if( not i in discard ):
                    cands_new.append( cand )
            
            dcands_new[key] = cands_new
            
        return( dcands_new )

    def found( self, cand, tset ):
        ok = False
        
        for tg in tset.targets:
            if( tg.check() and (tg.model == cand.model_name) ):
                ok = True
                for (c1, c2) in zip(tg.cols, cand.cols):
                    ok = ok and (c1 == c2)
        
        return( ok )
