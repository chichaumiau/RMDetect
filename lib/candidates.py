import math
import sys

from pair import *
from rnatools import *
#
# auxiliar merge pair function
#
def merge_pairs( pairs ):
    pairs_dict = pairs_2_dict( pairs )
    
    # check conflicting base pairs
    keys = pairs_dict.keys()
    keys.sort()
    
    for (o1, c1) in keys:
        v1 = pairs_dict.get( (o1, c1), None )
        if( v1 is None ):
            continue
        
        for (o2, c2) in keys:
            v2 = pairs_dict.get( (o2, c2), None )
            if( v2 is None ):
                continue

            # self reference
            if( (o1 == o2) and (c1 == c2) ):
                continue
            
            # overlapping 
            ok1 = not( (o1 == o2) or (o1 == c2) or (c1 == o2) or (c1 == c2) )

            # pseudo knots
            ok2 = (c2 < o1) or (o2 > c1) or ((o2 > o1) and (c2 < c1)) or ((o2 < o1) and (c2 > c1))
            
            if( not ok1 or not ok2 ):
                if (v1 < v2):
                    k = (o1, c1)
                elif (v1 > v2):
                    k = (o2, c2)
                elif(o1 > o2):
                    k = (o1, c1)
                else:
                    k = (o2, c2)
                
                if( not pairs_dict.get( k, None ) is None ):
                    del pairs_dict[k]

    return( pairs_dict.keys() )

def pattern_2_pairs( pat ):
    pairs = []
    stack = []
    for (i, c) in enumerate(pat):
        if( c == "(" ):
            stack.append( i )
        if( c == ")" ):
            if( len(stack) > 0 ):
                j = stack.pop()
                pairs.append( (j, i) )
            else:
                return False
    
    if( len(stack) == 0 ):
        return( pairs )
    else:
        return( false )

def pairs_2_dict( pairs ):
    pairs_dict = {}

    for k in pairs:
        pairs_dict[k] = pairs_dict.get( k, 0 ) + 1
    
    return( pairs_dict )

class Candidates:
    def filter(cands, min_score=None, min_bpp=None, min_evalue=None, cback_after_cand=None):
        result = []

        for (i, cand) in enumerate(cands):
            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )

            if( (not min_score is None) and (cand.score <= min_score) ):
                continue
        
            if( (not min_bpp is None) and (cand.bpp <= min_bpp) ):
                continue
        
            if( (not min_evalue is None) and ((cand.evalue == 0.0) or (-math.log(cand.evalue) <= min_evalue)) ):
                continue
            
            result.append( cand )
                    
        result.sort()
        return( result )

    def top_model(cands, N, cback_after_cand=None):
        cands.sort()
        
        result = []
        dict = {}

        for (i, cand) in enumerate(cands):
            key = cand.model.full_name()
            
            dict[key] = dict.get( key, 0 ) + 1
            
            if( dict[key] <= N ):
                result.append( cand ) 

            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )    

        return( result )

    def top_seq(cands, N, cback_after_cand=None):
        cands.sort()
        
        result = []
        dict = {}

        for (i, cand) in enumerate(cands):
            key = cand.key
            
            dict[key] = dict.get( key, 0 ) + 1
            
            if( dict[key] <= N ):
                result.append( cand ) 
                
            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )
                    
        return( result )
    
    def top_seq_model(cands, N, cback_after_cand=None):
        cands.sort()
        
        result = []
        dict = {}

        for (i, cand) in enumerate(cands):
            key = cand.key + "_" + cand.model.full_name()
            
            dict[key] = dict.get( key, 0 ) + 1
            
            if( dict[key] <= N ):
                result.append( cand ) 

            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )
    
        return( result )

    def get_model(cands, model_name, cback_after_cand=None):
        cands.sort()
        
        result = []
        
        for (i, cand) in enumerate(cands):
            if( cand.model.full_name().upper() == model_name.upper() ):
                result.append( cand )

            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )
 
    
        return( result )

    def nooverlap(cands, cback_after_cand=None):
        for i in range(len(cands)):
            if( cands[i] is None ):
                continue
            
            for j in range(i+1, len(cands)):
                if( cands[j] is None ):
                    continue

                dup = True
                (s1, e1) = cands[i].coords
                (s2, e2) = cands[j].coords
                
                if( (s1 < e2) and (s2 < e1) and (cands[i].key == cands[j].key) ):
                    # to be overlap all hte strands must overlap
                    for (pi, pj) in zip(cands[i].chains_pos, cands[j].chains_pos):
                        # pi[0] ->  first position of the first strand
                        # pi[-1] -> last position of the first strand
                        # pj[0] ->  first position of the second strand
                        # pj[-1] -> last position of the second strand
                        if( (pi[0] > pj[-1]) or (pi[-1] < pj[0]) ):
                            dup = False
                            break
                else:
                    dup = False
                
                if( dup ):
                    if( (cands[i].score > cands[j].score) and (cands[i].bpp > 0) ):
                        cands[j] = None
                    else:
                        cands[i] = None
                        break
                
            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )
        
        return( filter( lambda x: not x is None, cands ) )

    def compute_bpp(cands, seqs, constraint="", cback_after_cand=None):
        for (i, cand) in enumerate(cands):
            if( cand.bpp is None ):
                sequence = seqs.get_sequence( cand.key )
                
                bpp = cand.compute_bpp( sequence.seq, constraint )
                
                if( bpp == 0.0 ):
                    cands[i] = None
                else:
                    cand.bpp = bpp
                    
            if( not cback_after_cand is None ):
                cback_after_cand( i, len(cands) )
                        
        return( filter( lambda x: not x is None, cands ) )

    # get the consensus pattern for all candidates of the list
    # the parameter alignment will force the function to work with
    # alignment coordinates instead of sequence coordinates 
    def get_pattern_2d(cands, seqs, constraint="", align=False, strict=True):
        pairs = []
        singles = []
        max_len = 0
        
        # collects and counts all base pairs
        for cand in cands:
            sequence = seqs.get_sequence( cand.key )
            
            if( align ):
                cand.set_pairing_align()
                max_len = max( max_len, len(sequence.seq_gapped) )
                pairs.extend( cand.pairs_align )
                singles.extend( cand.singles_align )
            else:
                cand.set_pairing( strict=strict )
                max_len = max( max_len, len(sequence.seq) )
                pairs.extend( cand.pairs )
                singles.extend( cand.singles )
        
        # add the constraints pairs if there are constraints
        if( constraint != "" ):
            pairs.extend( pattern_2_pairs( constraint ) )

        pairs = merge_pairs( pairs )
        
        # get the brackets notation pattern
        bn_pattern = ["."] * (max_len)
        
        for (o, c) in pairs:
            bn_pattern[o] = "("
            bn_pattern[c] = ")"
        
        for i in singles:
            if( bn_pattern[i] == "." ):
                bn_pattern[i] = "x"
        
        if( constraint != "" ):
            for (i, c) in enumerate(constraint):
                if( (bn_pattern[i] == ".") and (not c in "()") ):
                    bn_pattern[i] = c

        return( "".join( bn_pattern ) )

    def get_postscript_2d(pattern_2d):
        ps_cmd = ""
        colors = ["1 0 0", "0 1 0", "0 0 1", "1 1 0", "1 0 1", "0 1 1", "0.5 0 0", "0 0.5 0", "0 0 0.5", "0.5 0.5 0", "0.5 0 0.5", "0 0.5 0.5"]
        
        last_i = 0
        label_i = 0
        pre_stack = None
        stack = []
        
        (stem_label_i, stem_label_j) = (-1, -1)
        stem = None
        count = 0
        
        for (i, c) in enumerate(pattern_2d):
            if( c in "()x" ):
                if( i > (last_i + 1) ):
                    label_i += 1
                last_i = i
            
            if( c == "(" ):
                stack.append( (i, label_i) )
                
            if( c == ")" ):
                (j, label_j) = stack.pop()
                
                # saves the previous stem?
                if( (not stem is None) and ((stem_label_i != label_i) or (stem_label_j != label_j)) ):
                    ncolor = colors[count % len(colors)]
                    
                    ps_cmd += "%d %d 10 %s omark\n" %(stem[0]+1, stem[1]+1, ncolor)
                    ps_cmd += "%d %d 10 %s omark\n" %(stem[2]+1, stem[3]+1, ncolor)
                    stem = None
                    count += 1

                # creates a new stem
                if( stem is None ):
                    (stem_label_i, stem_label_j) = (label_i, label_j)
                    stem = [j, j, i, i]
                    
                # updates the current stem
                else:
                    stem[0] = j
                    stem[3] = i

        # add the last stem
        if( not stem is None ):
            ncolor = colors[count % len(colors)]
            
            ps_cmd += "%d %d 10 %s omark\n" %(stem[0]+1, stem[1]+1, ncolor)
            ps_cmd += "%d %d 10 %s omark\n" %(stem[2]+1, stem[3]+1, ncolor)
        
        return( ps_cmd )
        
    #
    #
    filter = staticmethod( filter )
    top_seq = staticmethod( top_seq )
    top_model = staticmethod( top_model )
    top_seq_model = staticmethod( top_seq_model )
    get_model = staticmethod( get_model )
    nooverlap = staticmethod( nooverlap )
    compute_bpp = staticmethod( compute_bpp )
    get_pattern_2d = staticmethod( get_pattern_2d )
    get_postscript_2d = staticmethod( get_postscript_2d )


class Candidate:
    def __init__(self, model, key, coords):
        self.model = model
        self.key = key
        self.coords = coords
        
        self.chains_seqs = []
        self.chains_pos = []
        self.chains_cols = []
        
        self.score = None
        self.evalue = None
        self.bpp = None
        
        self.pairs = None
        self.pairs_align = None
        self.singles = None

    #
    # sets the candidate data 
    #
    def set_data_build(self, prob_model, prob_rand, chains_seqs, chains_pos, evalues):
        self.chains_seqs = chains_seqs
        self.chains_pos = chains_pos
        self.chains_cols = []
        
        self.score = (prob_model - prob_rand) / math.log(2.0)
        self.evalue = evalues.compute( self.score )

    def set_data_load(self, score, evalue, bpp, chains_pos, chains_seqs, chains_cols ):
        self.score = score
        self.evalue = evalue
        self.bpp = bpp

        self.chains_pos = chains_pos
        self.chains_seqs = chains_seqs
        self.chains_cols = chains_cols

    def __str__(self):
        return( CandidateParser.code( self ) )
    
    def __cmp__(self, other):
        if( self.key > other.key ):
            return( 1 )
        elif( self.key < other.key ):
            return( -1 )
        elif( self.evalue > other.evalue ):
            return( 1 )
        elif( self.evalue < other.evalue ):
            return( -1 )
        elif( self.score > other.score ):
            return( -1 )
        elif( self.score < other.score ):
            return( 1 )
        else:
            return( 0 )

    def compute_bpp(self, seq, constraint=""):
        self.set_pairing()
        
        # gets the candidate sequence to fold
        seq = seq[self.coords[0]:self.coords[1]]
        
        # folds with constraints and get the ensemble MFE
        # add additional imposed constraints
        pat_2d = self.get_pattern_2d( constraint )
        
        fold_c = RNATools().fold(seq, constraint=pat_2d, ensemble=True)
        
        # folds without constraints (or with the additional cnstrains only) and get the ensemble MFE
        fold_a = RNATools().fold(seq, constraint=constraint, ensemble=True)
        
        # computes the new joint bpp (thanks to Rolf Backofen)
        # kT = (T + 273.5 x 1.98717) / 1000; T = 37 
        kT = 0.61702 
        bpp = math.e ** ((fold_a.mfe - fold_c.mfe) / kT)
        
        if( bpp > 1.0 ):
            bpp = 0.0
        
        return( bpp )
    
    def get_pattern_2d(self, constraint=""):
        self.set_pairing()

        # gets the pairing constraints
        constraint = constraint[self.coords[0]:self.coords[1]]
        
        if( constraint != "" ):
            self.pairs.extend( pattern_2_pairs( constraint ) )
            self.pairs = merge_pairs( self.pairs )
        
        # get the brackets notation pattern
        bn_pattern = ["."] * (self.coords[1]-self.coords[0])
        
        for (p1, p2) in self.pairs:
            bn_pattern[p1] = "("
            bn_pattern[p2] = ")"

        for p in self.singles:
            bn_pattern[p] = "x"

        if( constraint != "" ):
            for (i, c) in enumerate(constraint):
                if( (bn_pattern[i] == ".") and (not c in "()") ):
                    bn_pattern[i] = c

        return( "".join( bn_pattern ) )

    def get_pattern_2d_align(self):
        self.set_pairing_align()
        
        # get the brackets notation pattern
        bn_pattern = []
        
        for (p1, p2) in self.pairs_align:
            diff = (p2 - len(bn_pattern)) + 1
            
            if( diff > 0 ):
                bn_pattern.extend( ["."] * diff )
                
            bn_pattern[p1] = "("
            bn_pattern[p2] = ")"

        for p in self.singles_align:
            diff = (p - len(bn_pattern)) + 1
            
            if( diff > 0 ):
                bn_pattern.extend( ["."] * diff )

            bn_pattern[p] = "x"

        return( "".join( bn_pattern ) )

    def get_postscript_2d(self, ncolor):
        # get the post script command
        ps_cmd = ""
        colors = ["1 0 0", "0 1 0", "0 0 1", "1 1 0", "1 0 1", "0 1 1", "0.5 0 0", "0 0.5 0", "0 0 0.5", "0.5 0.5 0", "0.5 0 0.5", "0 0.5 0.5"]
        
        ncolor = ncolor % len(colors)
        
        for (i, chain_pos) in enumerate(self.chains_pos):
            start = chain_pos[0] - self.coords[0] + 1
            end = chain_pos[-1] - self.coords[0] + 1
            
            ps_cmd += "%d %d 10 %s omark\n" %(start, end, colors[ncolor])  
        
        return( ps_cmd )

    
    def update_chains_cols(self, cols):
        self.chains_cols = []
        
        for chain_pos in self.chains_pos:
            chain_cols = []
            
            for pos in chain_pos:
                
                if( pos is None ):
                    col = None
                elif( cols is None ):
                    col = pos
                else:
                    col = cols[pos]
                    
                chain_cols.append( col )
            
            self.chains_cols.append( chain_cols )
            
    def update_chains_pos(self, offset):
        self.coords = (self.coords[0] + offset , self.coords[1] + offset)
        
        for chain_pos in self.chains_pos:
            for (i, pos) in enumerate(chain_pos):
                if( not pos is None ):
                    chain_pos[i] = pos + offset
      
    def get_pos(self, c, ndx):
        return( self.chains_pos[c][ndx] )

    def get_col(self, c, ndx):
        return( self.chains_cols[c][ndx] )
    
    def get_nt(self, c, ndx):
        return( self.chains_seqs[c][ndx] )
    
    def get_seq(self, c):
        return( "".join( self.seqs[c] ) )

    def get_first_pos(self, c):
        return( self.get_N_pos( c, 1 ) )

    def get_last_pos(self, c):
        return( self.get_N_pos( c, -1 ) )

    def get_N_pos(self, c, step):
        pos = None
        for pos in self.chains_pos[c][::step]:
            if( pos != None ):
                break

        return( pos )
    
    # strict mode:
    #    base pairs of type MUST -> ( )
    #    base pairs of type CAN  -> x x
    #    singles                 -> x
    # not strict mode:
    #    base pairs of type MUST -> ( )
    #    base pairs of type CAN  -> ( )
    #    singles                 -> x
    def set_pairing(self, strict=True):
        self.pairs = []
        self.singles = []
        
        for pair in self.model.pairing:
            n1 = self.model.nodes[pair.id1]
            n2 = self.model.nodes[pair.id2]
            
            p1 = self.chains_pos[n1.chain][n1.ndx]
            p2 = self.chains_pos[n2.chain][n2.ndx]
            
            if( p1 is None or p2 is None ):
                continue
            
            p1 = p1 - self.coords[0]
            p2 = p2 - self.coords[0]
            
            nt1 = (self.chains_seqs[n1.chain][n1.ndx])
            nt2 = (self.chains_seqs[n2.chain][n2.ndx])
            
            if( strict ):
                if( pair.ptype == Pair.PTYPE_MUST ):
                    self.pairs.append( (min(p1, p2), max(p1, p2) ) )
                elif( pair.ptype == Pair.PTYPE_CAN ):
                    if( (nt1 + nt2) in ["AU", "CG", "GC", "GU", "UA", "UG"] ):
                        self.pairs.append( (min(p1, p2), max(p1, p2) ) )
                    else:
                        self.singles.append( min(p1, p2) )
                        self.singles.append( max(p1, p2) )
            else:
                self.pairs.append( (min(p1, p2), max(p1, p2) ) )

        for id in self.model.unpaired:
            n = self.model.nodes[id]
            p = self.chains_pos[n.chain][n.ndx]
            
            if( not p is None ):
                self.singles.append( p - self.coords[0] )

    def set_pairing_align(self):
        self.pairs_align = []
        self.singles_align = []
        
        for pair in self.model.pairing:
            n1 = self.model.nodes[pair.id1]
            n2 = self.model.nodes[pair.id2]
            
            p1 = self.chains_cols[n1.chain][n1.ndx]
            p2 = self.chains_cols[n2.chain][n2.ndx]
            
            if( p1 is None or p2 is None ):
                continue
            
            nt1 = (self.chains_seqs[n1.chain][n1.ndx])
            nt2 = (self.chains_seqs[n2.chain][n2.ndx])
            
            if( pair.ptype == Pair.PTYPE_MUST ):
                self.pairs_align.append( (min(p1, p2), max(p1, p2) ) )
            elif( pair.ptype == Pair.PTYPE_CAN ):
                if( (nt1 + nt2) in ["AU", "CG", "GC", "GU", "UA", "UG"] ):
                    self.pairs_align.append( (min(p1, p2), max(p1, p2) ) )
                else:
                    self.singles_align.append( min(p1, p2) )
                    self.singles_align.append( max(p1, p2) )

        for id in self.model.unpaired:
            n = self.model.nodes[id]
            p = self.chains_cols[n.chain][n.ndx]
            
            if( not p is None ):
                self.singles_align.append( p )

            
class CandidateParser:
    def parse( str, models ):
        # TODO: check (with a reg. exp.) the syntax of each the line
        data = str.split()
        
        # search for the model in the model list
        model = None
        for m in models:
            if( (m.name == data[0]) and (m.version == data[1]) ):
                model = m
                break
        
        if( model is None ):
            sys.stderr.write( "ERROR >> CandidateParser.parse >> Model '%s' '%s' not found!\n" %(data[0], data[1]) )
            quit()
        
        coords = map( int, data[3].split( "-" ) )
        
        chains_seqs = []
        for dt in data[8].split(";"):
            chains_seqs.append( list(dt) )
            
        chains_pos = CandidateParser.parse_chains( data[9] )
        chains_cols = CandidateParser.parse_chains( data[10] )

        cand = Candidate( model, data[2], coords )
        cand.set_data_load(float(data[4]), float(data[5]), float(data[6]), chains_pos, chains_seqs, chains_cols )
        
        return( cand )
    
    def code( cand ):
        aux_func = lambda l: filter( lambda x:x is not None, l )[0]

        s = "%s\t%s\t" %(cand.model.name, cand.model.version)
        s += "%s\t%d-%d\t" %(cand.key, cand.coords[0], cand.coords[1])
        s += "%.5E\t%.5E\t%.5E\t" %(cand.score, cand.evalue, cand.bpp)
        s += "%s\t" %( ";".join( ["%d/%d" %(aux_func(p), aux_func(c)) for (p, c) in zip(cand.chains_pos, cand.chains_cols)] ) )
        s += "%s\t" %( ";".join( ["".join( s ) for s in cand.chains_seqs] ) )
        s += "%s\t" %CandidateParser.code_chains( cand.chains_pos )
        s += "%s\t" %CandidateParser.code_chains( cand.chains_cols )

        #s += "%s\t" %( ";".join( [",".join( map( str, p ) ) for p in self.chains_pos] ) )
        #s += "%s\t" %( ";".join( [",".join( map( str, c ) ) for c in self.chains_cols] ) )
        
        return( s )

    def parse_chains( str ):
        data = str.split( "/" )
        data = map( lambda dt: dt.split( ";" ), data )
        
        for dt in data:
            for (i, n) in enumerate(dt):
                if( n == "-" ):
                    dt[i] = None
                elif( n == "" ):
                    dt[i] = dt[i-1] + 1
                else:
                    dt[i] = int(dt[i]) 

        return( data )
    
    def code_chains( chains ):
        s1 = []
        for chain in chains:
            s2 = []
            for (i, n) in enumerate(chain):
                if n is None:
                    s2.append( "-" )
                elif( (i>0) and (n-1 == chain[i-1]) ):
                    s2.append( "" )
                else:
                    s2.append( str(chain[i]) )
                    
            s1.append( ";".join( s2 ) )
            
        return( "/".join( s1 ) )
        
    
    parse = staticmethod( parse )
    code = staticmethod( code )
    code_chains = staticmethod( code_chains )
    parse_chains = staticmethod( parse_chains )
