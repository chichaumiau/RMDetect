import csv
import os
import sys

class Fasta:
    def __init__(self, fi):
        self.seqs = []
        self.index = {}
        
        self.parse( fi )
    
    def parse(self, fi):
        for line in fi:
            line = line.strip()
            
            if( line == "" ):
                continue
            elif( line.startswith( ">" ) ):
                key = line[1:].split()[0]
                self.index[key] = len(self.seqs)

                self.seqs.append( [key, ""] )
            else:
                self.seqs[-1][1] += line

    def write(self, fo):
        for (key, seq) in self.seqs:
            fo.write( ">%s\n%s\n" %(key, seq) )
        
    def add_seq(self, key, seq):
        i = self.index.get( key, None )
        
        if( i is None ):                
            self.index[key] = len(self.seqs)
            self.seqs.append( (key, seq) )
        else:
            self.seqs[i][1] = seq
            
class Stk:
    NUCLEOTIDES = "ACUG"
    GAP_STANDARD = "."
    GAP_NON_STANDARD = "-_:,~"
    GAP_ALLOWED = GAP_STANDARD + GAP_NON_STANDARD
    CHAR_ALLOWED = NUCLEOTIDES + GAP_ALLOWED
    
    OPEN_STRAND = "(<{[ABCDEFGHIJ"
    CLOSE_STRAND = ")>}]abcdefghij"
    VALID_PAIR = ["()", "<>", "{}", "[]", "Aa", "Bb", "Cc", "Dd", "Ee", "Ff", "Gg", "Hh", "Ii", "Jj"]
    #MY_DATA_KEYS = ["CONS", "COEV", "BPSF", "CROSS", "CROSS_CONS", "CROSS_COEV", "CROSS_BPSF", "INT2_PROT", "INT2_RNA"]
    
    KEY_REPLACE_CHARS = [(">", "_"), ("<", "_"), ("+", ""), (" ", "_"), ("#", ""), ("\"", ""), ("'", ""), ("[", "("), ("]", ")"), ("&", "_")]
    
    MAX_LOOP = 5

    def __init__(self, fi=None, char_alowed_extra=""):
        self.CHAR_ALLOWED_EXTRA = char_alowed_extra
        
        self.gc_mydata = {}
        if( (fi is None) or (not self.parse_file( fi )) ):
            self.max_key_len = 15
            self.secondary_structure = ""
            self.rf_structure = ""
            self.seqs = []
            self.index = {}
            self.dstrands = []
            self.sstrands = []
            self.pairs = []
            #for key in Stk.MY_DATA_KEYS:
            #    self.gc_mydata[key] = ""
        else:
            self.parse_ss()
        
    def add_seq(self, key, sequence):
        if( not key in self.index.keys() ):
            # if it's a new key

            # update the index
            self.index[key] = len(self.seqs)

            # insert the new sequence
            self.seqs.append( (key, sequence) )

            # update the max key length
            self.max_key_len = max( len(key), self.max_key_len )
        else:
            # if it's an existing key
            i = self.index[key]
            
            # update the sequence
            (k, s) = self.seqs[i]
            self.seqs[i] = (key, s + sequence)

    def add_secondary_structure(self, secondary_structure):
        self.secondary_structure = secondary_structure

    def add_rf_structure(self, rf_structure):
        self.rf_structure = rf_structure
    
    def add_gc_mydata(self, key=None, data=None, dict=None):
        if( dict is None ):
            self.gc_mydata[key] = data
        else:
            self.gc_mydata = dict

    def get_gc_mydata(self, key=None):
        if( key is None ):
            return( self.gc_mydata )
        else:
            return( self.gc_mydata[key] )
    
    def get_cross_list(self):
        result = {}
        
        if( self.gc_mydata.has_key( "CROSS" ) ):
            cross = self.gc_mydata["CROSS"]
            
            for (i, c) in E_(cross):
                if( not c in Stk.GAP_ALLOWED ):
                    if( not result.has_key( c ) ):
                        result[c] = []
                    result[c].append( i )
        
        return( result )
        
    def del_seq(self, key):
        # if the key exists
        if( self.index.has_key( key ) ):
            i = self.index[key]
            
            # delete from the sequence list
            del self.seqs[i]
            
            # update the index
            self.index = {}
            self.max_key_len = 0
            for i in range( 0, len(self.seqs) ):
                (k, s) = self.seqs[i]
                self.index[k] = i
                self.max_key_len = max( len(k), self.max_key_len )
        
    def normalize_secondary_structure(self, gap=GAP_STANDARD, open_strand="<", close_strand=">"):
        new_secondary_structure = ""
        
        for s in self.secondary_structure:
            if( s in self.GAP_ALLOWED ):
                new_secondary_structure += gap
            elif( s in self.OPEN_STRAND ):
                new_secondary_structure += open_strand
            elif( s in self.CLOSE_STRAND ):
                new_secondary_structure += close_strand
            else:
                sys.stderr.write( "Bad char in secondary structure '%s'\n" %s )
                return False
        
        self.secondary_structure = new_secondary_structure

    def normalize_sequences(self, gap=GAP_STANDARD, upper=False):
        for (i, (key, seq)) in enumerate(self.seqs):
            self.seqs[i] = (key, self.normalize_sequence( seq, gap, upper ))

    def normalize_sequence(self, seq, gap=GAP_STANDARD, upper=False):
        # Capitalizes sequence
        if( upper ):
            seq = seq.upper()
        
        # changes all gaps for the standard one
        for old_gap in self.GAP_ALLOWED:
            seq = seq.replace( old_gap, gap )
        
        # changes all 'T' by 'U'
        seq = seq.replace( "T", "U" )
        seq = seq.replace( "t", "u" )
        
        return( seq )
    
    def clean_sequence(self, seq):
        # changes all gaps for the standard one
        for old_gap in self.GAP_ALLOWED:
            seq = seq.replace( old_gap, "" )
        
        # changes all 'T' by 'U'
        seq = seq.replace( "T", "U" )
        
        return( seq )
        
    # TODO: add headers to the file
    def write_file(self, fname):
        fo = open( fname, "w" )
        fo.write( "# STOCKHOLM 1.0 MKL=%d\n" %self.max_key_len )
        
        for seq in self.seqs:
            fo.write( "%s %s\n" %(seq[0].ljust( self.max_key_len ), seq[1]) )
        
        if( self.secondary_structure != "" ):
            fo.write( "%s %s\n" %("#=GC SS_cons".ljust( self.max_key_len ), self.secondary_structure) )

        if( self.rf_structure != "" ):
            fo.write( "%s %s\n" %("#=GC RF".ljust( self.max_key_len ), self.rf_structure) )
            
        for (key, value) in self.gc_mydata.items():
            if( value != "" ):
                name = "#=GC %s" %key 
                fo.write( "%s %s\n" %(name.ljust( self.max_key_len ), value) )
            
        fo.write( "//\n" )
        
        fo.close()
        
    def parse_file(self, fi):
        self.secondary_structure = ""
        self.rf_structure = ""
        self.seqs = []
        self.index = {}
        self.gc_mydata = {}
        self.max_key_len = 15

        #for key in Stk.MY_DATA_KEYS:
        #    self.gc_mydata[key] = ""
        
        header_ok = False
        end_ok = False
        ss_ok = False
        
        # read the rest of the file
        count = 0
        for line in fi:
            count += 1
            line = line.strip()

            # ignore blank lines
            if( line == "" ):
                continue
            
            # get the first line
            if( line.startswith( "# STOCKHOLM 1.0" ) ):
                header_ok = True
                continue
            
            # ignore everything before the first line
            if( not header_ok ):
                continue
            
            # end file 
            if( line == "//" ):
                end_ok = True
                break
            
            # gets info from line
            data = line.split( )
            
            if( line.startswith( "#=GC SS_cons " ) ):
                self.secondary_structure += data[2]
                ss_ok = True
            elif( line.startswith( "#=GC RF" ) ):
                self.rf_structure += data[2]
            elif( line.startswith( "#=GC " ) ):
                if( len(data) == 3 ):
                    key = data[1]
                    self.gc_mydata[key] = self.gc_mydata.get( key, "" ) + data[2]
                else:
                    sys.stderr.write( "#=GC line was ignored: '%s'\n" %line )
            elif( not line.startswith( "#" ) ):
                if( len(data) > 2 ):
                    sys.stderr.write( "Error in line %d -> No spaces allowed in key or sequence.\n" %count )
                    return False
                else:
                    data[1] = data[1].upper()
                    
                    # check if the sequence has only the allowed chars
                    invalid_chars = [s for s in data[1] if not s in (Stk.CHAR_ALLOWED + self.CHAR_ALLOWED_EXTRA)]
                    if( len( invalid_chars ) > 0 ):
                        sys.stderr.write( "Warning (line %d) -> Invalid character(s): '%s'. Valid ones are: '%s' and '%s'.\n" %(count, ",".join(invalid_chars), Stk.CHAR_ALLOWED, self.CHAR_ALLOWED_EXTRA) )
                        sys.stderr.write( "Skipping line.\n" )
                        continue
                    
                    for c in Stk.GAP_NON_STANDARD:
                        data[1] = data[1].replace( c, Stk.GAP_STANDARD )
                    
                    # add sequence
                    self.add_seq( data[0], data[1] )

        # check balanced secondary structure
        stacks = [[] for x in Stk.OPEN_STRAND]
        count = 0
        for s in self.secondary_structure:
            count += 1
            
            i = Stk.OPEN_STRAND.find( s )
            if( i > -1 ):
                stacks[i].append( s )
                continue
            
            i = Stk.CLOSE_STRAND.find( s )
            if( i > -1 ):
                if( len(stacks[i]) == 0 ):
                    sys.stderr.write( "Error (sec. struct. pos %d) -> Unmatched closing char: %s\n" %(count, s) )
                    return False

                stacks[i].pop()
                continue
        
        for i in range(0, len(stacks)):
            if( len(stacks[i]) > 0):
                sys.stderr.write( "Error (sec. struct.) -> Not enough closing chars for '%s'.\n" %Stk.OPEN_STRAND[i] )
                return False

        if( not header_ok ):
            sys.stderr.write( "Warning -> '# STOCKHOLM 1.0' signature not found.\n" )
        
        if( not end_ok ):
            sys.stderr.write( "Warning -> End line string, '//', not found.\n" )
            
        if( not ss_ok ):
            sys.stderr.write( "Warning -> Secondary structure not found.\n" )
            
        return True

    def parse_ss(self):
        self.dstrands = []
        self.sstrands = []
        self.pairs = []

        stack = [[] for x in Stk.OPEN_STRAND]
        singles = []

        for (i, c) in enumerate( self.secondary_structure ):
            self.pairs.append( None )
            
            open_i = Stk.OPEN_STRAND.find( c )
            if( open_i > -1 ):
                stack[open_i].append( i )
                continue
            
            close_i = Stk.CLOSE_STRAND.find( c )
            if( close_i > -1 ):
                j = stack[close_i].pop()
                
                self.pairs[-1] = j 
                self.pairs[j] = i
                continue
        
        in_ss = False
        ss_first = None
        ss_last = None
        in_ds = False
        ds_first = (None, None)
        ds_last = (None, None)
         
        for (i, j) in enumerate(self.pairs):
            if( j is None ):
                if( in_ds ):
                    self.dstrands.append( (ds_first[0], ds_last[0], ds_last[1], ds_first[1]) )
                    in_ds = False
                    ds_first = (None, None)
                    ds_last = (None, None)
                
                if( not in_ss ):
                    in_ss = True
                    ss_first = i
                
                ss_last = i
            else:
                if( in_ss ):
                    self.sstrands.append( (ss_first, ss_last) )
                    in_ss = False
                    ss_first = None
                    ss_last = None
                    
                if( i < j ):
                    if( not in_ds ):
                        in_ds = True
                        ds_first = (i, j)
                
                    ds_last = (i, j)

        if( in_ds ):
            self.dstrands.append( (ds_first[0], ds_last[0], ds_last[1], ds_first[1]) )
        if( in_ss ):
            self.sstrands.append( (ss_first, ss_last) )


    def filter(self, limit=1.0):
        if( len(self.seqs) > 0 ):
            delete_list = []
            
            # add gaps to the end of all sequences if they are shorter than the secondary structure
            for (i, seq) in E_(self.seqs):
                diff = len(self.secondary_structure) - len(seq[1])
                if( diff > 0 ):
                    self.seqs[i] = (seq[0], seq[1] + (Stk.GAP_STANDARD * diff)) 
            
            for i in range(0, len(self.secondary_structure)):
                if( self.secondary_structure[i] in Stk.GAP_ALLOWED ):
                    gap_count = 0 
                    for seq in self.seqs:
                        if seq[1][i] in Stk.GAP_ALLOWED:
                            gap_count += 1
                    
                    gap_perc = float(gap_count)/float(len(self.seqs))
                    if( gap_perc >= limit ):
                        delete_list.append( i )
            
            if( len(delete_list) > 0 ):
                delete_list.sort(reverse=True)
                self.max_key_len = 15
                
                for i in delete_list:
                    self.secondary_structure = self.secondary_structure[:i] + self.secondary_structure[i+1:]
                    self.rf_structure = self.rf_structure[:i] + self.rf_structure[i+1:]
                    
                    for (key, value) in self.gc_mydata.items():
                        if( value != "" ):
                            self.gc_mydata[key] = value[:i] + value[i+1:]   
                    
                    for j in range( 0, len(self.seqs) ):
                        key = self.seqs[j][0]
                        seq = self.seqs[j][1] 
                        seq = seq[:i] + seq[i+1:]
                        
                        # update sequence and max_key
                        self.seqs[j] = (key, seq)
                        self.max_key_len = max( len( key ), self.max_key_len )
            
                self.parse_ss()

            # normalizes the keys (Rallee requirements)!
            for (i, seq) in E_(self.seqs):
                key_orig = seq[0]
                key_new = seq[0]
                
                for (o, n) in self.KEY_REPLACE_CHARS:
                    key_new = key_new.replace( o, n )
                
                if( key_new !=  key_orig ):
                    self.seqs[i] = (key_new, seq[1])
                    self.index[key_new] = self.index[key_orig]
                    del self.index[key_orig]
    
    def get_pair_index(self):
        return( self.pairs )
    
    def get_seqs(self):
        return( self.seqs )

    def get_seq(self, key):
        ndx = self.index.get( key, None )
        
        if( ndx is None ):
            return None
        else:
            return self.seqs[ndx][1]
    
    def get_secondary_structure(self):
        return( self.secondary_structure )

    def get_rf_structure(self):
        return( self.rf_structure )
    
    def get_strands(self):
        return( self.dstrands, self.sstrands )
    
class Clustal:
    NUCLEOTIDES = "ACGTU"
    GAP_STANDARD = "."
    GAP_NON_STANDARD = "-_:,"
    GAP_ALLOWED = GAP_STANDARD + GAP_NON_STANDARD
    CHAR_ALLOWED = NUCLEOTIDES + GAP_ALLOWED
    
    MAX_LOOP = 5

    def __init__(self, fi=None):
        self.max_key_len = 15
        
        if( (fi is None) or (not self.parse_file( fi )) ):
            self.seqs = []
            self.index = {}
        
    def add_seq(self, key, sequence):
        # if it's a new key
        if( not self.index.has_key( key ) ):
            # update the index
            self.index[key] = len(self.seqs)

            # insert the new sequence
            self.seqs.append( (key, sequence) )

            # update the max key length
            self.max_key_len = max( len(key), self.max_key_len )
        else:
            # if it's an existing key
            i = self.index[key]
            
            # update the sequence
            (k, s) = self.seqs[i]
            self.seqs[i] = (key, s + sequence)

    def del_seq(self, key):
        # if the key exists
        if( self.index.has_key( key ) ):
            i = self.index[key]
            
            # delete from the sequence list
            del self.seqs[i]
            
            # update the index
            self.index = {}
            self.max_key_len = 0
            for i in range( 0, len(self.seqs) ):
                (k, s) = self.seqs[i]
                self.index[k] = i
                self.max_key_len = max( len(k), self.max_key_len )

    def normalize_sequences(self, gap):
        for i in range( 0, len(self.seqs) ):
            (k, s) = self.seqs[i]
            
            for old_gap in self.GAP_ALLOWED:
                s = s.replace( old_gap, gap )
                
            self.seqs[i] = (k, s)
        
    def write_string(self):
        txt = "CLUSTAL X.X.X multiple sequence alignment\n\n\n"
        
        consensus = []
        for seq in self.seqs:
            if( len(consensus) == 0 ):
                consensus = [[] for i in range( 0, len(seq[1]) )] 
            
            txt += "%s %s\n" %(seq[0].ljust( self.max_key_len ), seq[1].upper())
            
            for i in range(0, len(seq[1])):
                c = seq[1].upper()[i]
                if( not c in consensus[i] ):
                    consensus[i].append( c )
        
        consensus = ""
        for i in range(0, len(consensus)):
            if( len(consensus[i]) == 1 ):
                consensus += "*"
            else:
                consensus += " "
        
        txt += "%s %s\n" %(" ".ljust( self.max_key_len ), consensus)

        return( txt )

    def write_file(self, fname):
        fo = open( fname, "w" )
        fo.write( self.write_string() )
        fo.close()
                
    def parse_file(self, fi):
        result = True
        clustal_ok = False
        
        self.seqs = []
        self.index = {}
        self.max_key_len = 15
        
        for line in fi:
            line = line.strip()
    
            if( line == "" ):
                continue
            elif( not clustal_ok and line.startswith( "CLUSTAL" ) ):
                clustal_ok = True
            elif( clustal_ok and line.find( "*" ) >= 0 ):
                continue
            elif( clustal_ok ):
                data = line.split()
                self.add_seq( data[0], data[1] )
            else:
                sys.stderr.write( "ERROR: Not a clustal file\n" )
                break
    
        return( result )

    def get_seqs(self):
        return( self.seqs )
    
    def write_stk(self, fname):
        stk = Stk()
        
        max_len = 0
        for seq in self.seqs:
            stk.add_seq( seq[0], seq[1] )
            max_len = max(max_len, len(seq[1]))
            
        stk.add_secondary_structure( "." * max_len)
        
        stk.write_file( fname )    
