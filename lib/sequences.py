import os
import sys

import interop
import file_io

class Sequences:
    #
    # Load all the sequences either from a file, stdin or a string
    #
    def __init__(self, fname=None, seq=None, both_strands=False):
        if( interop.is_unix() ):
            self.os_end_key = "^D"
        
        if( interop.is_windows() ):
            self.os_end_key = "^Z"

        self.fname = fname
        self.from_stdin = False
        self.alignment = False
        self.both_strands = both_strands
        
        self.sequence_list = []
        self.sequence_max_pos = 0
        self.sequence_max_col = 0
        self.sequence_index = {}
    
        if( not seq is None ):
            self.__load_str( seq )
            self.from_stdin = True
        elif( not fname is None ):
            if( os.path.isfile( fname ) ):
                fi = open( fname )
                
                if( fname.endswith( ".seq" ) ):
                    self.__load_seq( fi )
                    
                elif( fname.endswith( ".fasta" ) ):
                    self.__load_fasta( fi )
                    
                elif( fname.endswith( ".stk" ) or fname.endswith( ".stockholm" ) ):
                    self.__load_stk( fi )
                    self.alignment = True
    
                elif( fname.endswith( ".aln" ) ):
                    self.__load_clustal( fi )
                    self.alignment = True
            else:
                sys.stderr.write( "File not found: %s\n" %fname )
                quit()
        else:
            sys.stderr.write( "Type an RNA sequence (end it with %s):\n" %self.os_end_key )
            
            self.__load_str( sys.stdin.read() )
            self.from_stdin = True        
        
        # TODO: check all sequences
        for seq in self.sequence_list:
            self.sequence_max_pos = max( self.sequence_max_pos, len(seq.seq) )
            self.sequence_max_col = max( [self.sequence_max_col, len(seq.seq_gapped), len(seq.seq)] )
            
        if( self.both_strands ):
            self.__reverse_complement()

    def seq_count(self):
        if( self.both_strands ):
            return( len(self.sequence_list) * 2 )
        else:
            return( len(self.sequence_list) )
        
    def get_sequence(self, key):
        return( self.sequence_list[self.sequence_index[key]] )

    def __load_str(self, seq):
        self.sequence_index["unique"] = 0
        self.sequence_list.append( Sequence( "unique", seq ) )

    def __load_seq(self, fi):
        self.sequence_index["unique"] = 0
        self.sequence_list.append( Sequence( "unique", fi.read() ) )

    def __load_fasta(self, fi):
        fasta = file_io.Fasta( fi )

        for s in fasta.seqs:
            self.sequence_index[s[0]] = len( self.sequence_list )
            self.sequence_list.append( Sequence( s[0], s[1] ) )

    def __load_stk(self, fi):
        stk = file_io.Stk( fi )

        for s in stk.seqs:
            self.sequence_index[s[0]] = len( self.sequence_list )
            self.sequence_list.append( Sequence( s[0], s[1] ) )

    def __load_clustal(self, fi):
        clustal = file_io.Clustal( fi )

        for s in clustal.seqs:
            self.sequence_index[s[0]] = len( self.sequence_list )
            self.sequence_list.append( Sequence( s[0], s[1] ) )
    
    def __reverse_complement(self):
        reverse_sequence_list = []
        
        for sequence in self.sequence_list:
            reverse_sequence_list.append( sequence.reverse() )
            
        for rev_sequence in reverse_sequence_list:
            self.sequence_index[rev_sequence.key] = len(self.sequence_list)
            self.sequence_list.append( rev_sequence )


class Sequence:
    VALID_GAPS = ".-_"

    def __init__(self, key="", seq="", rev=False):
        self.key = key
        self.seq = self.__valid_sequence( seq )
        self.seq_gapped = self.seq
        
        self.rev = rev
        self.col_index = None
               
        if( not self.seq is None ):
            gapped = False
            
            for gap in self.VALID_GAPS:
                if( gap in self.seq ):
                    gapped = True
                    break
            
            if( gapped ):
                self.__build_col_index()
        
        self.seq_len = len(self.seq)
                
    def reverse(self):
        new_seq = self.seq[::-1]
        for (c1, c2) in ( ("A", "u"), ("C", "g"), ("G", "c"), ("U", "a") ):
            new_seq = new_seq.replace( c1, c2 )
        
        return( Sequence( self.key + "--", new_seq, True ) )
    
    def __valid_sequence(self, seq):
        # validate sequence
        #    - check invalid nucleotides
        seq = seq.strip().replace( "\n", "" ).upper().replace( "T", "U" )
 
        return( seq )

    def __valid_key_DELETE(self, key):
        rep = "/?*."

        for c in rep:
            key = key.replace( c, "_" )
        
        return( key )
            
    def __build_col_index(self):
        self.seq = ""
        self.col_index = []
        
        for (n, c) in enumerate(self.seq_gapped):
            if( not c in self.VALID_GAPS ):
                self.seq += c
                self.col_index.append( n )
