import math
import os
import sys

import rnatools
import sequences
import candidates

class Outputs:
    def write_fold(fname, cands, seqs, group, constraint):
        rna = rnatools.RNATools()
        
        ps_2d = ""
        cands_group = []
        
        # for each specified group
        for (i, c) in enumerate(group):
            if( c >= len(cands) ):
                sys.stderr.write( "ERROR: Invalid candidate %d.\n\tRemeber: the list is 0-based.\n" )
                quit()
            
            cands_group.append( cands[c] )
            
            # get the corresponding sequence
            (start, end) = (cands[c].coords[0], cands[c].coords[1]+1)
            ndx = seqs.sequence_index[cands[c].key]
            seq = seqs.sequence_list[ndx].seq[start:end]
            
            if( i == 0 ):
                first_seq = seq
                first_len = len( seq )
            else:
                if( first_seq != seq ):
                    sys.stderr.write( "WARNING: sequence %d is different from the first sequence.\n" %i )
                if( first_len != len( seq ) ):
                    sys.stderr.write( "ERROR: secondary structures have different lengths! Can't merge\n" )
                    quit()
        
        pattern_2d_strict = candidates.Candidates.get_pattern_2d( cands_group, seqs, constraint )
        pattern_2d_loose = candidates.Candidates.get_pattern_2d( cands_group, seqs, constraint, strict=False )
        
        ps_2d = candidates.Candidates.get_postscript_2d( pattern_2d_loose )
        
        # get the fold forcing the mandatory base pairs
        fold = rna.fold( first_seq, pattern_2d_strict, constraint )
        rna.plot(seq, fold, fname, ps_2d)

    def write_alifold(fname, cands, seqs, group, constraint):
        rna = rnatools.RNATools()
        
        ps_2d = ""
        cands_group = []
        seqs_align = {}
        
        # for each specified group
        for (i, c) in enumerate(group):
            if( c >= len(cands) ):
                sys.stderr.write( "ERROR: Invalid candidate %d.\n\tRemeber: the list is 0-based.\n" )
                quit()
            
            cands_group.append( cands[c] )
            
            # get the corresponding sequence
            (start, end) = (cands[c].coords[0], cands[c].coords[1]+1)
            ndx = seqs.sequence_index[cands[c].key]
            seq = seqs.sequence_list[ndx].seq_gapped
            
            seqs_align["%s_%d" %(cands[c].key, i)] = seq
            
            if( i == 0 ):
                first_len = len( seq )
            else:
                if( first_len != len( seq ) ):
                    sys.stderr.write( "ERROR: sequences have different lengths! Can't merge\n" )
                    quit()
        
        pattern_2d_strict = candidates.Candidates.get_pattern_2d( cands_group, seqs, constraint, align=True )
        pattern_2d_loose = candidates.Candidates.get_pattern_2d( cands_group, seqs, constraint, align=True, strict=False )
        
        ps_2d = candidates.Candidates.get_postscript_2d( pattern_2d_loose )
        
        # get the fold forcing the mandatory base pairs
        fold = rna.alifold( seqs_align, pattern_2d_strict, constraint )
        rna.plot(seq, fold, fname, ps_2d)


    def write_info(fname, cands, seqs):
        # set header
        table = [["    #", "sequence", "model", "ver.", "score", "bpp", "pos", "cols", "motif"]]

        # set data
        for (i, cand) in enumerate(cands):
            row = []
            row.append( "%5d" %(i+1) ) 
            row.append( cand.key )
            row.append( cand.model.name )
            row.append( cand.model.version )
            row.append( "%6.2f" %cand.score )
            row.append( "%5.2f" %cand.bpp )
            row.append( "-".join( map( lambda x: "%d" %x[0], cand.chains_pos ) ) )
            row.append( "-".join( map( lambda x: "%d" %x[0], cand.chains_cols ) ) )
            row.append( "-".join( map( lambda x: "".join( x ), cand.chains_seqs ) ) )
            
            table.append( row )
       
        # compute column lengths
        col_len = [0] * len(table[0])
        
        for row in table:
            for (i, col) in enumerate(row):
                col_len[i] = max( col_len[i], len(col) )
                
        # print header
        fo = Outputs.__open_out_file( fname )
        fo.write( "\n%s\n" %("- " * 20 ) )
        fo.write( "Found %d candidates:\n\n" %(len(cands)) )
        
        for row in table:
            for (i, col) in enumerate(row):
                fo.write( "%s " %col.ljust( col_len[i] ) )
            fo.write( "\n" )
        
        fo.write( "%s\n" %("- " * 20 ) )
        fo.close()


    def write_seqs(fname, cands, seqs):
        max_key_len = reduce( lambda a, c: max( a, len(c.key) ), cands, 0 )

        # write the data
        fo = Outputs.__open_out_file( fname )

        for (i, cand) in enumerate(cands):
            # get the 2D pattern according to the mandatory base pairs
            pattern_2d = cand.get_pattern_2d()
            
            # get the corresponding sequence
            (start, end) = (cand.coords[0], cand.coords[1]+1)
            ndx = seqs.sequence_index[cand.key]
            seq = seqs.sequence_list[ndx].seq[start:end]
                
            fo.write( "%s\t%5d\t%s  %s\n" %( cand.model.full_name().ljust( 10 ), (i+1), cand.key.ljust( max_key_len ), seq ) )
            fo.write( "%s\t%s\t%s  %s\n" %( " " * 10, " " * 5, " " * max_key_len, pattern_2d ) )
            
        fo.close()

    def write_stk(fname, cands, seqs):
        max_key_len = reduce( lambda a, c: max( a, len(c.key) ), cands, 0 )
        patterns_2d= []

        # write the data
        fo = Outputs.__open_out_file( fname )
        fo.write( "# STOCKHOLM 1.0\n" )

        for (i, cand) in enumerate(cands):
            # get the corresponding sequence
            ndx = seqs.sequence_index[cand.key]
            seq = seqs.sequence_list[ndx].seq_gapped

            # get the 2D pattern according to the mandatory base pairs
            # append all until the end
            pat = cand.get_pattern_2d_align()
            pat += "." * (len(seq) - len(pat))
            
            txt = "#=GR %5d %s" %((i+1), cand.model.full_name())
                
            fo.write( "%s%s%s\n" %( cand.key, " " * (max_key_len - len(cand.key) + 2), seq ) )
            fo.write( "%s%s%s\n" %( txt, " " * (max_key_len - len(txt) + 2), pat ) )
        
        # merge structures
        pattern_2d = candidates.Candidates.get_pattern_2d( cands, seqs, align=True )
        
        txt = "#=GC SS_cons"
        fo.write( "%s%s%s\n" %( txt, " " * (max_key_len - len(txt) + 2), pattern_2d ) )
        fo.write( "//" )
        fo.close()

    def write_gff(fname, cands, seqs):
        # build the entries list
        entries = []
        for (i, cand) in enumerate(cands):
            (start, end) = (None, None)
            
            for chain_pos in cand.chains_pos:
                start = min( filter( lambda x: not x is None, chain_pos + [start] ) )
                end = max( chain_pos + [end] )
            
            sequence = seqs.sequence_list[ seqs.sequence_index[cand.key] ]
            
            if( sequence.rev ):
                start = sequence.seq_len - start - 1
                end = sequence.seq_len - end - 1
                (start, end) = (end, start)

            entries.append( (i, start, end, sequence.rev and "-" or "+", Outputs.__escape_sequence_key( sequence ) ) )
        
        # sort the entries by start coordinate
        entries.sort( cmp=lambda x, y:x[1] - y[1] )
            
        # write the data
        fo = Outputs.__open_out_file( fname )
        
        fo.write( "##gff-version   3\n" )
        for entry in entries:
            cand = cands[entry[0]]
            
            fo.write( "%s\tsmscan\tSM_%s\t%d\t%d\t%.5E\t%s\t.\t" %(entry[4], cand.model.name, entry[1], entry[2], cand.evalue, entry[3] ) )
            fo.write( "ID=sm_%06d;Name=%s_%s;Score=%.5E;BPP=%.5E;EValue=%.5E\n" %((entry[0]+1), cand.model.name, cand.model.version, cand.score, cand.bpp, cand.evalue) )
        
        fo.close()
    
    def __escape_sequence_key(sequence):
        key = list(sequence.key)
        
        if( sequence.rev ):
            key = key[:-1]

        for (i, c) in enumerate(key):
            if( not c.isalnum() and not c in ".:^*$@!+_?-|" ):
                key[i] = "_"
        
        return( "".join( key ) )
    
    def __open_out_file(out_name, attr="w"):
        if( not out_name is None ):
            fo = open( out_name, attr )
        else:
            fo = sys.stdout
        
        return( fo )

    write_info = staticmethod( write_info )
    write_seqs = staticmethod( write_seqs )
    write_stk = staticmethod( write_stk )
    write_fold = staticmethod( write_fold )
    write_alifold = staticmethod( write_alifold )
    write_gff = staticmethod( write_gff )
    
    __escape_sequence_key = staticmethod( __escape_sequence_key )
    __open_out_file = staticmethod( __open_out_file )
