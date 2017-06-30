# changed in version 0.0.2

import random
import os
import re
import sys

import interop
import file_io

from folds import *
from bpprobs import *


# TODO: transform all RNAtools' methods in static

class RNATools:
    PS_CMD_TAG = "% ADD COMMANDS HERE"
    
    def __init__(self):
        # decides in which system to operate
        if( interop.is_unix() ):
            self.os_cmd_rm = "rm"
            self.os_cmd_cp = "cp"
        if( interop.is_windows() ):
            self.os_cmd_rm = "del"
            self.os_cmd_cp = "copy"
        
    def bpprobs(self, sequence, constraint=""):
        key = self.get_sequence_key( sequence )
        
        ifile = self.set_input( key, [sequence, constraint] )
        ofile = self.set_output( key )
        bpfile = self.set_file_name( "dot", "ps" )
        
        opt = (constraint != "") and " -C" or ""
        
        cmd = "RNAfold -noPS %s -p < %s > %s" %(opt, ifile, ofile)
        os.system( cmd )
        
        bpprobs = self.parse_bpprobs( bpfile )

        self.clean_up( ifile, ofile, bpfile )
        
        return( bpprobs )

    def fold(self, sequence, constraint="", ensemble=False):
        key = self.get_sequence_key( sequence )
        
        ifile = self.set_input( key, [sequence, constraint] )
        ofile = self.set_output( key )
        bpfile = ""
        
        if( ensemble ):
            bpfile = self.set_file_name( "dot", "ps" )
        
        opt = (constraint != "") and " -C" or ""
        opt += ensemble and " -p" or ""
        
        cmd = "RNAfold -noPS %s < %s > %s" %(opt, ifile, ofile)
        os.system( cmd )
        
        fold = self.parse_fold( ofile, ensemble )

        self.clean_up( ifile, ofile, bpfile )
        
        return( fold )

    def alifold(self, sequences, constraint="", ensemble=False):
        key = self.get_sequence_key( "alitemp" )
        
        # build clustal
        clustal = file_io.Clustal()
        for (sk, ss) in sequences.items():
            clustal.add_seq( sk, ss )
        
        ifile = self.set_input( key, [constraint, clustal.write_string()] )
        ofile = self.set_output( key )
        bpfile = ""
        
        if( ensemble ):
            bpfile = self.set_file_name( "dot", "ps" )
        
        opt = (constraint != "") and " -C" or ""
        opt += ensemble and " -p" or ""
        
        cmd = "RNAalifold -noPS %s < %s > %s" %(opt, ifile, ofile)
        os.system( cmd )
        
        fold = self.parse_fold( ofile, ensemble )

        self.clean_up( ifile, ofile, bpfile )
        
        return( fold )
              
    def subopt(self, sequence, num, temp):
        key = self.get_sequence_key( sequence )
        
        ifile = self.set_input( key, [sequence] )
        ofile = self.set_output( key )
        
        cmd = "RNAsubopt -p %d -T %d < %s > %s" %(num, temp, ifile, ofile)
        os.system( cmd )
        
        folds = self.parse_subopt( ofile )
        
        self.clean_up( ifile, ofile )

        return( folds )

    def eval(self, sequence, folds):
        key = self.get_sequence_key( sequence )
        
        lines = []
        for fold in folds:
            lines += [sequence, fold.struct] 
        
        ifile = self.set_input( key, lines )
        ofile = self.set_output( key )
        
        cmd = "RNAeval < %s > %s" %(ifile, ofile)
        os.system( cmd )
        
        raw_folds = self.parse_eval( ofile )
        
        # update the folds with the MFEs
        for (fold, raw_fold) in zip(folds, raw_folds):
            # just check if we are matching the right fold
            if( fold.struct == raw_fold[0] ):
                fold.mfe = raw_fold[1]
        
        self.clean_up( ifile, ofile )

    def plot(self, sequence, fold, fname, ps_cmds=""):
        key = self.get_sequence_key( sequence )
        
        ifile = self.set_input( key, [sequence, fold.struct] )
        ofile = self.set_file_name( "rna", "ps" )

        opt = (ps_cmds != "") and "--pre \"%s\"" %RNATools.PS_CMD_TAG or ""

        cmd = "RNAplot %s < %s" %(opt, ifile)
        os.system( cmd )
        
        if( ps_cmds != "" ):
            self.plot_replace( ofile, ps_cmds )

        cmd = "ps2pdf %s %s" %(ofile, fname)
        res = os.system( cmd )
        
        if( res != 0 ):
            fname_alt = fname.replace( ".pdf", ".ps" )
            
            cmd = "%s %s %s" %(self.os_cmd_cp, ofile, fname_alt)
            os.system( cmd )
            
            sys.stdout.write( "WARNING: ps2pdf program not found!\n\tPostscript file '%s' was saved instead of '%s'\n" %(fname_alt, fname) )

        self.clean_up( ifile, ofile )
    
    # changed in version 0.0.2
    def get_sequence_key(self, sequence):
        return( "%s%x" %(sequence[:10], random.getrandbits(48)) )
    
    def set_file_name(self, name, ext):
        return( "./%s.%s" %(name, ext) )
    
    def set_input(self, key, lines):
        ifile = self.set_file_name( key, "input" )
        
        fi = open( ifile, "w" )
        fi.write( "\n".join( lines ) )
        fi.close()
        
        return( ifile )

    def set_output(self, key):
        return( self.set_file_name( key, "output" ) )
    
    def clean_up(self, *files):
        for f in files:
            if( (f != "") and os.path.isfile( f ) ):
                # deletes the file
                cmd = "%s %s" %(self.os_cmd_rm, os.path.normpath( f ))
                os.system( cmd )

    def parse_fold(self, ofile, ensemble=False):
        fo = open( ofile )
        
        # read, and ignore, the sequence line
        dummy = fo.readline()   

        # parse the secondary structure line
        (struct, mfe) = self.parse_line( fo.readline().strip() )
        
        if( ensemble ):
            # read, and ignore, the second structure
            (dummy, mfe) = self.parse_line( fo.readline().strip(), ensemble )
                 
        fo.close()
        
        return( Fold( struct, mfe ) )
    
    def parse_bpprobs(self, bpfile):
        bpprobs = BPProbs()
        
        fi = open( bpfile )
        start = False
        
        for line in fi:
            line = line.strip()
            
            if( not start ):
                start = (line == "%data starts here")
            elif( line == "showpage" ):
                break
            elif( "ubox" in line ):
                dt = line.split()
                
                bp1 = int(dt[0]) - 1
                bp2 = int(dt[1]) - 1
                
                bpprobs.add( bp1, bp2, float(dt[2])**2.0 )
        
        # changed in version 0.0.2
        # if start is still false then we have a problem!
        if( not start ):
            sys.stderr.write( "\nERROR: Can't parse 'rna.ps'. Check RNAFold version please.\n" )
            quit()
        
        fi.close()
        
        return( bpprobs )

    def parse_subopt(self, ofile):
        folds = []
        
        fo = open( ofile )

        for line in fo:
            (struct, mfe) = self.parse_line( line.strip() )
            folds.append( Fold( struct, mfe ) )
        
        fo.close()
        
        return( folds )

    def parse_eval(self, ofile):
        result = []
        expect_sequence = True 
        
        fo = open( ofile )
        
        for line in fo:
            if( line.startswith( "WARNING:" ) ):
                print line
            elif( expect_sequence ):
                # ignore, the sequence line
                expect_sequence = False
            else:
                # next time we wait for a sequence
                expect_sequence = True
        
                # read the secondary structure line with the MFE at the end
                (mfe, struct) = self.parse_line( line )
                result.append( (mfe, struct) )
                
        fo.close()
                
        return( result )
    
    def parse_line(self, line, ensemble=False):
        if( ensemble ):
            pat = "([\.\(\)\{\}\,\|]+)( \[\s*(\-?\d+\.\d+)\])"
        else:
            pat = "([\.\(\)]+)( \(\s*(\-?\d+\.\d+)\))?"
        
        match = re.match( pat, line )
        
        if( match != None ):
            struct = match.groups()[0]
            
            if( not match.groups()[2] is None ):
                mfe = float(match.groups()[2])
            else:
                mfe = 0.0
        else:
            print "ERROR in RNAeval output:\n'%s'" %line
            quit()
                
        return( struct, mfe )
        
    def plot_replace( self, fname, ps_cmds ):
        fi = open( fname )
        txt = fi.read().replace( RNATools.PS_CMD_TAG, ps_cmds )
        fi.close()
        
        fo = open( fname, "w" )
        fo.write( txt )
        fo.close()
        
        
if __name__ == '__main__':
    seq  = "ACGUACGUAAAAAACGUACGU"
    seqa = "AGGAACGUAAAAAACGUACCU"
    seqb = "AGGUUCGUAAAAGUCGUAUCU"
    seqc = "ACGUACGUCGCCGGCGUACGU"
    cnt  = "x(((.............)))x"
    
    aln = { '0':seq, 'A':seqa, 'B':seqb, 'C':seqc }
    
    rna = RNATools( )

    print "\nBase pair probabilities (w/o constraints):"
    bpps = rna.bpprobs( seq )
    print "\t0-20: ", bpps.prob_pair( 0, 20 )
    print "\t1-19: ", bpps.prob_pair( 1, 19 )
    print "\t1- 3: ", bpps.prob_pair( 1, 3 )
    
    print "\nBase pair probabilities (with constraints):"
    bpps = rna.bpprobs( seq, cnt )
    print "\t0-20: ", bpps.prob_pair( 0, 20 )
    print "\t1-19: ", bpps.prob_pair( 1, 19 )
    print "\t1- 3: ", bpps.prob_pair( 1, 3 )    
    
    print "\nFold (w/o constrains, w/o ensemble data):"
    fold = rna.fold( seq )
    print fold.struct

    print "\nFold (with constrains, w/o ensemble data):"
    fold = rna.fold( seq, cnt )
    print fold.struct
    print cnt

    print "\nFold (w/o constrains, with ensemble data):"
    fold = rna.fold( seq, ensemble=True )
    print fold.struct
    print fold.mfe

    print "\nFold (with constrains, with ensemble data):"
    fold = rna.fold( seq, cnt, ensemble=True )
    print fold.struct
    print cnt
    print fold.mfe

    """
    print "\nAlifold (w/o constrains, w/o ensemble data):"
    fold = rna.alifold( aln )
    print fold.struct

    print "\nAlifold (with constrains, w/o ensemble data):"
    fold = rna.alifold( aln, cnt )
    print fold.struct
    print cnt

    print "\nAlifold (w/o constrains, with ensemble data):"
    fold = rna.alifold( aln, ensemble=True )
    print fold.struct
    print fold.mfe

    print "\nAlifold (with constrains, with ensemble data):"
    fold = rna.alifold( aln, cnt, ensemble=True )
    print fold.struct
    print cnt
    print fold.mfe
    """
    
    print "\nSuboptimal folds:"
    folds = rna.subopt( seq, 5, temp=20.0 )
    for fold in folds:
        print fold.struct, fold.mfe
    
    print "\nEvaluate folds:"
    for s in aln.values():
        f = rna.fold( s )
        print s
        print "Predicted:", f.struct, f.mfe, 
        rna.eval( s, [f] )
        print "Evaluated:", f.mfe
        
    print "\nPlot"
    rna.plot( seq, rna.fold( seq ), "test.pdf")