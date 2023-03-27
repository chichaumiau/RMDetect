# recoded in version 0.0.3
# Now we use pipes to communicate with RNAfold software

import random
import os
import re
import sys

import interop
import file_io

import subprocess as sp

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
        # build the arguments
        #~args = ["--noPS", "-p"]
        args = ["-p"]
        
        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAfold", args=args, pre=[sequence], const=constraint, bpps=True )
        
        return( bpprobs )

    def fold(self, sequence, constraint="", ensemble=False):

        # build the arguments
        args = ["--noPS"]
        
        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAfold", args=args, pre=[sequence], const=constraint, ensemble=ensemble )

        # parse results
        fold = self.parse_fold( out, ensemble )

        return( fold )

    def alifold(self, sequences, constraint="", ensemble=False):
        # build the arguments
        args = ["--noPS"]
        
        # build clustal entries
        clustal = file_io.Clustal()
        
        for (sk, ss) in sequences.items():
            clustal.add_seq( sk, ss )
        
        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAalifold", args=args, const=constraint, post=[clustal.write_string()], ensemble=ensemble )
        
        # parse results
        fold = self.parse_fold( out, ensemble )

        return( fold )
              
    def subopt(self, sequence, num, temp):
        # build the arguments
        args = ["-p", str(num), "-T", str(temp)]

        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAsubopt", args=args, pre=[sequence] )
        
        # parse results
        folds = self.parse_subopt( out )

        return( folds )

    def eval(self, sequence, folds):
        # prepares input
        lines = []
        for fold in folds:
            lines += [sequence, fold.struct] 

        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAeval", [], lines )
        
        # parse output
        raw_folds = self.parse_eval( out )
        
        # update the folds with the MFEs
        for (fold, raw_fold) in zip(folds, raw_folds):
            # just check if we are matching the right fold
            if( fold.struct == raw_fold[0] ):
                fold.mfe = raw_fold[1]
        
    def plot(self, sequence, fold, fname, ps_cmds=""):
        rna_ps = "./rna.ps"
        
        # build the arguments
        args = ["--pre", "%s" %RNATools.PS_CMD_TAG]
        
        # call the tool
        (out, bpprobs) = self.vienna_tool( "RNAplot", args, [sequence, fold.struct] )
              
        # replace the commands in the "rna.ps" file
        self.plot_replace( rna_ps, ps_cmds )

        # convert from postscript to pdf 
        cmd = "ps2pdf %s %s" %(rna_ps, fname)
        res = os.system( cmd )
        
        if( res != 0 ):
            fname_alt = fname.replace( ".pdf", ".ps" )
            
            cmd = "%s %s %s" %(self.os_cmd_cp, rna_ps, fname_alt)
            os.system( cmd )
            
            sys.stdout.write( "WARNING: ps2pdf program not found!\n\tPostscript file '%s' was saved instead of '%s'\n" %(fname_alt, fname) )
    
    #
    # cmd: command to execute
    # args: arguments to use
    # lines: lines to use as input
    #
    def vienna_tool(self, cmd, args=[], pre=[], const="", post=[], ensemble=False, bpps=False):
        cmds = [cmd]
        cmds.extend( args )
        #~ cmds.extend( (const != "") and ["-C"] or [] )
        cmds.extend( ensemble and ["-p"] or [] )
        
        input = "\n".join( filter( lambda x: x != "", pre + [const] + post ) )

        #~sys.stderr.write( input )
        p = sp.Popen( cmds, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE )
        (out, err) = p.communicate( input.encode() )

        #~sys.stderr.write( out )
        #~sys.stderr.write( err )
        
        # if we produce a 'rna.ps' file parses and cleanup
        bpprobs = None
        if( bpps ):
            bpprobs = self.parse_bpprobs( "rna.ps" )
            self.clean_up( "rna.ps" )
        
        return( out, bpprobs )
       
    def clean_up(self, *files):
        for f in files:
            if( (f != "") and os.path.isfile( f ) ):
                # deletes the file
                cmd = "%s %s" %(self.os_cmd_rm, os.path.normpath( f ))
                os.system( cmd )

    #
    # parsing routines
    #
    def parse_fold(self, txt, ensemble=False):
        lines = txt.split( "\n")
        
        # parse the secondary structure line
        (struct, mfe) = self.parse_line( lines[1].strip() )
        
        if( ensemble ):
            # read, and ignore, the second structure
            (dummy, mfe) = self.parse_line( lines[2].strip(), ensemble )
        
        return( Fold( struct, mfe ) )
    
    def parse_bpprobs(self, bpfile):
        bpprobs = BPProbs()
        #~sys.stderr.write( bpfile)
        fi = open( bpfile )
        start = False
        
        for line in fi:
            line = line.strip()
            
            if( not start ):
                start = (line == "% data start here")
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

    def parse_subopt(self, out):
        folds = []
        
        for line in out.strip().split( "\n" ):
            (struct, mfe) = self.parse_line( line.strip() )
            folds.append( Fold( struct, mfe ) )
        
        return( folds )

    def parse_eval(self, out):
        result = []
        expect_sequence = True 
        
        for line in out.strip().split( "\n" ):
            if( line.startswith( "WARNING:" ) ):
                print(line)
            elif( expect_sequence ):
                # ignore, the sequence line
                expect_sequence = False
            else:
                # next time we wait for a sequence
                expect_sequence = True
        
                # read the secondary structure line with the MFE at the end
                (mfe, struct) = self.parse_line( line )
                result.append( (mfe, struct) )
                
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
            print("ERROR in output:\n\t'%s'" %line)
            quit()
                
        return( struct, mfe )
        
    def plot_replace( self, fname, ps_cmds="" ):
        if( ps_cmds != "" ):
            txt = open( fname ).read().replace( RNATools.PS_CMD_TAG, ps_cmds )
            
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

    print("\nBase pair probabilities (w/o constraints):")
    bpps = rna.bpprobs( seq )
    print("\t0-20: ", bpps.prob_pair( 0, 20 ) )
    print("\t1-19: ", bpps.prob_pair( 1, 19 ) )
    print("\t1- 3: ", bpps.prob_pair( 1, 3 ) )
    
    print("\nBase pair probabilities (with constraints):")
    bpps = rna.bpprobs( seq, cnt )
    print("\t0-20: ", bpps.prob_pair( 0, 20 ) )
    print("\t1-19: ", bpps.prob_pair( 1, 19 ) )
    print("\t1- 3: ", bpps.prob_pair( 1, 3 ) )    
    
    print("\nFold (w/o constrains, w/o ensemble data):")
    fold = rna.fold( seq )
    print(fold.struct)

    print("\nFold (with constrains, w/o ensemble data):")
    fold = rna.fold( seq, cnt )
    print(fold.struct)
    print(cnt)

    print("\nFold (w/o constrains, with ensemble data):")
    fold = rna.fold( seq, ensemble=True )
    print(fold.struct)
    print(fold.mfe)

    print("\nFold (with constrains, with ensemble data):")
    fold = rna.fold( seq, cnt, ensemble=True )
    print(fold.struct)
    print(cnt)
    print(fold.mfe)

    print("\nAlifold (w/o constrains, w/o ensemble data):")
    fold = rna.alifold( aln )
    print(fold.struct)

    print("\nAlifold (with constrains, w/o ensemble data):")
    fold = rna.alifold( aln, cnt )
    print(fold.struct)
    print(cnt)

    print("\nAlifold (w/o constrains, with ensemble data):")
    fold = rna.alifold( aln, ensemble=True )
    print(fold.struct)
    print(fold.mfe)

    print("\nAlifold (with constrains, with ensemble data):")
    fold = rna.alifold( aln, cnt, ensemble=True )
    print(fold.struct)
    print(cnt)
    print(fold.mfe)
    
    print("\nSuboptimal folds:")
    folds = rna.subopt( seq, 5, temp=20.0 )
    for fold in folds:
        print(fold.struct, fold.mfe)
    
    print("\nEvaluate folds:")
    for s in aln.values():
        f = rna.fold( s )
        print(s)
        print("Predicted:", f.struct, f.mfe, rna.eval( s, [f] ))
        
        print("Evaluated:", f.mfe)
        
    print("\nPlot")
    rna.plot( seq, rna.fold( seq ), "test.pdf")
