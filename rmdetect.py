#!/usr/bin/env python

#
# TODO: all internal class members should be named as '__<name>'
# TODO: the scan process should have a log.
#

import math
import optparse
import os
import sys
import time

import lib.candidates as candidates
import lib.config as config
import lib.results as results
import lib.scanner as scanner
import lib.sequences as sequences

last_cand = 0   

# mode can be one of: 'n', 'a', 'w' (see 'results.Results.write' for more information).
def cands_save( cands, tries, mode ):
    global options
    global seqs

    global scan_time_start
    global scan_time_end
    
    global scan_ctime_start
    global scan_ctime_end
    
    #
    # computes the joint bpp of the candidates
    #
    cands = candidates.Candidates.filter( cands, min_score=options.min_score ) #, min_evalue=options.min_evalue )
    cands = candidates.Candidates.compute_bpp( cands, seqs, config.constraint, cback_after_cand )
    cands = candidates.Candidates.filter( cands, min_bpp=options.min_bpp )
    cands.sort()
    
    if( mode == "a" ):
        comments = ""
    else:
        scan_time_end = time.time()
        scan_ctime_end = time.ctime()

        comments = "#- - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
        comments += "#\n"
        comments += "# script: %s\n" %(sys.argv[0])
        comments += "# parameters:\n"
        
        for s in sys.argv[1:]:
            comments += "#\t%s\n" %s
        
        comments += "#\n"
        comments += "# scan started at:     %s\n" %scan_ctime_start
        comments += "# scan ended at:       %s\n" %scan_ctime_end
        comments += "# total running time:  %d s\n" %(int(scan_time_end - scan_time_start))
        comments += "#- - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
        comments += "\n# meta data:\n"
        
    results.Results.write( options.out_file, cands, seqs, mode, tries, config.constraint, comments )

def cback_after_cand( n, total ):
    if( n > 0 ):
        back_str = "\b" * 52 
    else:
        back_str = ""
        
    sys.stderr.write( "%s Computing bpp: %6d of %6d%s" %(back_str, n+1, total, " " * 20 ))
    sys.stderr.flush()

    if( n+1 == total ):
        sys.stderr.write( "\n" )

        
# create the new (empty) result file
def cback_before_scan( total ):
    global options
    
    sys.stderr.write( "Scanning %d sequences...\n" %total )
    sys.stderr.flush()
    
    if( not options.out_file is None ):
        cands_save( [], tries=[0, 0], mode="n" )

# reset the candidate numbering for the next sequence
def cback_before_sequence( total, n, key ):
    global last_cand
    global seq_str
    
    seq_str = "Sequence %d of %d '%s'" %(n+1, total, key)
    sys.stderr.write( seq_str )
    sys.stderr.flush()
    
    last_cand = 0 

def cback_before_window( total, n ):
    global win_str
    
    win_str = "slide pos %d of %d%s" %(n+1, total, " " * (79 - len(win_str) - len(seq_str)))
    sys.stderr.write( win_str )
    sys.stderr.flush()
    
# append the last window results in to the file
def cback_after_window( cands, tries ):
    global last_cand
    global options
    global win_str
    
    sys.stderr.write( "%s" %( "\b" * len(win_str) ) )
    sys.stderr.flush()
    
    len_cands = len( cands )
    
    if( (len_cands - last_cand) > 1000 and (not options.out_file is None) ):
        cands_save( cands[last_cand:], tries, mode="a" )
        last_cand = len_cands
    
# write all data in to the result file (no append)
def cback_after_sequence( cands, tries ):
    global last_cand
    global options
    global seq_str
    
    sys.stderr.write( "%s" %( "\b" * len(seq_str) ) )
    sys.stderr.flush()

    if( not options.out_file is None ):
        cands_save( cands[last_cand:], tries, mode="a" )

# write all data in to the result file (no append)
def cback_after_scan( cands, tries ):
    global seqs
    
    cands_save( cands, tries, mode="w" )
    
#
# *** MAIN ***
#
config = config.Config()

scan_time_start = time.time()
scan_time_end = time.time()

scan_ctime_start = time.ctime()
scan_ctime_end = "still running..."

seq_str = ""
win_str = ""

#
# Read and parse the command line options
#
usage = "usage: %prog [options] [<in_file>]\n\t<in file> can be any fasta, stockholm, clustal or raw text file."
version = "%prog " + config.VERSION 

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--data-path",    action="store",      dest="data_path",    default=None,   type="string", help="Directory containing all models." )

parser.add_option( "--win-len",      action="store",      dest="win_len",      default=150,    type="int",    help="Window length for sliding window search." )
parser.add_option( "--win-step",     action="store",      dest="win_step",     default=None,   type="int",    help="Window step for sliding window search." )
parser.add_option( "--exc",          action="append",     dest="exclude",      default=[],                    help="List of models to exclude. By default exclude no model." )
parser.add_option( "--inc",          action="append",     dest="include",      default=[],                    help="List of models to include. By default include all models." )

parser.add_option( "--both-strands", action="store_true", dest="both_strands", default=False,                 help="Scans both strands of each sequence ." )

parser.add_option( "--out",          action="store",      dest="out_file",     default=None,   type="string", help="Output file. By default writes to standard output" )

parser.add_option( "--cutoff",       action="store",      dest="cutoff",       default=0.1,    type="float",  help="Cutoff value for halting search." )
parser.add_option( "--unknown-prob", action="store",      dest="unknown_prob", default=0.0001, type="float",  help="Minimum probability value for unknown nucleotides." )

parser.add_option( "--min-score",    action="store",      dest="min_score",    default=5.0,    type="float",  help="Minimum score to accept a candidate." )
parser.add_option( "--min-bpp",      action="store",      dest="min_bpp",      default=0.001,  type="float",  help="Minimum bpp to accept a candidate." )
#parser.add_option( "--min-evalue",   action="store",      dest="min_evalue",   default=None,   type="float",  help="Minimum -log(e-value) to accept candidate." )

parser.add_option( "--constraint",   action="store",      dest="constraint",   default=None,   type="string", help="File containing secondary structure constraints." )


(options, args) = parser.parse_args()

#
# Read the input sequences from wherever they came
#
seqs = sequences.Sequences( fname=((len(args) > 0) and args[0] or None), both_strands=options.both_strands )


#
# Configure and initialize the search class
#

# set the default values for 'win_len' and 'win_step'
config.win_len = options.win_len

if( options.win_step is None ):
    config.win_step = int(config.win_len / 2)
else:
    config.win_step = options.win_step

# Read the secondary structure constraints if they exist
if( not options.constraint is None ):
    if( os.path.isfile( options.constraint ) ):
        config.constraint = open( options.constraint ).read()
                    
        if( candidates.pattern_2_pairs( config.constraint ) == False ):
            sys.stderr.write( "ERROR: Unbalanced parenthesis in constraint!\n" )
            quit()
        
    else:
        sys.stderr.write( "WARNING: Constrain file '%s' not found!" %options.constraint )

config.data_path = options.data_path
config.exclude = options.exclude
config.include = options.include
config.cutoff = math.log(options.cutoff)
config.unknown_prob = math.log(options.unknown_prob)

config.configure()

# check the final window length 
if( config.constraint != "" and config.win_len < len(config.constraint) ):
    sys.stderr.write( "ERROR: Window length (%d) can't be smaller than constraint length (%d) !\n" %(config.win_len, len(config.constraint)) )
    quit()

#
# initializes the scanner object and configures all the callbacks
#
scnr = scanner.Scanner( config )

scnr.cback_before_scan = cback_before_scan
scnr.cback_before_sequence = cback_before_sequence
scnr.cback_before_window = cback_before_window
scnr.cback_after_window = cback_after_window
scnr.cback_after_sequence = cback_after_sequence
scnr.cback_after_scan = cback_after_scan

#
# Do the search
#
cands = scnr.scan( seqs )

sys.stderr.write( "\nDONE%s\n" %(" " * 80) )