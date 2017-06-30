#!/usr/bin/env python

import optparse
import sys

import lib.candidates as candidates
import lib.config as config
import lib.results as results
import lib.sequences as sequences

# reset the candidate numbering for the next sequence
def cback_after_cand( n, total ):
    if( n % 100 == 0 ):
        sys.stderr.write( "%sFiltering: %d of %d" %("\b" * 50, n, total) )
        sys.stderr.flush()
    
#
# *** MAIN ***
#
config = config.Config()

#
# Read command line options
#
usage = "usage: %prog [options]"
version = "%prog 1.0"

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--data-path",     action="store", dest="data_path",     default=None, type="string", help="Directory containing all models." )

parser.add_option( "--in",            action="store", dest="in_file",       default=None, type="string", help="Input file." )
parser.add_option( "--out",           action="store", dest="out_file",      default=None, type="string", help="Output file." )

parser.add_option( "--min-score",     action="store", dest="min_score",     default=None, type="float",  help="Minimum score to accept candidate." )
parser.add_option( "--min-bpp",       action="store", dest="min_bpp",       default=None, type="float",  help="Minimum b.p. probability to accept candidate." )

parser.add_option( "--model",         action="store", dest="model",         default=None, type="string", help="Retain only the candidates of the specified model." )
parser.add_option( "--top-seq",       action="store", dest="top_seq",       default=None, type="int",    help="For each sequence retain the top N candidates, all models confounded." )
parser.add_option( "--top-model",     action="store", dest="top_model",     default=None, type="int",    help="For each sequence retain the top N candidates, all models confounded." )
parser.add_option( "--top-seq-model", action="store", dest="top_seq_model", default=None, type="int",    help="For each sequence retain the top N candidates, all models confounded." )
parser.add_option( "--no-overlap",       action="store_true", dest="no_overlap",       default=False,          help="Remove overlapping candidates" )


(options, args) = parser.parse_args()

#
# Configure and initialize the search class
#
config.data_path = options.data_path
config.configure()

#
# Read the data from the results file
#
(seqs, cands, constraint, count, tries) = results.Results.load( options.in_file, config.models )

#
# Do the filtering
#
if( (not options.min_score is None) or (not options.min_bpp is None) ):
    cands = candidates.Candidates.filter( cands, min_score=options.min_score, min_bpp=options.min_bpp, cback_after_cand=cback_after_cand )
    sys.stderr.write( "%sDONE: 'score' and 'bpp' filtering%s\n" %("\b" * 50, " " * 50) )

if( not options.model is None ):
    cands = candidates.Candidates.get_model( cands, options.model, cback_after_cand )
    sys.stderr.write( "%sDONE: model filtering%s\n" %("\b" * 50, " " * 50) )

if( not options.top_seq is None ):
    cands = candidates.Candidates.top_seq( cands, options.top_seq, cback_after_cand )
    sys.stderr.write( "%sDONE: 'top-seq' filtering%s\n" %("\b" * 50, " " * 50) )

if( not options.top_model is None ):
    cands = candidates.Candidates.top_model( cands, options.top_model, cback_after_cand )
    sys.stderr.write( "%sDONE: 'top-model' filtering%s\n" %("\b" * 50, " " * 50) )

if( not options.top_seq_model is None ):
    cands = candidates.Candidates.top_seq_model( cands, options.top_seq_model, cback_after_cand )
    sys.stderr.write( "%sDONE: 'top-seq-model' filtering%s\n" %("\b" * 50, " " * 50) )

if( options.no_overlap ):
    cands = candidates.Candidates.nooverlap( cands, cback_after_cand )
    sys.stderr.write( "%sDONE: no overlap filtering%s\n" %("\b" * 50, " " * 50) )

results.Results.write( options.out_file, cands, seqs, "w", tries, config.constraint )
