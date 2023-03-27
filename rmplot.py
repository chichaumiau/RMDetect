#!/usr/bin/env python

import optparse
# ~import psyco

import cmotif_scanner3.lib.candidates as candidates
import cmotif_scanner3.lib.config as config
import cmotif_scanner3.lib.analysis as analysis
import cmotif_scanner3.lib.results as results
import cmotif_scanner3.lib.sequences as sequences

    
#
# *** MAIN ***
#

# ~psyco.full()

#
# Read command line options
#
usage = "usage: %prog [options] <'sesp'/'sts'/'align'/'plot'/'data'> <model> [<column pairs> <positives>]"
version = "%prog 1.0"

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--data-path", action="store", dest="data_path",  default=None,     type="string", help="Directory containing all models." )

parser.add_option( "--in",        action="store", dest="in_file",    default=None,     type="string", help="Input file." )
parser.add_option( "--out",       action="store", dest="out_file",   default=None,     type="string", help="Output file/directory." )

parser.add_option( "--no-fig",    action="store_true", dest="no_fig", default=False                 , help="Don't show the picture." )

parser.add_option( "--x-feat",    action="store", dest="x_feature",  default="score",  type="string", help="'evalue', 'bpp', ('score')" )
parser.add_option( "--y-feat",    action="store", dest="y_feature",  default="bpp",    type="string", help="'evalue', ('bpp'), 'score'" )

parser.add_option( "--x-min",     action="store", dest="x_min",      default=None,     type="float",  help="Smalest 'x' value for a predicted positive candidate." )
parser.add_option( "--y-min",     action="store", dest="y_min",      default=None,     type="float",  help="Smalest 'y' value for a predicted positive candidate." )

parser.add_option( "--threshold", action="store", dest="threshold",  default=None,     type="float",  help="Threshold for selecting candidates (matrix option)" )
parser.add_option( "--steps",     action="store", dest="steps",      default=None,     type="int",    help="Parameters step (sts option)" )
parser.add_option( "--signal",    action="store", dest="signal",     default=1,        type="int",    help="Area to consider the positives 1/-1 (sts option)" )
parser.add_option( "--stats",     action="store", dest="stats",      default="mcc",    type="string", help="Statistics to analyse 'mcc', 'tpr', 'fpr' (sts option)" )


(options, args) = parser.parse_args()

if( len(args) < 2 ):
    parser.print_usage()
    quit()

analysis_type = args[0].lower()
model_name = args[1].upper()

#
# Configure and initialize the search class
#
config = config.Config()
config.data_path = options.data_path
config.configure()

#
# Read the data from the results file
#
(seqs, cands, constraint, seq_count, tries) = results.Results.load( options.in_file, config.models )

#
# Output
#
if( analysis_type in ['sesp', 'data', 'sts'] ):
    if( len(args) < 4 ):
        print "Please indicate the positive columns (col_a1-col_b1,col_a2-col_b2,...,col_an-col_bn) and number of sequences"
        quit()

    cols = []
    col_pairs = args[2].split( "," )
    for col_pair in col_pairs:
        cols.append( map( int, col_pair.split( "-" ) ) )
        
    positives = int( args[3] ) * len(col_pairs)

    print "-----------------------------------"
    print "Results for %s:" %options.in_file
    print "\nParameters:"
    print "\tsignal == %d" %options.signal
    print "\tsteps == ", options.steps
    print "\tx_min == ", options.x_min
    print "\ty_min == ", options.y_min
    
    print ""

    if( analysis_type == 'sesp' ):
        analysis.Analysis.sens_spec( cands, model_name, cols, options.x_feature, options.y_feature, positives, tries, options.signal, options.x_min, options.y_min, options.out_file, options.no_fig )
    elif( analysis_type == "data" ):
        analysis.Analysis.data( cands, model_name, cols, options.x_feature, options.y_feature, positives )
    else:
        if( not options.stats in ["mcc", "tpr", "fpr"] ):
            print "ERROR: Invalid statistic option '%s'." %options.stats
            quit()
            
        analysis.Analysis.statistics( cands, model_name, cols, options.x_feature, options.y_feature, positives, tries, options.signal, options.stats, options.steps, options.out_file, options.no_fig )
    
    print "-----------------------------------"
             
elif( analysis_type == "align" ):
    cands = analysis.Analysis.matrix( cands, model_name, seq_count, options.x_feature, options.y_feature, options.x_min, options.y_min, options.threshold, options.out_file, options.no_fig )
    
elif( analysis_type == "plot" ):
    analysis.Analysis.plot( cands, model_name, options.x_feature, options.y_feature, options.out_file, options.no_fig )
    
else:
    print "Unknown type of analysis: '%s'" %analysis_type
    quit()
