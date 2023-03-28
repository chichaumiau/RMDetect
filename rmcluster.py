#!/usr/bin/env python
import sys
sys.path.append('./lib/')
import optparse
import os

import lib.cluster as cluster
import lib.candidates as candidates
import lib.config as config
import lib.interop as interop
import lib.results as results
import lib.sequences as sequences

   
#
# *** MAIN ***
#

config = config.Config()

#
# Read command line options
#
usage = "usage: %prog [options]"
version = "%prog " + config.VERSION

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--data-path", action="store",      dest="data_path", default=None, type="string", help="Directory containing all models." )

parser.add_option( "--in",        action="store",      dest="in_file",   default=None, type="string", help="Input file." )
parser.add_option( "--out",       action="store",      dest="out_file",  default=None, type="string", help="Output file/directory." )

parser.add_option( "--dlimit",    action="store",      dest="dlimit",    default=5,    type="int",    help="maximun column distance for two clusters to be merged." )

parser.add_option( "--out-dir",   action="store",      dest="out_dir",   default=None, type="string", help="Outputs each cluster in its own result file." )

parser.add_option( "--fig",       action="store_true", dest="fig",       default=False,               help="Display the cluster matrix." )
parser.add_option( "--fig-size",  action="store",      dest="fig_size",  default=50,                  help="Size of the cluster matrix." )

parser.add_option( "--min-count", action="store", dest="min_count",      default=3,    type="int",    help="Minimum absolute number of sequences in the cluster." )
parser.add_option( "--min-occur", action="store", dest="min_occur",      default=0.05, type="float",  help="Minimum percentage of cluster." )
parser.add_option( "--min-score", action="store", dest="min_score",      default=0.0,  type="float",  help="Minimum score to accept candidate." )
parser.add_option( "--min-bpp",   action="store", dest="min_bpp",        default=0.0,  type="float",  help="Minimum b.p. probability to accept candidate." )
parser.add_option( "--min-mi",    action="store", dest="min_mi",         default=0.0,  type="float",  help="Minimum b.p. probability to accept candidate." )
parser.add_option( "--min-h",     action="store", dest="min_h",          default=0.0,  type="float",  help="Minimum b.p. entropy to accept candidate." )

parser.add_option( "--filter-code", action="store", dest="filter_code",  default=None, type="string", help="Name of a file containing a filtering expression." )


(options, args) = parser.parse_args()

#
# Configure and initialize the search class
#
config.data_path = options.data_path
config.configure()

#
# Read the data from the results file
#
(seqs, cands, constraint, seq_count, tries) = results.Results.load( options.in_file, config.models )

#
# Output
#

# standard filter code
if( options.filter_code is None ):
    filter_code  = "ok = (count >= %d and " %options.min_count
    filter_code += "occur >= %f and " %options.min_occur
    filter_code += "score >= %f and " %options.min_score
    filter_code += "bpp >= %f and " %options.min_bpp
    filter_code += "mi >= %f and " %options.min_mi
    filter_code += "h >= %f)" %options.min_h
else:
    filter_code = open( options.filter_code ).read()
    
clust_algo = cluster.ClusterAlgorithm( cands, seq_count, seqs, filter_code ) 
clust_algo.go_cluster( options.dlimit )

# writes the summary of the results
clust_algo.write( options.out_file )

# writes one separate result file for each cluster
if( not options.out_dir is None ):
    out_dir = os.path.normpath( os.path.abspath( options.out_dir ) )

    if( not os.path.isdir( out_dir ) ):
        if( interop.is_unix() ):
            os.system( "mkdir %s" %out_dir )
        elif( interop.is_windows() ):
            os.system( 'mkdir "%s"' %out_dir )
        else:
            sys.stderr.write( "ERROR: Unknown system '%s'" %(os.name) )
    
    for (i, cluster) in enumerate(clust_algo.clusters):
        out_file = "%s/cluster_%03d_%s.res" %(out_dir, (i+1), cluster.model_name)
        results.Results.write( out_file, cluster.cands, seqs, "w", tries, constraint )

# plots the grafics
if( options.fig ):
    clust_algo.draw_matrix( options.fig_size )