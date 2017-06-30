#!/usr/bin/env python

import optparse
import sys

import lib.candidates as candidates
import lib.config as config
import lib.outputs as outputs
import lib.results as results
import lib.sequences as sequences
    
#
# *** MAIN ***
#

config = config.Config()

#
# Read command line options
#
usage =  "usage: %prog [options] [fold | alifold | gff | seqs | align]\n"
usage += "Ex:\n"
usage += "\tOutputs a list of module candidates with the respective sequences and secondary structures:\n"
usage += "\t$ %prog seqs < xpto.res\n"
usage += "\t$ %prog seqs < xpto.res > xpto.txt\n"
usage += "\t$ %prog --in=xpto.res --out=xpto.txt seqs\n"
usage += "\t$ %prog --in=xpto.res seqs > xpto.txt\n"
usage += "\t$ %prog --out=xpto.txt seqs < xpto.res\n\n"

usage += "\tOutputs a secondary structure diagram with the indicated module candidates:\n"
usage += "\t%prog --in=xpto.res --out=xpto.pdf --cands='1 2 3 4 5' fold\n"
usage += "\t%prog --out=xpto.pdf --cands='1 2 3 4 5' fold < in.res\n\n"
usage += "\t%prog --out=xpto.pdf --cands=all fold < in.res\n\n"

usage += "\tOutputs a list of module candidates in the GFF format:\n"
usage += "\t$ %prog gff < xpto.res\n"
usage += "\t$ %prog gff < xpto.res > xpto.txt\n"
usage += "\t$ %prog --in=xpto.res --out=xpto.txt gff\n"
usage += "\t$ %prog --in=xpto.res gff > xpto.txt\n"
usage += "\t$ %prog --out=xpto.txt gff < xpto.res\n\n"


version = "%prog " + config.VERSION 

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "-d", "--data-path", action="store", dest="data_path", default=None, type="string", help="Directory containing all models." )
parser.add_option( "-i", "--in",        action="store", dest="in_file",   default=None, type="string", help="Input file." )
parser.add_option( "-o", "--out",       action="store", dest="out_file",  default=None, type="string", help="Output file." )
parser.add_option( "-c", "--cands",     action="store", dest="cands",     default="0",  type="string", help="List of candidates to show ('all'=for all)." )

(options, args) = parser.parse_args()

#
# Configure and initialize the search class
#
config.data_path = options.data_path
config.configure()

#
# Read the data from the results file
#
sys.stderr.write( "Reading data from results file: %s\n" %options.in_file )

(seqs, cands, constraint, count, tries) = results.Results.load( options.in_file, config.models )

#
# Output
#
if( len(args) == 0 ):
    outputs.Outputs.write_info( options.out_file, cands, seqs )
else:
    if( args[0] == "seqs" ):
        outputs.Outputs.write_seqs( options.out_file, cands, seqs )
    
    elif( args[0] == "align" ):
        outputs.Outputs.write_stk( options.out_file, cands, seqs )
        
    elif( args[0] in ("fold", "alifold") ):
        if( options.out_file is None ):
            sys.stderr.write( "ERROR: Specify the output PDF file, please.\n" )
            quit()
        
        group = []
        if( options.cands == "all" ):
            group = range( len( cands ) )
        else:
            for d in options.cands.split(" "):
                if( d.isdigit() ):
                    group.append( int(d)-1 )
        
        if( len(group) > 0 ):
            if( not options.out_file.endswith( ".pdf" ) ):
                sys.stderr.write( "\nWARNING: '.pdf' extension added to '%s'\n" %(options.out_file) )
                options.out_file += ".pdf"

            sys.stderr.write( "\nFolding %d candidates: '%s'\n" %(len(group), " ".join( map( lambda x: str(x+1), group ) ) ))
            
            if( args[0] == "fold" ):
                outputs.Outputs.write_fold( options.out_file, cands, seqs, group, constraint )
            else:
                outputs.Outputs.write_alifold( options.out_file, cands, seqs, group, constraint )
                
            sys.stderr.write( "\nDONE\n" )
        else:
            sys.stderr.write( "\nNo candidates to show\n" )
        
    elif( args[0] == "gff" ):
        outputs.Outputs.write_gff( options.out_file, cands, seqs )
    
