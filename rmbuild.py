#!/usr/bin/env python

import optparse
import os
import re
import sys

import lib.bayes as bayes
import lib.def_parser as def_parser
import lib.file_io as file_io
import lib.model as model
import lib.model_parser as model_parser
import lib.node as node
import lib.pair as pair
import lib.pdb_utils as pdb_utils

#
# LOCAL FUNCTIONS
#
def get_columns_from_sources( ref_seqs, patterns, data_ref_seq ):
    columns = []
    
    for (ref_seq, pattern) in zip(ref_seqs, patterns):
        match = ""
        for (c, p) in zip(ref_seq, pattern):
            if( p == "*" ):
                match += (c in "ACGU" and c) or "."
                
        pos = data_ref_seq.find( match )
        
        if( pos < 0 ):
            return( None )
        
        for p in pattern:
            if( p == "*" ):
                columns.append( pos )
                pos += 1
            else:
                columns.append( None )
    
    return( columns )

def data_save( fname, data ):
    if( len(data) == 0 ):
        return
    
    fo = open( fname, "w" )
    
    
    fo.write( "key\tweight" )
    for (k, v) in data[0].items():
        if( k not in ["key", "weight"] ):
            fo.write( "\t%s" %k )
    
    fo.write( "\n" )
    
    for entry in data:
        fo.write( "%s\t%9.5f" %(entry["key"], entry["weight"]) )
        
        for (k, v) in entry.items():
            if( k not in ["key", "weight"] ):
                fo.write( "\t%s" %v )
        
        fo.write( "\n" )
    
    fo.close()

def analyse_pdb( pdb_file, coords ):
    pdb = pdb_utils.PDBStruct( )
    pdb.load( pdb_file, coords )
    
    ref_seq = pdb.get_sequences()
    interactions = pdb.get_interactions()
    
    return( ref_seq, interactions )

def find_seq( ref_seq, data_seq ):
    aux1 = ref_seq.replace( ".", "" ).replace( "-", "" )
    aux2 = data_seq.replace( ".", "" ).replace( "-", "" )
    
    pos = aux2.find( aux1 )
    
    for (i, c) in enumerate(data_seq):
        if( c not in ".-" ):
            if( pos == 0 ):
                return i
            else:
                pos -= 1
    
    return i 

def ref_seq_sync( ref_seqs, data_source_ref_seq ):
    # i - pointer for 'data_source_ref_seq'
    # j - pointer for ref_seqs
    ref_seqs_new = []
       
    for ref_seq in ref_seqs:
        ref_seq_new = ""
        
        i = find_seq( ref_seq, data_source_ref_seq )
        j = 0
        while( j < len(ref_seq) ):
            ci = data_source_ref_seq[i]
            cj = ref_seq[j]
            
            if( ci == cj ):
                ref_seq_new += ci
                (i, j) = (i+1, j+1)
            elif( ci in ".-" ):
                ref_seq_new += "."
                i += 1
            else:
                j += 1
        
        ref_seqs_new.append( ref_seq_new )
    
    return( ref_seqs_new ) 

def ref_seq_patterns( ref_seqs, data_source_ref_seq ):
    # i - pointer for 'data_source_ref_seq'
    # j - pointer for ref_seqs
    patterns = [] 
    
    for ref_seq in ref_seqs:
        pattern = ""
        
        i = find_seq( ref_seq, data_source_ref_seq )        
        j = 0
        while( j < len(ref_seq) ):
            ci = data_source_ref_seq[i]
            cj = ref_seq[j]
            
            if( ci == cj ):
                pattern += "*"
                (i, j) = (i+1, j+1)
            elif( cj in ".-" ):
                pattern += "-"
                j += 1
            else:
                print "ERROR: Error building patterns!"
                quit()
        
        patterns.append( pattern )
    
    return( patterns ) 

def translate_id( id, ref_seqs ):
    ref_seq = "".join( ref_seqs )
    count = id
    
    for (i, c) in enumerate(ref_seq):
        if( c not in ".-" ):
            if( count == 0 ):
                return i

            count -= 1
    
    print "ERROR: id %d not found in %s." %(id, ref_seq)
    quit() 

def check_def_model( fname, def_model ):
    # check if the ref seq exists
    if( len(def_model.ref_seqs) == 0 ):
        print "ERROR: REF_SEQ not found in the definition file '%s'" %fname
        quit()
    
    # check if the data sources exist
    if( len(def_model.data_sources) == 0 ):
        print "ERROR: data sources not found in the definition file '%s'" %fname
        quit()

    # check if the ref_seq length is equal to the number of nodes
    if( len("".join(def_model.ref_seqs)) != len(def_model.nodes) ):
        print "ERROR: %d nodes found %d expected." %(len(def_model.nodes), len("".join(def_model.ref_seqs)))
        print "       The number of nodes must be equal to REF_SEQS length: '%s'." %(", ".join(def_model.ref_seqs))
        quit()

    # check if each pattern has the same length of the respective ref_seq
    for data_source in def_model.data_sources:
        if( len("".join(def_model.ref_seqs)) != len("".join(data_source.patterns)) ):
            print "ERROR: The pattern '%s' and the REF_SEQ '%s' must have the same length." %("".join(def_source.patterns), "".join(def_model.ref_seqs))
            quit()

def check_ref_seqs( ref_seqs, data_source_ref_seq, data_source ):
    model_ref_seq = "".join( ref_seqs ).replace( ".", "" ).replace( "-", "" )
    data_source_ref_seq = data_source_ref_seq.replace( ".", "" ).replace( "-", "" )
    
    if( model_ref_seq != data_source_ref_seq ):
        print "ERROR: REF_SEQ in file '%s' must contain the declared REF_SEQ '%s'" %(data_source, model_ref_seq)

#
# *** MAIN ***
#
stk_files = []

#
# Read command line options
#
usage = "usage: %prog [options]"
version = "%prog 1.0"

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--name",   action="store", dest="name",     default=None, type="string", help="Module name and version: <NAME_V1.V2>." )
parser.add_option( "--pdb",    action="store", dest="pdb_file", default=None, type="string", help="PDB model." )
parser.add_option( "--coords", action="store", dest="coords",   default=None, type="string", help="Coodinates of the strands in the PDB model." )
parser.add_option( "--def",    action="store", dest="def_file", default=None, type="string", help="Module definition." )
parser.add_option( "--out",    action="store", dest="out_file", default=None, type="string", help="Model file." )

(options, args) = parser.parse_args()

if( (options.pdb_file is None) and (options.def_file is None) ):
    print "ERROR: Please define either PDB file or a definition file"
    quit()

def_file = options.def_file

# if we have a pdb file:
if( options.pdb_file is not None ):
    if( options.coords is None ):
        print "ERROR: Please define the coordinates of the strands in the pdb file"
        print "       CHAIN:POS:LEN,CHAIN:POS:LEN"
        quit()

    if( len(args) == 0 ):
        print "ERROR: Please define at least on STK file"
        quit()

    m = re.match( "^([A-Za-z]+)_([0-9]+\.[0-9]+)$", options.name )
    if( m is None ):
        print "WARNING: No (or invalid) model name and version."
        print "         Default model name and version used: NEWMODEL_1.0"
        
        name = "NEWMODEL"
        version = "1.0"
    else:
        name = m.groups()[0]
        version = m.groups()[1]
    
    # check the def_file
    if( def_file is None ):
        def_dir = os.path.dirname(options.pdb_file)
        def_file = "%s/%s_%s.def" %(def_dir != "" and def_dir or ".", name, version)
        
    # get the chains, sequences and the interactions from the pdb
    (ref_seqs, interactions) = analyse_pdb( options.pdb_file, options.coords )
    
    data_sources = []

    # for each STK file gets the reference sequence and build a super set of columns for each chain
    for stk_file in args:
        stk = file_io.Stk( open( stk_file ) )
        
        # get data ref seq
        data_source_ref_seq = stk.get_seq( "REF_SEQ" )
        check_ref_seqs( ref_seqs, data_source_ref_seq, stk_file )
        
        ref_seqs = ref_seq_sync( ref_seqs, data_source_ref_seq )
        
        data_sources.append( model.ModelDefinitionDataSource(stk_file, ref_seq=data_source_ref_seq) )
    
    # builds the patterns given the super ref_seq
    for data_source in data_sources:
        data_source.patterns = ref_seq_patterns( ref_seqs, data_source.ref_seq )
        
    # updates the interactions id's based on the patterns
    for inter in interactions:
        inter.id1 = translate_id( inter.id1, ref_seqs )
        inter.id2 = translate_id( inter.id2, ref_seqs )
    
    # setup the definition with the sequences, chains and interactions
    nodes = []
    pairing = []
    unpaired = []
    paired = set()
    symmetric = (len(ref_seqs) == 2) and (ref_seqs[0] == ref_seqs[1]) 

    for inter in filter( lambda x: x.wc, interactions ):
        pairing.append( pair.Pair(inter.id1, inter.id2, pair.Pair.PTYPE_MUST) )
        paired.add( inter.id1 )
        paired.add( inter.id2 )
    
    for (chain, ref_seq) in enumerate(ref_seqs):
        for (ndx, c) in enumerate(ref_seq):
            id = len(nodes)
            
            conds = map( lambda x: x.id1, filter( lambda x: x.id2 == id, interactions ) )
            nodes.append( node.Node( id, chain, ndx, conds ) )
            
            if( id not in paired ):
                unpaired.append( id )
    
    # save a definition file and set the definition input to that file
    model = model.Model( name, version, ref_seqs, data_sources, nodes, None, pairing, unpaired, None, sep_min=4, sep_max=0, symmetric=symmetric )
    model_parser.ModelParser().write( def_file, model )

# open the definition file
def_model = model_parser.ModelParser().parse( def_file )

check_def_model( def_file, def_model )

data = []

# for each source file
for def_data_source in def_model.data_sources:
    stk = file_io.Stk( open( def_data_source.align ) )
    
    # get data ref seq
    data_source_ref_seq = stk.get_seq( "REF_SEQ" )
    
    if( data_source_ref_seq is None ):
        print "ERROR: REF_SEQ not found in stk file '%s'" %def_data_source.align
        quit()
    
    # get the columns according to the defined chain pattern
    check_ref_seqs( def_model.ref_seqs, data_source_ref_seq, def_data_source.align )
    
    columns = get_columns_from_sources( def_model.ref_seqs, def_data_source.patterns, data_source_ref_seq )
    
    if( columns is None ):
        print "ERROR: Can't match the REF_SEQ from '%s' and the pattern from the definition file '%s'" %(def_data_source.align, def_file)
        quit()
    
    # for each sequence gets a "data" entry (contains an entry for each position of the pattern and the number of rows of the stk)
    weight = 100.0 / float(len(stk.seqs))
    for (key, seq) in stk.seqs:
        row = {'key': key, 'weight': weight}
        
        for (i, c) in enumerate( columns ):
            if( c is not None ):
                nt = seq[c]
            else:
                nt = "."
            
            row["N%d" %i] = nt
            
        data.append( row ) 

# save the data for training
data_file = def_file.replace( ".def", "" ) + ".data"
data_save( data_file, data )

# using the defined interactions build the Bayesian model and compute the probabilities
pmodel = bayes.Model( data )

for (i, node) in enumerate(def_model.nodes):
    # if the node has no dependency
    if( node.conds is None ):
        node.probs = def_model.prob_joint( pmodel, "N%d" %i )
    else:
        node.probs = def_model.prob_cond( pmodel, "N%d" %i, ["N%d" %i for i in node.conds] )

# save the model
if( options.out_file is None ):
    out_file = def_file.replace( ".def", "" ) + ".model"
    
    print "WARNING: No model file name was defined."
    print "         Default model file name used: '%s'" %out_file
else:
    out_file = options.out_file

def_model.data_sources = None
def_model.ref_seqs = None

model_parser.ModelParser().write( out_file, def_model )