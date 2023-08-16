#
# TODO: syntax verification of each line using regular expressions.
# TODO: at least one node at each strand should not allow a gap
#

import math
import os

from lib.model import *
from lib.node import *
from lib.pair import *

class ModelParser:
    def __init__(self):
        self.ptypes =  { "MUST":         Pair.PTYPE_MUST,
                         "NO":           Pair.PTYPE_NO,
                         "CAN":          Pair.PTYPE_CAN }
        
        self.ptypes_inv =  { Pair.PTYPE_MUST:      "MUST",
                             Pair.PTYPE_NO:        "NO",
                             Pair.PTYPE_CAN:       "CAN" }
    
    def parse(self, ifile):
        self.ifile = ifile
        
        self.name = ""
        self.version = ""
        self.ref_seqs = []
        self.data_sources = []
        self.nodes = []
        self.pairing = []
        self.unpaired = []
        self.strands = []
        self.order = None
        self.sep_min = None
        self.sep_max = None
        self.symmetric = False
        
        fi = open( ifile )
        
        count = 0
        for line in fi:
            count += 1
            line = line.strip()
            
            if( (line == "") or line.startswith( "#" ) ):
                continue

            if( line.startswith( "NAME" ) ):
                self.parse_name( line, count )
            elif( line.startswith( "REF_SEQS" ) ):
                self.parse_ref_seqs( line, count )
            elif( line.startswith( "DATA_SOURCE" ) ):
                self.parse_data_source( line, count )
            elif( line.startswith( "NODE" ) ):
                self.parse_node( line, count )
            elif( line.startswith( "PROB" ) ):
                self.parse_prob( line, count )
            elif( line.startswith( "PAIRS" ) ):
                self.parse_pairs( line, count )
            elif( line.startswith( "STRANDS" ) ):
                self.parse_strands( line, count )
            elif( line.startswith( "NO_PAIRING" ) ):
                self.parse_no_pairing( line, count )
            elif( line.startswith( "ORDER" ) ):
                self.parse_order( line, count )
            elif( line.startswith( "SEPARATION" ) ):
                self.parse_separation( line, count )
            elif( line.startswith( "SYMMETRIC" ) ):
                self.symmetric = True
        
        fi.close()
        
        return( Model( self.name, self.version, self.ref_seqs, self.data_sources, self.nodes, self.order, self.pairing, self.unpaired, self.strands, self.sep_min, self.sep_max, self.symmetric ) )

    def parse_name(self, line, count):
        data = line.split()

        if( len(data) == 3 ):
            self.name = data[1]
            self.version = data[2]
        else:
            print("ERROR: syntax error parsing NAME in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()

    def parse_ref_seqs(self, line, count):
        data = line.split()
        
        if( len(data) > 1 ):
            self.ref_seqs = data[1:]
        else:
            print("ERROR: syntax error parsing REF_SEQ in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()

    def parse_data_source(self, line, count):
        data = line.split()

        if( len(data) > 1 ):
            if( not os.path.isfile( data[-1] ) ):
                print("ERROR: alignment file %s not found in line %d" %(data[-1], count))
                print("       in file %s" %self.ifile)
                quit()
            
            self.data_sources.append( ModelDefinitionDataSource( data[-1], data[1:-1] ) )
        else:
            print("ERROR: syntax error parsing DATA_SOURCE in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()

    def parse_node(self, line, count):
        data = line.split()
        
        if( len(data) == 5 ):
            (id, chain, pos) = map(int, data[1:4])
            
            if( data[4] == "-" ):
                conds = None
            else:
                # ~conds = map( int, data[4].split( "," ) )
                conds = list(map( int, data[4].split( "," ) ))
            
            # enforces the order and sequentiality of node ids
            next_id = len(self.nodes)
            
            if( next_id == id ):
                self.nodes.append( Node( id, chain, pos, conds ) )
            else:
                print("ERROR: missing node id %d in line %d" %(next_id, count))
                print("       in file %s" %self.ifile)
                quit()
        else:
            print("ERROR: syntax error parsing NODE in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()
            
    def parse_prob(self, line, count):
        data = line.split()

        if( len(data) >= 2 ):
            id = int(data[1])
            
            dict = {}
            for dt in data[3:]:
                if( dt == "GC" ):
                    dict = "GC"
                    break
                else:
                    (k, v) = dt.split( ":" )

                    if( (v == "-") or (float(v) == 0.0) ):
                        continue
                    
                    # we only work with the logarithm of prob.
                    dict[k] = math.log(float(v))
                
            self.nodes[id].probs[data[2]] = dict
            
        else:
            print("ERROR: syntax error parsing PROB in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()
        
    def parse_pairs(self, line, count):
        data = line.split()

        if( len(data) == 4 ):
            (id1, id2) = map( int, data[1:3] )
            self.pairing.append( Pair( id1, id2, self.ptypes[data[3]] ) )
        else:
            print("ERROR: syntax error in parsing PAIRS line %d" %count)
            print("       in file %s" %self.ifile)
            quit()

    def parse_strands(self, line, count):
        data = line.split()

        if( (len(data) == 3) and ("-" in data[1]) and ("-" in data[2]) ):
            (A, B) = map( int, data[1].split( "-" ) )
            (C, D) = map( int, data[2].split( "-" ) )

            self.strands.append( Strands( A, B, C, D ) )
        else:
            print("ERROR: syntax error parsing STRANDS in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()

    def parse_order(self, line, count):
        data = line.split()

        if( len(data) == (len(self.nodes) + 1) ):
            self.order = map( int, data[1:] )
        else:
            print("ERROR: ORDER declaration has %d nodes. Expected %d nodes. In line %d" %(len(data)-1, len(self.nodes), count))
            print("       in file %s" %self.ifile)
            quit()

    def parse_no_pairing(self, line, count):
        data = line.split()

        for id in map(int, data[1:] ):
            self.unpaired.append( id )
    
    def parse_separation(self, line, count):
        data = line.split()

        if( len(data) == 3 ):
            (self.sep_min, self.sep_max) = (None, None)
            
            if( data[1] != "-" ):
                self.sep_min = int(data[1])
                
            if( data[2] != "-" ):
                self.sep_max = int(data[2])
        else:
            print("ERROR: syntax error parsing SEPARATION in line %d" %count)
            print("       in file %s" %self.ifile)
            quit()


    def write(self, ofile, model):
        self.fo = open( ofile, "w" )
        self.model = model
        
        self.__write_name( )
        self.__write_ref_seqs( )
        self.__write_data_sources( )
        self.__write_nodes( )
        self.__write_probs( )
        self.__write_pairs( )
        self.__write_no_pairing( )
        self.__write_order( )
        self.__write_separation( )
        
        self.fo.close()

    def __write_name(self):
        self.fo.write( "NAME\t%s\t%s\n\n" %(self.model.name, self.model.version) )
        
    def __write_ref_seqs(self):
        if( (self.model.ref_seqs is not None) and (self.model.ref_seqs != [] )):
            self.fo.write( "REF_SEQS" )
            
            for ref_seq in self.model.ref_seqs:
                self.fo.write( "\t%s" %(ref_seq) )
            
            self.fo.write( "\n\n" )

    def __write_data_sources(self):
        if( (self.model.data_sources is not None) and (self.model.data_sources != [] )):
            for data_source in self.model.data_sources:
                self.fo.write( "DATA_SOURCE" )
                self.fo.write( "\t%s" %"\t".join( (data_source.patterns) ) )
                self.fo.write( "\t%s\n" %data_source.align )
            
            self.fo.write( "\n" )
        
    def __write_nodes(self):
        for node in self.model.nodes:
            if( (node.conds is None) or (node.conds == []) ):
                conds = "-"
            else:
                conds = ",".join( map( str, node.conds ) )
            
            self.fo.write( "NODE\t%2d\t%d\t%d\t%s\n" %(node.id, node.chain, node.ndx, conds) )
        
        self.fo.write( "\n" )
            
    def __write_probs(self):
        for node in self.model.nodes:
            if( "*" in node.probs ):
                self.fo.write( "PROB\t%2d\t*" %(node.id) )
                self.__write_prob_line( node.probs["*"] )
            else:
                keys = node.probs.keys()
                sorted(keys)
                
                for k in keys:
                    self.fo.write( "PROB\t%2d\t%s" %(node.id, k) )
                    self.__write_prob_line( node.probs[k] )
            
            self.fo.write( "\n" )
            
    def __write_prob_line(self, probs):
        for k in "ACGU.":
            if( probs.__contains__( k ) ):
                s = "\t%s:%.3f" %(k, probs[k]) 
            else:
                s = "\t%s:-   " %(k)
                
            self.fo.write( s )
        self.fo.write( "\n" )
    
    def __write_pairs(self):
        for pair in self.model.pairing:
            self.fo.write( "PAIRS\t%2d\t%2d\t%s\n" %(pair.id1, pair.id2, self.ptypes_inv[pair.ptype]) )
        
        self.fo.write( "\n" )

    def __write_order(self):
        self.fo.write( "ORDER" )
        for node in self.model.order:
            self.fo.write( "\t%2d" %node.id )

        self.fo.write( "\n\n" )

    def __write_no_pairing(self):
        if( len(self.model.unpaired) > 0 ):
            self.fo.write( "NO_PAIRING" )
            for i in self.model.unpaired:
                self.fo.write( "\t%2d" %i )
    
            self.fo.write( "\n\n" )    
        
    def __write_separation(self):
        self.fo.write( "SEPARATION" )
        self.fo.write( "\t%s" %(self.model.sep_min is None and "-" or str(self.model.sep_min) ) )
        self.fo.write( "\t%s\n\n" %(self.model.sep_max is None and "-" or str(self.model.sep_max) ) )
