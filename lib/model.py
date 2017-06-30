import copy
import math
import sys

import bayes
#import model_parser

from node import *
from pair import *


class Model:
    def __init__(self, name, version, ref_seqs, data_sources, nodes, order, pairing, unpaired, strands, sep_min=0, sep_max=0, symmetric=False):
        # standard model definitions
        self.name = name
        self.version = version
        self.ref_seqs = ref_seqs
        self.data_sources = data_sources
        self.nodes = nodes
        self.pairing = pairing
        self.unpaired = unpaired
        self.strands = strands
        self.sep_min = sep_min
        self.sep_max = sep_max
        self.symmetric = symmetric
        
        if( order is None ):
            self.__init_order()
        else:
            self.order = []
            for id in order:
                self.order.append( filter( lambda n: n.id == id, self.nodes )[0] )
        
        # TODO: remove this line!!
        #self.__init_order()

        # initialize derived variables
        self.chains_length = {}
        self.__init_chains()

        # initialize search structures and data
        self.OFFSETS = []      # offset matrix
        self.precomp = []      # positions that can take advantage of pre computed probabilities
        self.count = 0
        self.pairs_list = []
        
        self.__init_offsets()
        self.__init_pairing_list()
        self.__init_search_nodes()
        
        # TODO: validate if all probabilities are normalized
        # TODO: validate that the id and the position of the node in the list are the same
    
    def full_name(self):
        fn = "%s_%s" %(self.name, self.version)
        
        return( fn.upper() )

    def explode(self, config, n=None, seq=[], prob=0.0):
        limit = math.log( 0.25 ) * len(self.nodes)
        
        if( prob < limit ):
            return
        
        if( n is None ):
            n = 0
            seq = [None] * len(self.nodes)
            
            self.explode( config, n, seq, prob )
        elif( n == len(self.order) ):
            print "".join( seq ),  (prob - limit) / math.log(2.0)
        else:
            node = self.order[n]
            
            if( node.conds is None ):
                prob_key = "*"
            else:
                prob_key = ""
                    
                for cond in node.conds:
                    prob_key += seq[cond]
                
            prob_dict = node.probs.get( prob_key, {} )
                
            for nt in "ACGU.":
                if( prob_dict == "GC" ):
                    prob_node = math.log( 0.25 )
                else:
                    prob_node = prob_dict.get( nt, config.unknown_prob )
                
                new_seq = copy.deepcopy( seq )
                new_seq[node.id] = nt
                
                self.explode( config, n+1, new_seq, prob + prob_node )                     

    def prob_joint( self, pmodel, var ):
        jp = bayes.JointProb( pmodel, [var] )
        
        result = {}
        dict = {}

        for (k, v) in jp.probs.items():
            dict[k[0]] = round(v, 3)
        
        self.prob_normalize( dict )
        
        result["*"] = dict
        
        return( result )
        
    def prob_cond( self, pmodel, var_a, vars_b ):
        cp = bayes.CondProb( pmodel, [var_a], vars_b )
        
        result = {}

        for (k, v) in cp.probs.items():
            # remember: P(k0|k1)
            k0 = k[0][0]              
            k1 = "".join( k[1] )
            
            v = round(v, 3)
            
            if( (k0 == "." and v >= 0.01) or (v >= 0.005 ) ):
                if( not result.has_key( k1 ) ):
                    result[k1] = {}
                    
                result[k1][k0] = v
        
        # normalize individually each conditional distribution 
        for (k, v) in result.items():
            self.prob_normalize( v )
        
        return( result )
    
    def prob_normalize(self, dict):
        residue = 1.0
        
        for (k, v) in dict.items():
            residue -= v
        
        for (k, v) in dict.items():
            dict[k] += residue / float(len(dict.keys()))

    def __init_chains(self):
        self.chains_length = {}
        
        for node in self.nodes:
            self.chains_length[node.chain] = self.chains_length.get( node.chain, 0 ) + 1
                    
        self.chains_count = len(self.chains_length)

    def __init_offsets(self):
        # build the first offset
        self.OFFSETS = [[range(self.chains_length[i]) for i in xrange(self.chains_count)]]
        
        # go through each node positions
        for node in self.nodes:
            # if the node position can contain a gap
            if( node.has_gap() ):
                new_offsets = []
                
                for offset in self.OFFSETS:
                    # shift all positions to the right of the current node
                    new_offset = copy.deepcopy( offset )
                    new_offsets.append( new_offset )
                    
                    for i in xrange( len(new_offset[node.chain])-1, node.ndx, -1 ):
                        new_offset[node.chain][i] = new_offset[node.chain][i-1]
                    
                    new_offset[node.chain][node.ndx] = None
                
                # extend the original list
                self.OFFSETS.extend( new_offsets )
        
        self.count = len(self.OFFSETS)
        
    # prepare a list of all mandatory pairing positions in all gap combinations allowed by the model
    def __init_pairing_list(self):
        self.pairs_list = set()
        
        for offset in self.OFFSETS:
            pairs = []
            
            for pair in self.pairing:
                n1 = self.nodes[pair.id1]
                n2 = self.nodes[pair.id2]
                
                p = (n1.chain, offset[n1.chain][n1.ndx], n2.chain, offset[n2.chain][n2.ndx], pair.ptype)
                
                if( p[1] is None or p[3] is None ):
                    print "\nERROR: Canonical pairs should not allow gaps."
                    print "       Check the pairing (%d, %d)" %(pair.id1, pair.id2)
                    print "       Check the probabilities of pairing nodes"
                    quit()

                pairs.append( p )
            
            self.pairs_list.add( tuple(pairs) )
        
        # reconvert to list
        self.pairs_list = map( lambda x: list(x), list(self.pairs_list) )
            

    def __init_order(self):
        graph = []

        # a) prepare the initial graph
        # b) count the number of chains
        chains = set()
        for node in self.nodes:
            chains.add( node.chain )
            #             node pointer, preferential order, parents 
            graph.append( [node, []] )

        # order the graph by probabilities
        graph.sort( cmp=cmp_order )
        
        # prepare all the combinations of directions
        gap_dirs = [[]]
        for chain in chains:
            new_gap_dirs = []
            for gap_dir in gap_dirs:
                new_gap_dirs.append( gap_dir + [-1] )
                new_gap_dirs.append( gap_dir + [1] )
            gap_dirs = new_gap_dirs

        # compute the best order
        order_best = None
        div_best = sys.maxint
        circles = ""
        
        for gap_dir in gap_dirs:
            # prepare the graph
            for g in graph:
                if( g[0].conds is None ):
                    g[1] = []
                else:
                    g[1] = copy.deepcopy( g[0].conds )
                
            for node in filter( lambda n: n.has_gap(), self.nodes ):
                dir = gap_dir[node.chain]
                
                for g in graph:
                    if( g[0].chain == node.chain ):
                        if( ((dir > 0) and (g[0].id > node.id)) or ((dir < 0) and (g[0].id < node.id)) ):
                            g[1].append( node.id )
            
            #print gap_dir
            #for g in graph:
            #    print g[0].id, g[0].has_gap(), g[1]
            
            # get circular references
            circles += str(gap_dir) + "\n"
            circles += self.__check_circular(graph)
             
            order = []
            ids = []

            go = True
            while go:
                go = False
                for (i, g) in enumerate(graph):
                    if( not i in order ):
                        insert = True
                        for id in g[1]:
                            if( not id in ids ):
                                insert = False
                                break
                        
                        if( insert ):
                            go = True
                            order.append( i )
                            ids.append( g[0].id )
                            break
                            
            if( len(order) == len(graph) ):
                # evaluate the order
                div = 0
                for (i, n) in enumerate(order):
                    div += abs(i - n)
                
                if( div_best > div ):
                    div_best = div
                    order_best = copy.deepcopy( order )
        
        if( not order_best is None ):
            # append the current node
            self.order = []
            for n in order_best:
                self.order.append( graph[n][0] )
        else:
            sys.stderr.write( "Model.__init_order() > ERROR > Bad node dependencies!" )
            sys.stderr.write( "Circular dependencies detected:\n%s" %circles )
            quit()
    
    def __check_circular(self, graph):
        result = ""
        for g in graph:
            result += self.__check_circular_aux( graph, g, [])
        
        return result
    
    def __check_circular_aux(self, graph, g, circle):
        result = ""

        new_circle = copy.deepcopy(circle)
        new_circle.append( g[0].id )
        
        if( g[0].id in circle ):
            result = ", ".join( map(str, circle) ) + ", " + str( g[0].id ) + "\n"
        else:
            for id in g[1]:
                for gaux in graph:
                    if( gaux[0].id == id ):
                        result += self.__check_circular_aux( graph, gaux, new_circle )
                        break
        
        return( result )
        
        
    def __init_order_original(self):
        self.order = []
        
        # build the first order list
        order = []
        for node in self.nodes:
            order.append( [node.id, node.max_prob(), node.conds, node.has_gap(), node.chain] )
        
        # sort by probability (inverted)
        order.sort( cmp=cmp_order )
        
        # adjust by parent/child info
        for i in xrange(len(order)):
            # if the node was already inserted, thus is None in the original list, skip it
            if( order[i] is None ):
                continue

            # add all the children that occur before the current node
            for j in xrange(i+1, len(order)):
                if( (not order[j] is None) and (not order[j][2] is None) and (order[i][0] in order[j][2]) ):
                    self.order.insert( 0, self.nodes[order[j][0]] )
                    order[j] = None

            # append the current node
            self.order.insert( 0, self.nodes[order[i][0]] )
            order[i] = None

    def __init_search_nodes(self):
        self.root = NodeSearch( self, self.OFFSETS )

class ModelDefinitionDataSource:
    def __init__(self, align, patterns=None, ref_seq=None):
        self.align = align
        self.patterns = patterns
        self.ref_seq = ref_seq
             
def cmp_order( o1, o2 ):
    p1 = o1[0].max_prob()
    p2 = o2[0].max_prob()
    
    if( p1 > p2 ):
        return -1
    elif( p1 < p2 ):
        return 1
    else:
        return 0
    
def cmp_order_original( o1, o2 ):
    if( o1[1] > o2[1] ):
        return 1
    elif( o1[1] < o2[1] ):
        return -1
    elif( o1[0] > o2[0] ):
        return 1
    elif( o1[0] < o2[0] ):
        return -1
    else:
        return 0