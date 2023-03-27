import math, sys

class Node:
    def __init__(self, id, chain, ndx, conds, probs=None):
        self.id = id
        self.chain = chain
        self.ndx = ndx
        self.conds = conds
        
        if( probs is None ):
            self.probs = {}
        else:
            self.probs = probs
        
    def max_prob(self):
        max_p = math.log( 1e-100 )
        
        for (k1, v1) in self.probs.items():
            if( v1 == "GC" ):
                max_p = max( max_p, math.log(0.25) )  # round estimate
            else:
                for (k2, v2) in v1.items():
                    max_p = max( max_p, v2 )
        
        return( max_p )
    
    def has_gap(self):
        gap = False
        
        for (k, v) in self.probs.items():
            if( (v != "GC") and ( "." in v.keys()) ):
                gap = True
                break

        return( gap )

class NodeSearch( Node ):
    def __init__(self, model, offsets, pred=None, level=-1):
        self.root = (level < 0)

        self.model = model
        self.pred = pred # predecessor node
        self.succs = []
        self.parents = []
        self.min_dists = None

        # search data
        self.NT = ""
        self.pos = None
        self.prob_accum = None
        self.prob_rand = None
        self.best_succ = None
        
        self.level = level
        
        #
        # ONLY FOR NON ROOT NODES
        #
        if( not self.root ):
            # initialize node data
            node = self.model.order[level]
            Node.__init__(self, node.id, node.chain, node.ndx, node.conds)
            
            self.probs = node.probs 
            
            # OFFSET MUST NOT BE AMBIGOUS IN THIS CHAIN:NDX POSITION!
            self.oft = offsets[0][self.chain][self.ndx]
            
            # if we have only one offset left then compute minimum distances between chains 
            if( len(offsets) == 1 ):
                self.min_dists = []
                
                for i in range(len(offsets[0])):
                    
                    max_i = max([j for j in offsets[0][i] if j!=None])
                    
                    for j in range(len(offsets[0])):
                        if( i != j ):
                            self.min_dists.append( (i, j, max_i + self.model.sep_min ) )
            
            # if this node depends on parents

            if( not self.conds is None ):

                self.parents = [None] * len(self.conds)

                self.get_parents( self.conds, self.parents )
                # ~print(self.parents)
                if( None in self.parents ):
                    sys.stderr.write( "NodeSearch.__init__() > ERROR > Unresolved parents in node '%d'\n" %(self.id) )
                    quit()

        #
        # FOR ALL NODES NOW
        #

        # prepare the child nodes
        level_new = level + 1
 
        if( level_new < len(self.model.order) ):
            node_new = self.model.order[level_new]
            
            # split the offsets in an non ambiguous way 
            dict = {}
            for (i, offset) in enumerate(offsets):
                oft = offset[node_new.chain][node_new.ndx]
                
                dict[oft] = dict.get( oft, [] )
                dict[oft].append( i )
            
            for (k, v) in dict.items():
                offsets_new = []
                
                for i in v:
                    offsets_new.append( offsets[i] )
                
                self.succs.append( NodeSearch( self.model, offsets_new, self, level_new ) )

    def ident(self):
        return( "  |" * (self.level+1) )

    def show_all(self):
        if( self.level < 0 ):
            print("ROOT")
        else:
            self.show_node()
        
        for sn in self.succs:
            sn.show_all() 

    def show_node(self):
        if( self.level < 0 ):
            print("ROOT")
        else:
            print(self.ident(), "-+ NODE -> id: %d, chain:%d, ndx:%d, oft:%s" %(self.id, self.chain, self.ndx, str(self.oft)))
            
            if( len(self.parents) > 0 ):
                print(self.ident(), "  ** PARENTS: ",)
                for p in self.parents:
                    print(p.id,)
                print
    
    def get_parents(self, conds, parents):
        # check if this node is root
        if( not self.pred is None ):
            # ~print(self.id, list(conds))
            # ~print(parents)
            if( self.id in list(conds) ):
                parents[conds.index( self.id )] = self

            self.pred.get_parents( conds, parents )
        
    def get_solution(self):
        # prepare the result
        # ~seqs = [[""] * self.model.chains_length[c] for c in xrange(self.model.chains_count)]
        # ~poss = [[0] * self.model.chains_length[c] for c in xrange(self.model.chains_count)]
        seqs = [[""] * self.model.chains_length[c] for c in range(self.model.chains_count)]
        poss = [[0] * self.model.chains_length[c] for c in range(self.model.chains_count)]

        # start iteration
        snode = self
        
        while not snode is None:
            if( not snode.root ):
                # get NT
                seqs[snode.chain][snode.ndx] = snode.NT
                poss[snode.chain][snode.ndx] = snode.pos
                
                # get pos
                poss[snode.chain][snode.ndx] = snode.pos
            
            snode = snode.succ_best
        
        return( seqs, poss, self.prob_best, self.rand_best )
    
    def eval(self, sequence, pointers, GC, unknown_prob, cutoff_correction, verbose):
        # start with generic definitions
        self.pos = None
        self.NT = ""
        self.prob_node = 0.0
        self.prob_rand = 0.0
        
        # check for chain overlapping
        if( not self.min_dists is None ):
            for (i, j, dist) in self.min_dists:
                d = pointers[j] - pointers[i]

                if( (d >= 0) and (d < dist) ):
                    if( verbose ):
                        print("Quit by overlap! expected dist (%d) >= real dist (%d)" %(dist, d))

                    return( False ) 
        
        #
        # ONLY FOR NON ROOT NODES
        #
        if( not self.root ):
            # get the nucleotide corresponding to this node
            if( self.oft is None ):
                self.pos = None
                self.NT = "."
            else:
                self.pos = pointers[self.chain] + self.oft
                self.NT = sequence[self.pos]       
            
            # compute the probability associated with this node
            if( self.conds is None ):
                prob_key = "*"
            else:
                prob_key = ""
                
                for parent in self.parents:
                    prob_key += parent.NT
            
            prob_dict = self.probs.get( prob_key, {} )
            
            if( prob_dict == "GC" ):
                prob_node = GC[self.NT]
            else:
                prob_node = prob_dict.get( self.NT, unknown_prob )                    
    
            # compute the probability of the random model
            prob_rand = GC.get( self.NT, math.log(0.25) )    # 0.25 is for the gap case
            #self.prob_rand = GC.get( self.NT, self.prob_node )      # 'prob' is for the gap case
    
            # compute the model probability
            self.prob_node = self.pred.prob_node + prob_node
            self.prob_rand = self.pred.prob_rand + prob_rand

        # final validation
        ok = False
        
        if( verbose and self.root ):
            print("START SCAN")
            print("Model: ", self.model.name, self.model.version )
        if( verbose and not self.root ):
            ident =  "\t" * self.level
            print("")
            print(ident, "ID:", self.id)
            print(ident, "Pos:", self.pos, "nt:", self.NT, "prob_key: ", prob_key, "conditional: ", self.conds)
            print(ident, "Local  (log[p])> p_node:", prob_node, ", p_rand:", prob_rand, ", score:", (prob_node - prob_rand) / math.log(2.0))
            print(ident, "Global (log[p])> p_node:", self.prob_node, ", p_rand:", (self.prob_rand + cutoff_correction), ", *SCORE*:", (self.prob_node - self.prob_rand) / math.log(2.0))
            print(ident, "continue?", self.prob_node >= (self.prob_rand + cutoff_correction), "by: ", self.prob_node - (self.prob_rand + cutoff_correction))
            if( not self.prob_node >= (self.prob_rand + cutoff_correction) ):
                print(ident, "JUMP OFF\n\n")
        
        # if we're above the random model then continue
        if( self.prob_node >= (self.prob_rand + cutoff_correction) ):
            self.succ_best = None
            self.prob_best = -500   # in the absence of a MIN_INT -500 should be a sufficiently small value
            self.rand_best = -500   # in the absence of a MIN_INT -500 should be a sufficiently small value

            if( len(self.succs) == 0 ):
                if( verbose and not self.root ):
                    print(ident, "FINISHING RECURSION!!")

                self.succ_best = None
                self.prob_best = self.prob_node
                self.rand_best = self.prob_rand
                ok = True
            else:
                if( verbose and not self.root ):
                    print(ident, "We still have %d succs" %(len(self.succs)))

                for succ in self.succs:
                    succ_ok = succ.eval( sequence, pointers, GC, unknown_prob, cutoff_correction, verbose )
                    
                    # is this node the last one?
                    if( succ_ok and (succ.prob_best > self.prob_best) ):
                        self.succ_best = succ
                        self.prob_best = succ.prob_best
                        self.rand_best = succ.rand_best
                        
                        # we just need one succ_ok!
                        ok = True

        return( ok )
