import copy
import math
import random
import time

#import models
from candidates import *
#  import gcx 
from rnatools import *
from gcx import *

print(GC())

#  print(gcx.GC)

class Scanner:
    def __init__(self, config):
        self.config = config
        
        self.rna_tools = RNATools()
        
        self.cback_before_scan = lambda i: 0            # before starting the search
        self.cback_before_sequence = lambda i, j, s: 0  # before each sequence
        self.cback_before_window = lambda i, j: 0       # before each window inside the sequence
        self.cback_after_window = lambda l, t: 0        # after each window, returns all the candidates found in the sequence until now 
        self.cback_after_sequence = lambda l, t: 0      # after each sequence, returns all the candidates found in the sequence
        self.cback_after_scan = lambda l, t: 0          # after all the search, returns all the candidates
    
    #
    # 'scan' goes through all the sequences to search in
    #
    def scan(self, seqs, gc=None):
        cands = []
        tries = [0, 0]
        
        sequences = seqs.sequence_list

        # --- CALLBACK ---
        self.cback_before_scan( len(sequences) )

        for (i, sequence) in enumerate(sequences):

            # --- CALLBACK ---
            self.cback_before_sequence( len(sequences), i, sequence.key )

            if( (self.config.win_len == 0) or (self.config.win_step == 0) ):
                (cands_seq, tries_seq) = self.search( sequence.key, sequence.seq, gc )

                # update the columns of the candidates
                for cand in cands_seq:
                    cand.update_chains_cols( sequence.col_index )
            else:
                (cands_seq, tries_seq) = self.slide( sequence.key, sequence.seq, sequence.col_index, gc )

            # remove duplicate candidates
            cands_seq = self.__filter_duplicates(cands_seq)
            
            # --- CALLBACK ---
            self.cback_after_sequence( cands_seq, tries_seq )
            
            # add the candidates of the last sequence to the full list of candidates
            cands.extend( cands_seq )

            tries[0] += tries_seq[0]
            tries[1] += tries_seq[1]
            
        # --- CALLBACK ---
        self.cback_after_scan( cands, tries )

        # return the full list of candidates
        return( cands, tries )

    #
    # 'slide' goes through all the sliding windows of a single sequence
    #
    def slide(self, key, seq, col_index, gc=None):
        if( gc is None ):
            gc = GC.compute( seq )
        
        cands = []
        tries = [0, 0]
        
        ptr = 0
        ok = True
        
        while (ok):
            # --- CALLBACK ---
            self.cback_before_window( len(seq), ptr )

            end = ptr + self.config.win_len
            
            if( end >= len(seq) ):
                ptr = max( 0, len(seq) - self.config.win_len )
                end = len(seq)
                ok = False
            
            seq_win = seq[ptr:ptr + self.config.win_len]
            
            (cands_win, tries_win) = self.search( key, seq_win, gc )
                
            # update the position and columns of the candidates
            for cand in cands_win:
                cand.update_chains_pos( ptr )
                cand.update_chains_cols( col_index )

            cands.extend( cands_win )
            
            tries[0] += tries_win[0]
            tries[1] += tries_win[1]
            
            # --- CALLBACK ---
            self.cback_after_window( cands, tries )
            
            ptr += self.config.win_step
                               
        return( cands, tries )

    #
    # 'search' applies each model to the given sliding window
    #
    def search(self, key, seq, gc=None):
        tries = [0, 0]
        cands = []
        
        if( len(seq) > 0 ):
            seq_len = len(seq)
            coords = (0, seq_len)

            if( gc is None ):
                gc = GC.compute( seq )
    
            # get the sequence base pair probabilities subject to the eventual constraints
            bpprobs = self.rna_tools.bpprobs( seq, self.config.constraint )
            
            # for each model
            for (model, evalues) in zip(self.config.models, self.config.evalues):
                evalues.select_gc( gc )
                
                # starts one pointer for each chain in the model
                pointers = [0] * model.chains_count

                while True:
                    # get the sequences for each chain
                    tries[0] += 1
                    
                    if( self.__get_bpp_mean( model, bpprobs, pointers ) >= self.config.min_bpp_mean ):
                        tries[1] += 1
                        
                        verbose = False
                        #verbose = (pointers[0] == 59 and pointers[1] == 106)

                        if( model.root.eval( seq, pointers, gc, self.config.unknown_prob, self.config.cutoff, verbose ) ):
                            (seqs, poss, prob_model, prob_rand ) = model.root.get_solution()
                            
                            cand = Candidate(  model, key, coords )
                            cand.set_data_build(prob_model, prob_rand, seqs, poss, evalues )
                            cands.append( cand )
                    
                    if( not self.pointer_inc( model, seq_len, pointers ) ):
                        break

        return( cands, tries )

    def pointer_inc(self, model, seq_len, pointers):
        result = True
        
        # number of chains to search
        i = (model.chains_count-1)
        
        # increment sequentially the pointer for each chain
        while True:
            # increment the pointer for chain 'i'
            pointers[i] += 1
            
            # if it reached the end of the sequence 
            if( (pointers[i] + model.chains_length[i]) >= seq_len ):
                # increments the next strand
                i -= 1

                # if there's no next strand ends the incrementing
                if( i < 0 ):
                    result = False
                    break

                # now reset the original pointer
                if( model.symmetric ):
                    # if its symmetric do not cross the strands
                    pointers[i+1] = pointers[i] + 1
                else:
                    # otherwise put it back to 0
                    pointers[i+1] = 0
                
            else:
                break
        
        return( result )

    def __get_bpp_mean(self, model, bpprobs, pointers):
        sum_prob = 0.0
        sum_pairs = 0
        
        for pairs in model.pairs_list:
            for (c1, o1, c2, o2, ptype) in pairs:
                if( ptype == Pair.PTYPE_MUST ):
                    prob = bpprobs.prob_pair( pointers[c1]+o1, pointers[c2]+o2 )
                    
                    if( prob == 0 ): 
                        count = 0
                        break

                    sum_prob += prob
                    sum_pairs += 1

        if( sum_pairs == 0.0 ):
            return 0.0
        else:
            return (sum_prob / float(sum_pairs))

    def __filter_duplicates(self, cands):
        for i in xrange(len(cands)):
            if( cands[i] is None ):
                continue
            
            for j in xrange(i+1, len(cands)):
                if( cands[j] is None ):
                    continue

                dup = True
                (s1, e1) = cands[i].coords
                (s2, e2) = cands[j].coords
                
                if( (s1 < e2) and (s2 < e1) and (cands[i].key == cands[j].key) and (cands[i].model == cands[j].model) ):
                    for (pi, pj) in zip(cands[i].chains_pos, cands[j].chains_pos):
                        if( pi[0] != pj[0] ):
                            dup = False
                            break
                else:
                    dup = False
                
                if( dup ):
                    if( (cands[i].score > cands[j].score) and (cands[i].bpp > 0) ):
                        cands[j] = None
                    else:
                        cands[i] = None
                        break
        
        return( filter( lambda x: not x is None, cands ) )
