#!/usr/bin/env python

from rpy2.robjects import r

import math
import optparse
import random
import sys
import time

import lib.config as config
import lib.scanner as scanner
import lib.gcx as GC

def draw_seq( length ):
    return( "".join( ["ACGU"[random.randint( 0, 3 )] for i in range(length)] ) )

def next_seq( length, prev_seq ):
    if( prev_seq == "" ):
        next_list = ["A"] * length
    else:
        next_list = list(prev_seq)
    
        for (i, c) in enumerate(next_list):
            pos = ("ACGU".find( c ) + 1) % 4
            next_list[i] = "ACGU"[pos]
            
            if( pos > 0 ):
                break
            elif( i == (length-1) ):
                return( None )
    
    return( "".join( next_list ) )

def get_gc_weights( gc_classes, seq, L ):
    result = []
    
    K = L * math.log( 4.0 )

    for gc in gc_classes:
        p = 0.0
        for c in seq:
            p += math.log(gc[c])
        
        result.append( math.e ** (p + K) )
    
    return( result )

def get_gc_scores( gc_classes, seqs, prob_model ):
    result = []
    
    lseqs = map(lambda s: "".join(s), seqs)
    sseqs = "".join( lseqs ).replace( ".", "" )
    
    for gc in gc_classes:
        p = 0.0
        for c in sseqs:
            p += math.log(gc[c])
        
        result.append( max((prob_model - p) / math.log(2.0), 0.0) )
    
    return( result )

def get_gc_evalues( gc_hists, evalue_min ):
    # get all scores from all gc classes and sort them
    score_max = -100.0
    score_min = 100.0
    
    # TODO: should be a parameter
    score_step = 0.1
    
    for hist in gc_hists:
        if( len(hist.keys()) > 0 ):
            score_max = max(max(hist.keys()), score_max)
            score_min = min(min(hist.keys()), score_min)
    
    score = score_min
    scores = []
    while( score <= score_max ):
        score_key = round(score, 1)
        scores.append( score_key )
        score += score_step

    # integrate all probabilities with score > than the current
    gc_count = []
    for hist in gc_hists:
        gc_count.append( float(sum(hist.values())) )
    
    gc_evalue = []
    for (i, hist) in enumerate( gc_hists ):
        gc_evalue.append( {} )
        
        last_score = ""
        for score in scores[::-1]:
            p = hist.get( score, 0.0 ) / gc_count[i]
            gc_evalue[i][score] = max( gc_evalue[i].get(last_score, 0.0) + p, evalue_min )
            last_score = score

    return (gc_evalue, scores)

def save_data( ofile, gc_classes, gc_value, scores ):
    fo = open( ofile, "w" )
    
    # build the gc_classes info and header
    header = "#SCORES"
    fo.write( "[GC_CLASSES]\n" )
    
    for (i, gc) in enumerate(gc_classes):
        gc_class_name = "GC_%04d" %i 
        
        fo.write( gc_class_name )
        header += "\t%10s" %gc_class_name
        
        for (k, v) in gc.items():
            fo.write( "\t%s\t%.4f" %(k, v) )
        fo.write( "\n")

    # go through all scores to save E-values
    fo.write( "[E_VALUES]\n%s\n" %header )

    # write one line for each score
    for score in scores:
        fo.write( "%s" %(score) )
            
        for (i, gc) in enumerate( gc_classes ):
            fo.write( "\t%10.3E" %(gc_evalue[i].get( score, 0.0 ) ))
            
            #fo.write( "; %.3e\t%.3e\t%.3e" %(value, value/gc_count[i], gc_evalue[i].get( score, 0.0 ) ) )
        fo.write( "\n")
    
    fo.close()


def diff_gc_evalues( gce1, gce2, scores ):
    gc = len(gce1) // 2
    diff = 0.0
    Atotal = 0.0
    
    for i in range(1, len(scores)):
        score_a = scores[i-1]
        score_b = scores[i]
        
        X = score_b - score_a
        Ya1 = gce1[gc].get(score_a, None)
        Yb1 = gce1[gc].get(score_b, None)

        Ya2 = gce2[gc].get(score_a, None)
        Yb2 = gce2[gc].get(score_b, None)

        if( None not in [Ya1, Yb1, Ya2, Yb2] ):
            A1 = X * (Yb1 + abs(Yb1-Ya1)/2.0)
            A2 = X * (Yb2 + abs(Yb2-Ya2)/2.0)
                      
            diff += abs(A1-A2)
            Atotal += A1

    return diff/Atotal            

#
# MAIN
#
# Change here to configure the internal parameters of the simulation
N_SAVE = 10000
SPACER = "CCCCCC"

#
# parameters read
#
usage = "usage: %prog [options] <model> <evalues file>"
version = "%prog 1.0"

parser = optparse.OptionParser( usage=usage, version=version )

parser.add_option( "--data-path", action="store", dest="data_path", default=None,    type="string", help="Directory containing all models." )
parser.add_option( "--samples",   action="store", dest="samples",   default=1000000, type="int",    help="Number of samples." )

(options, args) = parser.parse_args()

if (len(args) != 2):
    parser.error("incorrect number of arguments")
else:
    MODEL_IN = args[0]
    FILE_OUT = args[1]
    
#
# initialization
#
gc_dict = {}
evalues = {}

for p in range( 1, 10 ):
    gc_dict[p] = GC.GC.specific( float(p) / 10.0 )
    evalues[p] = {}

#
# Change here to configure the parameters of the model
#
config = config.Config()

config.data_path = options.data_path
config.verbose = False
#config.cutoff_correction = math.log(1E-20)

config.configure()

t0 = time.time()

found_model = False
for model in config.models:
    # we just care with the specified model
    if( model.full_name().upper() == MODEL_IN ):
        found_model = True
        
        print ("computing e-values table for '%s' model" %model.name)
        print ("saving results on '%s'" %FILE_OUT)
        
        L = len(model.nodes)
        
        # compute the number of simulations
        sample_space = 4 ** L
        sample_count = options.samples
        sample_diff_count = 0
        
        gc_classes = GC.GC.get_classes(resolution=5, min_p=15, max_p=85)
        
        gc_hists = [{} for i in range(len(gc_classes))]
        gc_evalue_last = None
        
        pointers = [0, model.chains_length[0] + len(SPACER)]
        
        print ("total sampling space: %d" %sample_space)
        print ("data points to sample: %d" %sample_count)
        
        s = ""
        for count in range(min(sample_count, sample_space) ):
            # for safety, save the entire data every 'n' iteractions
            if( (count % N_SAVE == 0) and (count > 0) ):
                # compute the evalues of the sample so far
                (gc_evalue, scores) = get_gc_evalues( gc_hists, 1.0 / float(sample_space) )
                
                # save data for safety
                save_data( FILE_OUT, gc_classes, gc_evalue, scores )
                
                diff = 0.0
                if( gc_evalue_last is not None ):
                    diff = diff_gc_evalues( gc_evalue_last, gc_evalue, scores )
                    
                    if( diff < 0.001 ):
                        # IF the difference between this sample and the last one is less than 1/1000
                        sample_diff_count += N_SAVE
                        # stop when, after adding more than 100000 data points the difference is less than 1/1000
                        if( sample_diff_count > 100000 ):
                            print ("\nConverged after %d simulations" %count)
                            break
                    else:
                        # ELSE reset the counter
                        sample_diff_count = 0

                telapsed = time.time() - t0
                sys.stderr.write( "%sdata saved after %d simulations, TTF: %ds, diff: %5.4f (%d)" %('\b'*100, count, (telapsed * sample_count / count)-telapsed, diff, sample_diff_count) )
                sys.stderr.flush()
                
                gc_evalue_last = gc_evalue
        
            # STEP 1: DRAW A RANDOM SEQUENCE
            # ALTERNATIVE 1 RANDOM SAMPLING
            s = draw_seq( L )
            
            # ALTERNATIVE 2 SEQUENTIAL
            #s = next_seq( L, s )

            #if( s is None ):
            #    break
            
            # build the connected sequence to pass to the model
            sequence = s[:model.chains_length[0]] + SPACER + s[model.chains_length[0]:]
            
            # STEP 2: EVALUATE SEQUENCE PROBABILITY ACCORDING TO THE MODEL
            if( model.root.eval( sequence, pointers, GC.GC.neutral(), config.unknown_prob, config.cutoff_correction, False ) ):
                # it should always enter here
                (seqs, poss, prob_model, prob_rand ) = model.root.get_solution()
                
                # STEP 3: EVALUATE SEQUENCE SCORE ACCORDING TO EACH GC_CLASS
                gc_scores = get_gc_scores( gc_classes, seqs, prob_model )
            else:
                gc_scores = [0.0] * len(gc_classes) 
            
            # STEP 4: EVALUATE SEQUENCE WEIGHT ACCORDING TO EACH GC_CLASS
            weights = get_gc_weights( gc_classes, s, L )
                   
            # STEP 5: STORE VALUES
            for (i, (weight, gc_score)) in enumerate( zip( weights, gc_scores ) ):
                gc_score_key = round(gc_score, 1)
                
                gc_hists[i][gc_score_key] = gc_hists[i].get( gc_score_key, 0.0 ) + weight
        
        (gc_evalue, scores) = get_gc_evalues( gc_hists, 1.0 / float(sample_space) )
        save_data( FILE_OUT, gc_classes, gc_evalue, scores )
    
if( not found_model ):
    print ("WARNING: model '%s' not found" %MODEL_IN)
