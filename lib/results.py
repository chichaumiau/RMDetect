import sys

import sequences
import candidates

class Results:
    #
    # mode:
    #    "n" - starts writing a new file
    #    "a" - appends data to an allready writen file
    #    "w" - rewrites all the data
    #
    def write(out_name, cands, seqs, mode, tries, constraint="", comments=""):
        if( not out_name is None ):
            fo = open( out_name, (mode == "a") and "a" or "w" )
        else:
            fo = sys.stdout
        
        if( comments != "" ):
            fo.write( "%s" %comments)
        
        if( (mode == "n") or (mode == "w") ):
            if( seqs.from_stdin ):
                fo.write( "##seq = %s\n" %(seqs.sequence_list[0].seq))
            else:
                fo.write( "##seq_file = %s\n" %(seqs.fname) )
                
            if( constraint != "" ):
                fo.write( "##constraint = %s\n" %(constraint) )
                
            fo.write( "##both_strands = %s\n" %(seqs.both_strands and "YES" or "NO") )
            fo.write( "##sequences = %s\n" %(seqs.seq_count()) )
            fo.write( "##tests_all = %d\n" %(tries[0]) )
            fo.write( "##tests_pairs = %d\n" %(tries[1]) )
        
        for cand in cands:
            fo.write( "%s\n" %cand )
            
        if( not out_name is None ):
            fo.close()

    def load(in_name, models):
        fname = None
        seq = None
        constraint = ""
        both_strands = False
        cands = []
        count = -1
        tries = [0,0]
        
        if( not in_name is None ):
            fi = open( in_name )
        else:
            fi = sys.stdin
        
        for line in fi:
            line = line.strip()
            
            if( line.startswith( "##seq = " ) ):
                seq = line.split( "=" )[1].strip()
            elif( line.startswith( "##seq_file = " ) ):
                fname = line.split( "=" )[1].strip()
            elif( line.startswith( "##constraint = " ) ):
                constraint = line.split( "=" )[1].strip()
            elif( line.startswith( "##both_strands = " ) ):
                both_strands = (line.split( "=" )[1].strip() == 'YES')
            elif( line.startswith( "##sequences = " ) ):
                count = int(line.split( "=" )[1].strip())
            elif( line.startswith( "##tests_all = " ) ):
                tries[0] = int(line.split( "=" )[1].strip())
            elif( line.startswith( "##tests_pairs = " ) ):
                tries[1] = int(line.split( "=" )[1].strip())
            elif( line.startswith( "#" ) or (line == "") ):
                continue
            else:
                cands.append( candidates.CandidateParser.parse( line, models ) )

        seqs = sequences.Sequences(fname=fname, seq=seq, both_strands=both_strands)

        if( not in_name is None ):
            fi.close()

        return( seqs, cands, constraint, count, tries )
    
    write = staticmethod( write )
    load = staticmethod( load )