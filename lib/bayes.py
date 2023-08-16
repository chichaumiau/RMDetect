import math
import sys
import time

#from rpy import r as R

# try to import psyco
try:
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    import numpy as np
    mpl_installed = True

except ImportError:
    mpl_installed = False


def pp( x, y, x_min=0.0, y_min=0.0, signal=1 ):
    if( signal > 0 ):
        return( (x >= x_min) and (y >= y_min) )    # positive
    else:
        return( (x >= x_min) or (y >= y_min) )    # negative
    
class Analysis:
    #
    #
    #
    def sens_spec(cands, model_name, col_pairs, x_feature, y_feature, positives, tries, signal, x_min=None, y_min=None, out_file=None, no_fig=False ):
        positive_real_keys = set()
        
        ltp  = [[], []]
        lfp  = [[], []]
        lfps = [[], []]
        ltn  = [[], []]
        ltns = [[], []]
        lfn  = [[], []]

        xlim = None
        ylim = None
        
        for cand in cands:
            # check if the candidate belongs to the model of interest
            if( cand.model.full_name() == model_name ):
                ignore = False
                
                # gets the (x, y) values according to an arbitrary definition  
                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                if( xlim is None ):
                    xlim = [x, x]
                    ylim = [y, y]
                else:
                    xlim = [min(xlim[0], x), max(xlim[1], x)]
                    ylim = [min(ylim[0], y), max(ylim[1], y)]

                positive_real = False
                shuffle = False

                if( not "SHUFFLE" in cand.key ):
                    for cols in col_pairs:
                        ok = True
                        for (i, col) in enumerate(cols):
                            # accept a tolerance of N nucleotides (for now N=2)
                            if ( abs(cand.chains_cols[i][0] - col) > 2 ):
                                ok = False
                                break
                        
                        if( ok ):
                            positive_real = True
                            positive_cols = "-".join( map( str, cols ) )
                            break
                else:
                    shuffle = True
                    
                # check if the candidate is inside the positive area
                positive_predicted = pp( x, y, x_min, y_min, signal )
                
                # 'pointer' will point to the list to be updated
                pointer = None
                if( shuffle ):
                    if( positive_predicted ):
                        pointer = lfps
                    else:
                        pointer = ltns
                else:
                    if( positive_real ):
                        positive_key = cand.key + ";" + positive_cols
                    
                        # check if it's a duplicate candidate and ignore it
                        if( not positive_key in positive_real_keys ):
                            positive_real_keys.add( positive_key )
                            
                            if( positive_predicted ):
                                pointer = ltp
                            else:
                                pointer = lfn
                    else:
                        if( positive_predicted ):
                            pointer = lfp
                        else:
                            pointer = ltn
                
                if( not pointer is None ):
                    pointer[0].append( x )
                    pointer[1].append( y )
        
        # with or without tries:
        if( tries[0] == 0 ):
            (tp, fp, tn, fn) = (len(ltp[0]), len(lfp[0]) + len(lfps[0]), len(ltn[0]) + len(ltns[0]), positives - len(ltp[0]))
        else:
            negatives = tries[0] - positives
            (tp, fp, tn, fn) = (len(ltp[0]), len(lfp[0]) + len(lfps[0]), negatives-(len(lfp[0]) + len(lfps[0])), positives - len(ltp[0]))

        (mcc, tpr, fpr, tnr, fdr) = Analysis.__calc_stats( tp, fp, tn, fn )
        
        print ("TP: %d, FP: %d, TN: %d, FN: %d" %(tp, fp, tn, fn))
        print ("MCC: %.3f" %(mcc))
        print ("TPR (sensitivity): %.3f" %tpr)
        print ("TNR (specificity): %.3f" %tnr)
        print ("FPR:               %.3f" %fpr)
        print ("FDR:               %.3f" %fdr)
        
        print ("\n#   TP\tFP\tTN\tFN\tMCC\tTPR\tTNR\tFPR\tFDR")
        print ("#OL: %d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" %(tp, fp, tn, fn, mcc, tpr, tnr, fpr, fdr))

        # plots the candidates values
        dx = (xlim[1]-xlim[0]) * 0.10
        dy = (ylim[1]-ylim[0]) * 0.10
        
        fig = plt.figure()
        
        ax = fig.add_subplot( 111 )

        ax.axhline( color="k" )
        ax.axvline( color="k" )
        ax.grid()
        ax.set_xlabel( x_feature )
        ax.set_ylabel( y_feature )
        
        ax.plot( np.array(lfp[0]),  np.array(lfp[1]),  "r1",
                 np.array(lfps[0]), np.array(lfps[1]), "ro",
                 np.array(ltn[0]),  np.array(ltn[1]),  "g1",
                 np.array(ltns[0]), np.array(ltns[1]), "go",
                 np.array(lfn[0]),  np.array(lfn[1]),  "rx",
                 np.array(ltp[0]),  np.array(ltp[1]),  "gx" )
        
        """
        ax.plot( np.array(ltp[0]),  np.array(ltp[1]), "gx", 
                 np.array(lfn[0]),  np.array(lfn[1]), "rx",
                 np.array(lfp[0]),  np.array(lfp[1]), "r1",
                 np.array(ltn[0]),  np.array(ltn[1]), "b1" )
        """
        
        ax.set_xlim( xmin=min(0-dx, xlim[0]-dx), xmax=xlim[1]+dx )
        ax.set_ylim( ymin=min(0-dy, ylim[0]-dy), ymax=ylim[1]+dy )

        if( not x_min is None ):
            ax.axvline( x=x_min, color="b" )
        if( not y_min is None ):
            ax.axhline( y=y_min, color="b" )
            
        if( not no_fig ):
            plt.show()
        
        if( not out_file is None ):
            plt.savefig()



    def data(cands, model_name, col_pairs, x_feature, y_feature, positives ):
        for cand in cands:
            # check if the candidate belongs to the model of interest
            if( cand.model.full_name() == model_name ):
                # gets the (x, y) values according to an arbitrary definition  
                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                positive = False
                if( not "SHUFFLE" in cand.key ):
                    for cols in col_pairs:
                        ok = True
                        for (i, col) in enumerate(cols):
                            # accept a tolerance of N nucleotides (for now N=2)
                            if ( abs(cand.chains_cols[i][0] - col) > 2 ):
                                ok = False
                                break
                        
                        if( ok ):
                            positive = True
                            break
                
                print ("D\t%10s\t%5s\t%7.3f\t%5.3f" %(model_name, positive, cand.bpp, cand.score))
                    
    #
    #
    #
    def statistics(cands, model_name, col_pairs, x_feature, y_feature, positives, tries, signal, stats="mcc", steps=10, out_file=None, no_fig=False ):
        xlim = None
        ylim = None
        
        clist = []
        positive_real_keys = set()
        
        scores_p = [[] for i in range(len(col_pairs))]
        scores_n = []
        bpps_p = [[] for i in range(len(col_pairs))]
        bpps_n = []
        
        print ("Positives (cols):")
        print ("\t%s" %", ".join( map( lambda cp: "(%d-%d)" %(cp[0], cp[1]), col_pairs ) ))
        print ("")

        for cand in cands:
            # check if the candidate belongs to the model of interest
            if( cand.model.full_name() == model_name ):
                ignore = False
                
                # gets the (x, y) values according to an arbitrary definition  
                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                if( xlim is None ):
                    xlim = [x, x]
                    ylim = [y, y]
                else:
                    xlim = [max(min(xlim[0], x), 0.001), max(xlim[1], x)]
                    ylim = [max(min(ylim[0], y), 0.001), max(ylim[1], y)]

                positive_real = False
                positive_cols = ""
                if( not "SHUFFLE" in cand.key ):
                    for (i, cols) in enumerate(col_pairs):
                        ok = True
                        for (j, col) in enumerate(cols):
                            # accept a tolerance of N nucleotides (for now N=2)
                            if ( abs(cand.chains_cols[j][0] - col) > 2 ):
                                ok = False
                                break
                        
                        # if it's a true candidate
                        if( ok ):
                            positive_real = True
                            positive_cols = "-".join( map( str, cols ) )
                            break

                if( positive_real ):
                    positive_key = cand.key + ";" + positive_cols
                    
                    # check if it's a duplicate candidate ignore it
                    if( positive_key in positive_real_keys ):
                        ignore = True
                    else:
                        scores_p[i].append( x )
                        bpps_p[i].append( y )
                        positive_real_keys.add( positive_key )
                else:
                    scores_n.append( x )
                    bpps_n.append( y )
                
                if( not ignore ):
                    clist.append( (positive_real, x, y) )

        if( xlim is None ):
            sys.stderr.write( "Analysis.statistics() > ERROR > No candidate found from the model '%s'" %cand.model.full_name() )
            quit()
            
        step_x = (xlim[1] - xlim[0]) / float(steps)
        step_y = (ylim[1] - ylim[0]) / float(steps)
        
        plot_data = np.array([[0.0] * steps for i in range( steps )])
        plot_xx = np.array([[0.0] * steps for i in range( steps )])
        plot_yy = np.array([[0.0] * steps for i in range( steps )])
        
        sts_max = 0.0
        sts_coords = (0.0, 0.0, 0.0, 0.0)
        stats_max = (0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)

        xx = xlim[0]
        xc = 0
        while xc < steps:
            sys.stderr.write( "%sStep: %d of %d" %("\b" * 80, xc, steps) )
            sys.stderr.flush() 
            
            yy = ylim[0]
            yc = 0
            
            while yc < steps:
                (tp, fp, tn, fn)  = (0, 0, 0, 0)
                for (positive_real, x, y) in clist:
                    # check if the candidate is inside the positive area
                    positive_predicted = pp( x, y, xx, yy, signal ) 
                    
                    if( positive_real ):
                        if( positive_predicted ):
                            tp += 1
                    else:
                        if( positive_predicted ):
                            fp += 1
                        else:
                            tn += 1
                
                
                # with or without tries:
                if( tries[0] > 0 ):
                    negatives = tries[0] - positives
                    tn = negatives - fp

                fn = positives - tp
                (mcc, tpr, fpr, tnr, fdr) = Analysis.__calc_stats( tp, fp, tn, fn )
                
                if( stats == "mcc" ):
                    sts = mcc
                elif( stats == "tpr" ):
                    sts = tpr
                elif( stats == "fpr" ):
                    sts = fpr
                    
                plot_data[xc][yc] = sts
                plot_xx[xc][yc] = xx
                plot_yy[xc][yc] = yy
                
                if( sts > sts_max ):
                    sts_max = sts
                    sts_coords = (xc, yc, xx, yy)
                    stats_max = (tp, fp, tn, fn, mcc, tpr, fpr, tnr, fdr)

                yy += step_y
                yc += 1
                
            xx += step_x
            xc += 1
        
        sys.stderr.write( "%sDONE%s\n" %("\b" * 80, " " * 80) )
        
        print ("Best %s: %.3f" %(stats, sts_max))
        print ("\t%s = %.3f" %(x_feature, sts_coords[2]))
        print ("\t%s = %.3f" %(y_feature, sts_coords[3]))
        
        (tp, fp, tn, fn, mcc, tpr, fpr, tnr, fdr) = stats_max
        print ("\nStatistics of the best MCC:" )
        print ("\tTP: %d, FP: %d, TN: %d, FN: %d" %(tp, fp, tn, fn))
        print ("\tMCC: %.3f" %(mcc))
        print ("\tTPR (sensitivity): %.3f" %tpr)
        print ("\tTNR (specificity): %.3f" %tnr)
        print ("\tFPR:               %.3f" %fpr)
        print ("\tFDR:               %.3f" %fdr)

        print ("\n#OL: %d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" %(tp, fp, tn, fn, mcc, tpr, tnr, fpr, fdr))

        print ("\n#Negative scores -> mean=%7.3f, sd=%7.3f" %(R.mean( scores_n ), R.sd( scores_n )),)
        print ("99%%=%7.3f" %R.quantile( scores_n, (0.01) )['1%'],)
        print ("95%%=%7.3f" %R.quantile( scores_n, (0.05) )['5%'],)
        print ("90%%=%7.3f" %R.quantile( scores_n, (0.10) )['10%'])
        
        for (i, cols) in enumerate(col_pairs):
            print ("#Positive scores %12s -> N=%5d, mean=%7.3f, sd=%7.3f" %(str(cols), len(scores_p[i]), R.mean( scores_p[i] ), R.sd( scores_p[i] )),)
            print ("99%%=%7.3f" %R.quantile( scores_p[i], (0.01) )['1%'],)
            print ("95%%=%7.3f" %R.quantile( scores_p[i], (0.05) )['5%'],)
            print ("90%%=%7.3f" %R.quantile( scores_p[i], (0.10) )['10%'])

        print ("\n#Negative bpps -> mean=%7.3f, sd=%7.3f" %(R.mean( bpps_n ), R.sd( bpps_n )))

        for (i, cols) in enumerate(col_pairs):
            print ("#Positive bpps %12s -> N=%5d, mean=%7.3f, sd=%7.3f" %(str(cols), len(bpps_p[i]), R.mean( bpps_p[i] ), R.sd( bpps_p[i] )),)
            print ("99%%=%7.3f" %R.quantile( bpps_p[i], (0.01) )['1%'],)
            print ("95%%=%7.3f" %R.quantile( bpps_p[i], (0.05) )['5%'],)
            print ("90%%=%7.3f" %R.quantile( bpps_p[i], (0.10) )['10%'])

        if( not no_fig ):
            R.split_screen( R.c(len(col_pairs), 1) )
            for (i, cols) in enumerate(col_pairs):
                R.screen( i+1 )
                R.hist( scores_p[i], main="", xlab="", ylab="" )
            raw_input()
        
            R.plot( R.density(scores_n), main="", xlab="", ylab="" )
            raw_input()

            R.split_screen( R.c(len(col_pairs), 1) )
            for (i, cols) in enumerate(col_pairs):
                R.screen( i+1 )
                R.hist( bpps_p[i], main="", xlab="", ylab="" )
            raw_input()
        
            R.plot( R.density(bpps_n), main="", xlab="", ylab="" )
            raw_input()

        # plots the candidates values
        fig = plt.figure()
        
        ax = axes3d.Axes3D(fig)
        
        ax.set_xlabel( x_feature )
        ax.set_ylabel( y_feature )
        ax.set_zlabel( stats )
        
        ax.plot_wireframe(plot_xx, plot_yy, plot_data, rstride=5, cstride=5)
        
        if( not out_file is None ):
            fig.savefig( out_file )
        
        if( no_fig ):
            plt.close()
        else:
            plt.show()

    def plot(cands, model_name, x_feature, y_feature, out_file=None, no_fig=False ):
        data  = [[], []]

        xlim = None
        ylim = None

        for cand in cands:
            # check if the candidate belongs to the model of interest
            if( cand.model.full_name() == model_name ):
                # gets the (x, y) values according to an arbitrary definition  
                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                if( xlim is None ):
                    xlim = [x, x]
                    ylim = [y, y]
                else:
                    xlim = [min(xlim[0], x), max(xlim[1], x)]
                    ylim = [min(ylim[0], y), max(ylim[1], y)]

                data[0].append( x )
                data[1].append( y )
        
        print ("%s: mean=%.3f; sd=%.3f" %(x_feature, R.mean(np.array(data[0])), R.sd(np.array(data[0]))))
        print ("%s: mean=%.3f; sd=%.3f" %(y_feature, R.mean(np.array(data[1])), R.sd(np.array(data[1]))))
       
        
        # plots the candidates values
        dx = (xlim[1]-xlim[0]) * 0.10
        dy = (ylim[1]-ylim[0]) * 0.10
        
        fig = plt.figure()
        
        ax = fig.add_subplot( 111 )

        # --- histogram ---
        from mpl_toolkits.axes_grid import make_axes_locatable
        
        divider = make_axes_locatable(ax)
        axHistx = divider.new_vertical(1.2, pad=0.1, sharex=ax)
        axHisty = divider.new_horizontal(1.2, pad=0.1, sharey=ax)
        fig.add_axes(axHistx)
        fig.add_axes(axHisty)
        axHistx.hist( np.array(data[0]), bins=np.arange( xlim[0], xlim[1], (xlim[1]-xlim[0])/25.0 ) )
        axHisty.hist( np.array(data[1]), bins=np.arange( ylim[0], ylim[1], (ylim[1]-ylim[0])/25.0 ), orientation='horizontal')
        
        for tl in axHistx.get_xticklabels():
            tl.set_visible(False)
        
        for tl in axHisty.get_yticklabels():
            tl.set_visible(False)
        # ---

        ax.axhline( color="k" )
        ax.axvline( color="k" )
        ax.grid()
        ax.set_xlabel( x_feature )
        ax.set_ylabel( y_feature )
        
        ax.plot( np.array(data[0]),  np.array(data[1]),  "go" )

        ax.set_xlim( xmin=min(0-dx, xlim[0]-dx), xmax=xlim[1]+dx )
        ax.set_ylim( ymin=min(0-dy, ylim[0]-dy), ymax=ylim[1]+dy )
           
        if( not no_fig ):
            plt.show()
        
        if( not out_file is None ):
            plt.savefig()



    def __calc_mi(seqs, model):
        pair_count = len(model.pairing)

        counts = [0.0] * pair_count
        dicts_x = [{} for n in range(pair_count)]
        dicts_y = [{} for n in range(pair_count)]
        dicts_xy = [{} for n in range(pair_count)]
        
        for seq in seqs:
            for (i, pair) in enumerate(model.pairing):
                x = seq[pair.id1]
                y = seq[pair.id2]
                xy = x+y
                
                if( x != "." and y != "." ):
                    counts[i] += 1.0
                    dicts_x[i][x] = dicts_x[i].get( x, 0.0 ) + 1.0
                    dicts_y[i][y] = dicts_y[i].get( y, 0.0 ) + 1.0
                    dicts_xy[i][xy] = dicts_xy[i].get( xy, 0.0 ) + 1.0

        total_mi = 0.0
        for i in range(pair_count):
            mi = 0.0
            
            for x in "ACGU":
                for y in "ACGU":
                    xy = x + y
                    
                    p_xy = dicts_xy[i].get( xy, 0.0 ) / counts[i]
                    
                    if( p_xy > 0.0 ):
                        p_x = dicts_x[i].get( x, 0.0 ) / counts[i]
                        p_y = dicts_y[i].get( y, 0.0 ) / counts[i]
                        
                        mi += p_xy * math.log( p_xy / (p_x * p_y), 2 )
                        
            total_mi += mi
            
        return( total_mi / (2.0 * float(pair_count)) )


    def __calc_h(seqs, model):
        pair_count = len(model.pairing)

        counts = [0.0] * pair_count
        dicts_xy = [{} for n in range(pair_count)]
        
        for seq in seqs:
            for (i, pair) in enumerate(model.pairing):
                xy = seq[pair.id1]+seq[pair.id2]
                
                if( xy in ("AU", "CG", "GC", "GU", "UA", "UG" ) ):
                    counts[i] += 1.0
                    dicts_xy[i][xy] = dicts_xy[i].get( xy, 0.0 ) + 1.0

        max_h = math.log( 1.0/6.0, 2.0 ) * float(pair_count)
        total_h = 0.0
        for i in range(pair_count):
            h = 0.0
            
            for (xy, freq) in dicts_xy[i].items():
                prob = freq / counts[i]
                h += prob * math.log( prob, 2 ) 
            
            total_h += h
            
        return( total_h / max_h )
            

    def matrix_old_1(cands, model_name, seq_count, x_feature="evalue", y_feature="bpp", x_min=7.5, y_min=0.2, threshold=None, out_file=None, no_fig=False):
        MARGIN = 5

        sys.stderr.write( "Processing %d candidates:\n" %len(cands) )

        min_col = sys.maxint
        max_col = 0

        mcands = {}         # dict. indexed by coordinates. In each entry stores all the cands. from that coords.
        clusters = []       # list of clusters
        
        for cand in cands:
            r = cand.chains_cols[0][0]
            c = cand.chains_cols[1][0]

            min_col = min(min_col, min(c, r))
            max_col = max(max_col, max(c+1, r+1))
            
            if( cand.model.full_name() == model_name ):
                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                # if the candidate is a "good" one
                if( True or x > x_min and y > y_min ):
                    k = (r, c)

                    # stores the candidate in the proper entry.
                    if( not mcands.has_key( k ) ):
                        mcands[k] = []
                        clusters.append( [k] )
            
                    mcands[k].append( cand )
        
        # clusters the positions
        # initially each point (r, c) belongs to its own cluster
        sys.stderr.write( "\tClustering %d points\n" %len(mcands) )
        
        for i in range(len(clusters)):
            if( clusters[i] is None ):
                continue
            
            for j in range(i+1, len(clusters)):
                if( clusters[j] is None ):
                    continue
                
                merge = False
                for (ri, ci) in clusters[i]:
                    for (rj, cj) in clusters[j]:
                        if( (abs(ri-rj) <= MARGIN) and (abs(ci-cj) <= MARGIN) ):
                            merge = True
                            break
                    if( merge ):
                        break
                
                # merge clusters
                if( merge ):
                    clusters[i].extend( clusters[j] )
                    clusters[j] = None

        # now each cluster is a list of related positions
        clusters = filter( lambda c: not c is None, clusters )

        sys.stderr.write( "\tFiltering %d clusters\n" %len(clusters) )
        
        for cluster in clusters:
            cluster_cands = {}
            
            for k in cluster:
                for cand in mcands[k]:
                    # if this sequence is not represented or if is not the 'best' one
                    if( (not cluster_cands.has_key( cand.key )) or (cluster_cands[cand.key].score < cand.score) ):
                        cluster_cands[cand.key] = cand

            # TODO: Compute the entropy here
            avg_score = 0.0
            avg_bpp = 0.0
            #seqs = []
            for (key, cand) in cluster_cands.items():
                avg_score += cand.score
                avg_bpp += cand.bpp
                #seqs.append( self.chains_seqs )
            
            perc = float(len(cluster_cands)) / float(seq_count)
            avg_score = avg_score / len(cluster_cands)
            avg_bpp = avg_bpp / len(cluster_cands)
            
            if( perc >= 0.01 ):
                print ("Cluster ->\tseqs:%6d\t%6.3f\t%%\tavg_score:\t%6.3f\tavg_bpp:\t%6.3f:" %(len(cluster_cands), perc*100.0, avg_score, avg_bpp))
                
                #for k in cluster:
                #    print k,
                
                #print 
                for cand in cluster_cands:
                    print (cand)
                    print ("\t%s" %(str(cand), ["".join( s ) for s in cand.chains_seqs]))
                    
        sys.stderr.write( "Building the matrix\n" )
        
        # Decide weather to show the graphic or not
        # Maximum matrix 500 x 500
        print (min_col, max_col, out_file, no_fig)
        
        if( (max_col < 1000) and (not ((out_file is None) and (no_fig))) ):
            plot_data = np.array([[0.0] * max_col for r in range( max_col )])
            plot_xx = np.array([[r for c in range( max_col )] for r in range( max_col )])
            plot_yy = np.array([[c for c in range( max_col )] for r in range( max_col )])

            #for cluster in clusters:
            #    total = 0
            #    for (r, c) in cluster:
            #        total += len(count[(r, c)])

            #for ((r, c), cand_keys) in count.items():
            #    perc = float(len(cand_keys)) / float(seq_count)
                
                #if( (threshold is None) or (perc > threshold) ):
            #    plot_data[r][c] = perc
            
            fig = plt.figure()
            
            ax = axes3d.Axes3D(fig)
            ax.plot_wireframe(plot_xx, plot_yy, plot_data, rstride=5, cstride=5)
    
            ax.set_xlabel( "cols (1st strand)" )
            ax.set_ylabel( "cols (2nd strand)" )
            ax.set_zlabel( "% of occurrence" )
            
            ax.set_zlim3d( (0.0, 1.0) )
    
            if( not out_file is None ):
                fig.savefig( out_file )
            
            if( no_fig ):
                plt.close()
            else:
                plt.show()

    def matrix_old_2(cands, model_name, seq_count, x_feature="score", y_feature="bpp", x_min=7.5, y_min=0.2, threshold=None, out_file=None, no_fig=False):
        MARGIN = 5
        
        sys.stderr.write( "Processing %d candidates:\n" %len(cands) )

        max_col = 0

        clusters = []       # list of clusters
        
        # builds the very first cluster list
        # all candidates with the same exact coordinates belong to the same cluster
        cluster_index = {}
        
        # only works with the candidates of the specified model
        for cand in filter( lambda c: c.model.full_name() == model_name, cands ):
            # get the coordinates of the candidate
            row = cand.chains_cols[0][0]
            col = cand.chains_cols[1][0]
            
            max_col = max( max_col, max(row+1, col+1) )

            key_coords = (row, col)
            
            # get the cluster for this coordinate or a new one
            ndx = cluster_index.get( key_coords, None )
            if( ndx is None ):
                cluster = Cluster( cand )
                cluster_index[key_coords] = len(clusters)
                clusters.append( cluster )
            else:
                cluster = clusters[ndx]

            cluster.add_candidate( cand )
        
        # merging clusters        
        merged = True
        while( merged ):
            sys.stderr.write( "\n\tMerging %d clusters (* == 100): " %len(clusters) )
            sys.stderr.flush()
            merged = False
            
            for i in range(len(clusters)):
                if( i % 100 == 0 ):
                    sys.stderr.write( "*" )
                    sys.stderr.flush()
                
                if( not clusters[i] is None ):
                    for j in range(i+1, len(clusters)):
                        if( (not clusters[j] is None) and (clusters[i].is_close( clusters[j], MARGIN )) ):
                            clusters[i].merge( clusters[j] )
                            clusters[j] = None
                            merged = True
            
            clusters = filter( lambda c: not c is None, clusters )
                
        sys.stderr.write( "\tDisplaying data\n" )
        
        # Displays data
        cluster_count = 0
        cluster_data = {}
        for cluster in clusters:
            (count, avg_score, avg_bpp, sd_score, sd_bpp, mi, rep_count, rep_row, rep_col, coords) = cluster.statistics()
            
            perc = float(count) / float(seq_count)

            for (key, count) in coords.items():
                cluster_data[key] = perc
            
            if( perc >= 0.05 ):
                cluster_count += 1
                #print "Cluster %d -> seqs: %6d %6.2f %% avg_score: %7.3f %7.3f avg_bpp: %5.3f %5.3f MI: %7.3f coords: %4d - %4d / %6d" %(cluster_count, len(cluster_cands), perc*100.0, avg_score, sd_score, avg_bpp, sd_bpp, mutual_info, max_coords[0], max_coords[1], max_count)
                print ("Cluster %d -> seqs: %6d %6.2f %% avg_score: %7.3f avg_bpp: %5.3f MI: %7.3f coords: %4d - %4d / %6d" %(cluster_count, count, perc*100.0, avg_score, avg_bpp, mi, rep_row, rep_col, rep_count))
                
                for (key, cand) in cluster.cands.items():
                    print ("\t%s" %str(cand))
                    
                print ("-----------------------------------------------")
                
        sys.stderr.write( "DONE with the clusters\n" )
            
        # Decide weather to show the graphic or not
        # Maximum matrix 500 x 500
        if( (max_col < 1000) and (not ((out_file is None) and (no_fig))) ):
            sys.stderr.write( "Building the matrix\n" )
            
            plot_data = np.array([[0.0] * max_col for r in range( max_col )])
            plot_xx = np.array([[r for c in range( max_col )] for r in range( max_col )])
            plot_yy = np.array([[c for c in range( max_col )] for r in range( max_col )])

            for ((row, col), perc) in cluster_data.items():
                plot_data[row][col] = perc
            
            # 3D
            fig = plt.figure()
            
            ax = axes3d.Axes3D(fig)
            ax.plot_wireframe(plot_xx, plot_yy, plot_data, rstride=5, cstride=5)
            #ax.plot_surface(plot_xx, plot_yy, plot_data, rstride=25, cstride=25)
    
            ax.set_xlabel( "cols (1st strand)" )
            ax.set_ylabel( "cols (2nd strand)" )
            ax.set_zlabel( "% of occurrence" )
            
            ax.set_zlim3d( (0.0, 1.0) )
    
            if( not out_file is None ):
                fig.savefig( out_file )
            
            if( no_fig ):
                plt.close()
            else:
                plt.show()

    def matrix(cands, model_name, seq_count, x_feature="evalue", y_feature="bpp", x_min=7.5, y_min=0.2, threshold=None, out_file=None, no_fig=False):
        MARGIN = 5
        
        sys.stderr.write( "Processing %d candidates:\n" %len(cands) )

        max_col = 0

        mcands = {}         # gets the keys of all sequences with candidates in a given position
        clusters = []       # list of clusters
        model = None        # pointer to the model
        
        # builds the very first cluster list
        for cand in cands:
            if( cand.model.full_name() == model_name ):
                model = cand.model
                
                row = cand.chains_cols[0][0]
                col = cand.chains_cols[1][0]

                max_col = max(max_col, max(col+1, row+1))
            
                key_coords = (row, col)

                if( not mcands.has_key( key_coords ) ):
                    mcands[key_coords] = []
                    clusters.append( [key_coords] )
        
                mcands[key_coords].append( cand )
        
        # clusters the positions
        sys.stderr.write( "\tClustering %d positions\n" %len(mcands) )
        
        for i in range(len(clusters)):
            if( not clusters[i] is None ):
                for j in range(i+1, len(clusters)):
                    if( not clusters[j] is None ):
                        drop_out = False
                        for (ri, ci) in clusters[i]:
                            for (rj, cj) in clusters[j]:
                                if( (abs(ri-rj) <= MARGIN) and (abs(ci-cj) <= MARGIN) ):
                                    # merge both clusters
                                    clusters[i].extend( clusters[j] )
                                    clusters[j] = None

                                    drop_out = True
                                    break
                            if( drop_out ):
                                break

        # now each cluster is a list of related positions
        clusters = filter( lambda c: not c is None, clusters )

        sys.stderr.write( "\tFiltering %d clusters\n" %len(clusters) )
        
        # build a list of the most representative candidates for each cluster
        cluster_data = {}
        cluster_count = 0
        
        for cluster in clusters:
            cluster_cands = {}
            
            for key_coords in cluster:
                for cand in mcands[key_coords]:
                    # if this sequence is not represented or if is best than the one represented
                    if( (not cluster_cands.has_key( cand.key )) or (cluster_cands[cand.key].score < cand.score) ):
                        cluster_cands[cand.key] = cand

            # computes cluster's statistics here 
            all_scores = []
            all_bpps = []
            seqs = []
            for (key, cand) in cluster_cands.items():
                all_scores.append( cand.score )
                all_bpps.append( cand.bpp )
                seqs.append( "".join( ["".join( s ) for s in cand.chains_seqs] ) )
            
            perc = float(len(cluster_cands)) / float(seq_count)

            (avg_score, sd_score) = (R.mean(all_scores), R.sd(all_scores))
            (avg_bpp, sd_bpp) = (R.mean(all_bpps), R.sd(all_bpps))
            
            mutual_info = Analysis.__calc_mi( seqs, model )
            
            if( perc >= 0.05 ):
                cluster_count += 1
                
                s = ""
                max_count = 0
                for (row, col) in cluster:
                    s += "\tst1:%4d\tst2:%4d\tcands:%d" %(row, col, len(mcands[(row, col)]))
                    if( len(mcands[(row, col)]) > max_count ):
                        max_coords = (row, col)
                        max_count = len(mcands[(row, col)])  

                #print "Cluster %d -> seqs: %6d %6.2f %% avg_score: %7.3f %7.3f avg_bpp: %5.3f %5.3f MI: %7.3f coords: %4d - %4d / %6d" %(cluster_count, len(cluster_cands), perc*100.0, avg_score, sd_score, avg_bpp, sd_bpp, mutual_info, max_coords[0], max_coords[1], max_count)
                print ("Cluster %d -> seqs: %6d %6.2f %% avg_score: %7.3f avg_bpp: %5.3f MI: %7.3f coords: %4d - %4d / %6d" %(cluster_count, len(cluster_cands), perc*100.0, avg_score, avg_bpp, mutual_info, max_coords[0], max_coords[1], max_count))
                
                print ("")
                print (s)
                print ("")
                
                for (key, cand) in cluster_cands.items():
                    print ("\t%s" %str(cand))
                    
                print ("-----------------------------------------------")
                
            for key_coords in cluster:
                cluster_data[key_coords] = perc
        
        sys.stderr.write( "DONE with the clusters\n" )
            
        # Decide weather to show the graphic or not
        # Maximum matrix 500 x 500
        if( (max_col < 1000) and (not ((out_file is None) and (no_fig))) ):
            sys.stderr.write( "Building the matrix\n" )
            
            plot_data = np.array([[0.0] * max_col for r in range( max_col )])
            plot_xx = np.array([[r for c in range( max_col )] for r in range( max_col )])
            plot_yy = np.array([[c for c in range( max_col )] for r in range( max_col )])

            for ((row, col), perc) in cluster_data.items():
                plot_data[row][col] = perc
                print (row, col, perc)
            
            # 3D
            fig = plt.figure()
            
            ax = axes3d.Axes3D(fig)
            ax.plot_wireframe(plot_xx, plot_yy, plot_data, rstride=5, cstride=5)
            #ax.plot_surface(plot_xx, plot_yy, plot_data, rstride=25, cstride=25)
    
            ax.set_xlabel( "cols (1st strand)" )
            ax.set_ylabel( "cols (2nd strand)" )
            ax.set_zlabel( "% of occurrence" )
            
            ax.set_zlim3d( (0.0, 1.0) )
    
            if( not out_file is None ):
                fig.savefig( out_file )
            
            if( no_fig ):
                plt.close()
            else:
                plt.show()
     
    #       
    #
    #
    def matrix_old_1(cands, model_name, x_feature="evalue", y_feature="bpp", x_min=7.5, y_min=0.2, threshold=None, out_file=None, no_fig=False):
        keys = set()

        sys.stderr.write( "Processing %d candidates\n" %len(cands) )
        
        max_col = 0
        for cand in cands:
            keys.add( cand.key )
            
            row = cand.chains_cols[0][0]
            col = cand.chains_cols[1][0]
            
            max_col = max( max_col, max( row, col ) )
        
        max_col += 1
        
        sys.stderr.write( "Initalize matrix\n" )
        count = np.array([[0.0] * max_col for r in range( max_col )])
        
        for (i, cand) in enumerate(cands):
            if( cand.model.full_name() == model_name ):
                r = cand.chains_cols[0][0]
                c = cand.chains_cols[1][0]

                x = Analysis.__feature( cand, x_feature )
                y = Analysis.__feature( cand, y_feature )

                if( x > x_min and y > y_min ):
                    count[r][c] += 1.0

        plot_data = np.array([[0.0] * max_col for r in range( max_col )])
        plot_xx = np.array([[r for c in range( max_col )] for r in range( max_col )])
        plot_yy = np.array([[c for c in range( max_col )] for r in range( max_col )])
        
        hot_coords = set()

        for row in range(max_col):
            sys.stderr.write( "%sStep: %d of %d" %("\b" * 80, row, max_col) )
            sys.stderr.flush() 
            
            for col in range(max_col):

                coords_aux = set()
                for o_row in range( max(0, row-4), min(max_col, row+4) ):
                    for o_col in range( max(0, col-4), min(max_col, col+4) ):
                        plot_data[row][col] += count[o_row][o_col]
                        
                        if( (not threshold is None) and (count[o_row][o_col] > 0) ):
                            coords_aux.add( (o_row, o_col) )
                
                plot_data[row][col] = plot_data[row][col] / float(len(keys))
                
                if( (not threshold is None) and (plot_data[row][col] > threshold) ):
                    hot_coords |= coords_aux
                
        sys.stderr.write( "%sDONE%s\n" %("\b" * 80, " " * 80) )

        for (r, c) in hot_coords:
            print (r, c)
        
        fig = plt.figure()
        
        ax = axes3d.Axes3D(fig)
        ax.plot_wireframe(plot_xx, plot_yy, plot_data, rstride=5, cstride=5)

        ax.set_xlabel( "cols (1st strand)" )
        ax.set_ylabel( "cols (2nd strand)" )
        ax.set_zlabel( "% of occurrence" )
        
        ax.set_zlim3d( (0.0, 1.0) )

        if( not out_file is None ):
            fig.savefig( out_file )
        
        if( no_fig ):
            plt.close()
        else:
            plt.show()


    def __calc_stats( tp, fp, tn, fn ):
        mcc = float((tp*tn) - (fp*fn)) / np.sqrt( float((tp+fp) * (tp + fn) * (tn + fp) * (tn + fn)) )
        tpr = float(tp) / float(tp + fn)
        fpr = float(fp) / float(fp + tn)
        tnr = float(tn) / float(fp + tn)
        
        if( (fp == 0.0) and (tp == 0.0) ):
            fdr = 1.0
        else:
            fdr = float(fp) / float(fp + tp)
        
        return( (mcc, tpr, fpr, tnr, fdr) )

    def __feature(cand, choice):
        if( choice == "bpp" ):
            result = cand.bpp
        elif( choice == "evalue" ):
            if( cand.evalue == 0.0 ):
                result = 20.0
            else:
                result = -np.log(cand.evalue)
        elif( choice == "score" ):
            result = cand.score
        else:
            print ("ERROR: wrong choice")
            quit()
        
        return( result )
    
    sens_spec = staticmethod( sens_spec )
    data = staticmethod( data )
    statistics = staticmethod( statistics )
    matrix = staticmethod( matrix )
    plot = staticmethod( plot )
    
    __calc_stats = staticmethod( __calc_stats )
    __calc_mi = staticmethod( __calc_mi )
    __feature = staticmethod( __feature )
