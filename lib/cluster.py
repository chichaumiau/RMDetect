import math
import sys

try:
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    MPL_INSTALLED = True
except ImportError:
    MPL_INSTALLED = False

try:
    import numpy as np
    NP_INSTALLED = True
except ImportError:
    NP_INSTALLED = False

class Cluster:
    def __init__(self, model_name, cand, row, col ):
        self.model_name = model_name
        self.cands = []
        self.coords = {}
        
        self.row = row
        self.col = col
        
        self.count = 0
        self.occur = 0.0
        self.score = 0.0
        self.bpp = 0.0
        self.mi = 0.0
        
        # initializes
        self.cands = [cand]
        self.coords[(row, col)] = 1
        
    def add_cand(self, cand):
        self.cands.append( cand )

    def merge(self, other, dlimit):
        result = False
        
        if( (self.model_name == other.model_name) and (abs(self.row-other.row) <= dlimit) and (abs(self.col-other.col) <= dlimit) ):
            # merge both candidate lists
            self.cands.extend( other.cands )
            
            # merge both coords lists
            for (k, v) in other.coords.items():
                self.coords[k] = self.coords.get( k, 0 ) + v

            # gets the most representative
            max_v = max(self.coords.values())
            for ((row, col), v) in self.coords.items():
                if( v == max_v ):
                    self.row = row
                    self.col = col
            
            result = True
        
        return( result )
    
    def calc_metrics(self, seq_count, filter_code):
        # calc. 'count', 'occur', 'score', 'bpp' and 'mi'
        seqs = set()
        for cand in self.cands:
            seqs.add( cand.key )
        
        self.count = len(self.cands)
        self.occur = float(len(seqs)) / float(seq_count)
        
        score_sum = 0.0
        bpp_sum = 0.0
        seqs = []
        
        for cand in self.cands:
            score_sum += cand.score
            bpp_sum += cand.bpp
            
            seqs.append( "".join( ["".join( s ) for s in cand.chains_seqs] ) )
        
        self.score = score_sum / float(self.count)
        self.bpp = bpp_sum / float(self.count)

        # gets the model from the first candidate
        # we call it model_instance to distinguish from the model in the next line 
        model_instance = self.cands[0].model    
        self.mi = Cluster.__calc_mi( seqs, model_instance )
        self.h = Cluster.__calc_h( seqs, model_instance )
        
        (model, count, occur, score, bpp, mi, h) = (self.model_name, self.count, self.occur, self.score, self.bpp, self.mi, self.h)
        exec( filter_code )
        
        return( ok )

    def __str__(self):
        return( "model: %s count: %6d occur_(%%): %6.2f score: %7.3f bpp: %5.3f MI: %7.3f H: %7.3f cols: %4d %4d" %(self.model_name.ljust(10), self.count, self.occur*100.0, self.score, self.bpp, self.mi, self.h, self.row, self.col) )

    def __calc_mi(seqs, model):
        pair_count = len(model.pairing)

        counts = [0.0] * pair_count
        dicts_x = [{} for n in xrange(pair_count)]
        dicts_y = [{} for n in xrange(pair_count)]
        dicts_xy = [{} for n in xrange(pair_count)]
        
        for seq in seqs:
            for (i, pair) in enumerate(model.pairing):
                x = seq[pair.id1]
                y = seq[pair.id2]
                xy = x+y
                
                if( x != "." and y != "." and xy in ("AU", "CG", "GC", "GU", "UA", "UG") ):
                    counts[i] += 1.0
                    dicts_x[i][x] = dicts_x[i].get( x, 0.0 ) + 1.0
                    dicts_y[i][y] = dicts_y[i].get( y, 0.0 ) + 1.0
                    dicts_xy[i][xy] = dicts_xy[i].get( xy, 0.0 ) + 1.0

        total_mi = 0.0
        for i in xrange(pair_count):
            mi = 0.0
            
            for x in "ACGU":
                for y in "ACGU":
                    freq_xy = dicts_xy[i].get( x + y, 0.0 )
                    
                    if( (freq_xy > 0.0) and (counts[i] > 0.0) ):
                        p_xy = freq_xy / counts[i]
                        p_x = dicts_x[i].get( x, 0.0 ) / counts[i]
                        p_y = dicts_y[i].get( y, 0.0 ) / counts[i]
                        
                        mi += p_xy * math.log( p_xy / (p_x * p_y), 2 )
                        
            total_mi += mi
            
        return( total_mi / (2.0 * float(pair_count)) )

    def __calc_h(seqs, model):
        # maximum possible H for 6 diferent base pairs
        MAX_H = math.log( 6.0, 2 )
        
        pair_count = len(model.pairing)

        counts = [0.0] * pair_count
        dicts_xy = [{} for n in xrange(pair_count)]
        
        for seq in seqs:
            for (i, pair) in enumerate(model.pairing):
                xy = seq[pair.id1] + seq[pair.id2]
                
                if( xy in ("AU", "CG", "GC", "GU", "UA", "UG") ):
                    counts[i] += 1.0
                    dicts_xy[i][xy] = dicts_xy[i].get( xy, 0.0 ) + 1.0

        total_h = 0.0
        for i in xrange(pair_count):
            h = 0.0

            for xy in ("AU", "CG", "GC", "GU", "UA", "UG"):
                freq_xy = dicts_xy[i].get( xy, 0.0 )
                
                if( (freq_xy > 0.0) and (counts[i] > 0.0) ):
                    p_xy = freq_xy / counts[i]
                    h += -p_xy * math.log( p_xy, 2 )
                       
            # multiplies by the probability of observing a real base pair (reduces H for positions with high mispairing)
            h = h * (counts[i]/float(len(seqs)))
             
            total_h += h
            
        return( total_h / (MAX_H * float(pair_count)) )


    __calc_mi = staticmethod( __calc_mi )
    __calc_h = staticmethod( __calc_h )
        

class ClusterAlgorithm:
    def __init__(self, cands, seq_count, seqs, filter_code):
        self.cands = cands
        self.seq_count = seq_count
        self.seqs = seqs
        self.filter_code = filter_code
        
        self.clusters = []
        
    def go_cluster(self, dlimit):
        sys.stderr.write( "Processing %d candidates:" %len(self.cands) )
        
        aux_func = lambda l: filter( lambda x:x is not None, l )[0]

        # builds the very first cluster list
        build_index = {}
        for cand in self.cands:
            row = aux_func( cand.chains_cols[0] ) # gets the first no None entry of cols
            col = aux_func( cand.chains_cols[1] ) # gets the first no None entry of cols
            
            key = (cand.model.full_name(), row, col)

            # if there's not yet a cluster for this position
            ndx = build_index.get( key, -1 ) 
            if( ndx < 0 ):
                build_index[key] = len(self.clusters)
                self.clusters.append( Cluster( cand.model.full_name(), cand, row, col ) )
            else:
                self.clusters[ndx].add_cand( cand )
        
        # clusters the positions
        merge = True
        while( merge ):
            sys.stderr.write( "\n\tMerging %d clusters:                " %len(self.clusters) )
            
            merge = False
        
            for i in xrange(len(self.clusters)):
                if( i % 100 == 0 ):
                    sys.stderr.write( "%s%10d" %( "\b" * 10, i) )
                    sys.stderr.flush()
                
                if( not self.clusters[i] is None ):
                    for j in xrange(i+1, len(self.clusters)):

                        if( (not self.clusters[j] is None) and (self.clusters[i].merge( self.clusters[j], dlimit)) ):
                            merge = True
                            self.clusters[j] = None

            # removes all NULL clusters
            self.clusters = filter( lambda c: not c is None, self.clusters )

        # filter clusters base on occurrence
        sys.stderr.write( "\n\tFiltering %d clusters\n" %len(self.clusters) )
        
        for (i, cluster) in enumerate(self.clusters):
            ok = cluster.calc_metrics( self.seq_count, self.filter_code )
            
            if( not ok ):
                self.clusters[i] = None

        # removes all NULL clusters
        self.clusters = filter( lambda c: not c is None, self.clusters )
        
        sys.stderr.write( "\tRetained %d clusters\n" %len(self.clusters) )

    def write(self, fname):
        if( not fname is None ):
            fo = open( out_name, (mode == "a") and "a" or "w" )
        else:
            fo = sys.stdout

        for (i, cluster) in enumerate(self.clusters):
            print("Cluster: %5d -> %s" %((i+1), str(cluster)))

    def draw_matrix(self, fig_size):
        if( not MPL_INSTALLED ):
            sys.stderr.write( "\tCan't draw: Matplot lib not installed\n" )
            return
        
        if( not NP_INSTALLED ):
            sys.stderr.write( "\tCan't draw: Numpy not installed\n" )
            return
        
        # gets the maximum dimensions of the sequences used
        sys.stderr.write( "Building the matrix to draw\n" )
        
        max_cols = 0
        for seq in self.seqs.sequence_list:
            max_cols = max( max_cols, len(seq.seq_gapped) )
        
        # compute the scaling factor    
        scale = int(round(float(max_cols) / float(fig_size)))
        
        max_cols = (max_cols / scale) + 1
        
        plot_data = np.array([[0.0] * max_cols for r in xrange( max_cols )])
        plot_xx = np.array([[int(r * scale) for c in xrange( max_cols )] for r in xrange( max_cols )])
        plot_yy = np.array([[int(c * scale) for c in xrange( max_cols )] for r in xrange( max_cols )])
        
        for cluster in self.clusters:
            x = int(cluster.row / scale)
            y = int(cluster.col / scale)
            plot_data[x][y] = max(cluster.occur, plot_data[x][y])

        
        # 3D
        fig = plt.figure()
        
        ax = axes3d.Axes3D(fig)
        ax.plot_wireframe(plot_xx, plot_yy, plot_data)

        ax.set_xlabel( "cols (1st strand)" )
        ax.set_ylabel( "cols (2nd strand)" )
        ax.set_zlabel( "% of occurrence" )
        
        ax.set_zlim3d( (0.0, 1.0) )

        plt.show()
        
