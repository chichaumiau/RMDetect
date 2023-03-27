# changed in version 0.0.2

import math
import os
import sys

from model_parser import *
from evalues import *

# try to import psyco
# ~try:
    # ~import psyco
    # ~psyco.full()

# ~except ImportError:
    # ~sys.stderr.write( "'Psyco' package not found. Using 'Psyco' is much faster!\n" )

class Config:
    VERSION = "0.0.6"
    ENV_DATA_PATH = "RMDETECT_DATA"
    
    MODEL_SUFFIX = "model"
    EVALUES_SUFFIX = "evalues"
    
    def __init__(self):
        #
        # define all the default values
        #
        self.verbose = False
        
        self.cutoff_correction = math.log(0.1)
        self.unknown_prob = math.log(0.0001)
        self.min_bpp_mean = 0.001
        
        self.win_len = 0
        self.win_step = 0
        
        self.data_path = None

        self.both_strands = False
        
        self.models = []
        self.evalues = []
        
        self.exclude = []
        self.include = []
        
        self.constraint = ""
        
    def configure(self):
        self.__set_data_path( )
        self.__load_models()

    #
    # Define the data path (i.e., model and evalues file path)
    #
    # changed in version 0.0.2
    def __set_data_path(self):
        env_data_path = os.getenv( self.ENV_DATA_PATH )
        scr_data_path = "%s/models" %(os.path.dirname(sys.argv[0]))
        
        # the order to follow in the search is:
        
        sys.stderr.write( "Searching for models in:\n" )

        # 1 - check if the 'data_path' parameters is defined.
        sys.stderr.write( "\t1. '--data-path' option: " )
        if( (not self.data_path is None) and (os.path.isdir( self.data_path )) ):
            pass
            
        else:
            sys.stderr.write( "NOT FOUND -> '%s'\n" %self.data_path )
            sys.stderr.write( "\t2. '%s' environment var: " %self.ENV_DATA_PATH )
            
            # 2 - check if the environment variable 'ENV_DATA_PATH' is defined
            if( (not env_data_path is None) and (os.path.isdir( env_data_path )) ):
                self.data_path = env_data_path
            else:
                sys.stderr.write( "NOT FOUND -> '%s'\n" %env_data_path )
                sys.stderr.write( "\t3. Default 'models' directory: " )
            
                # 3 - checks if the software directory contains a sub directory named 'models' at the same level of the script directory
                if( os.path.isdir( scr_data_path ) ):
                    self.data_path = scr_data_path
                else:
                    sys.stderr.write( "NOT FOUND -> '%s'\n" %scr_data_path )
                    sys.stderr.write( "\t4. Current dir: " )
    
                    self.data_path = os.curdir
                
        # normalizes the path name
        self.data_path = os.path.dirname( self.data_path  + "/" )
        sys.stderr.write( "OK! -> '%s'\n\n" %self.data_path )
                   
    def __load_models(self):
        #
        # Build model list from the files in the data path 
        # 
        mparser = ModelParser()
        eparser = EValuesParser()
        if( not os.path.isdir( self.data_path ) ):
            sys.stderr.write( "ERROR: defined data_path directory '%s' doesn't exist!\n" %self.data_path )
            quit()
        flist = os.listdir( self.data_path )
        sys.stderr.write( "Loading models from: '%s' " %self.data_path )
        sys.stderr.flush()
        
        included = []
        excluded = []
        for f in flist:
            if( f.endswith( self.MODEL_SUFFIX ) ):
                sys.stderr.write( "." )
                sys.stderr.flush()
                # ~print( "===%s/%s" %(self.data_path, f) )

                model = mparser.parse( "%s/%s" %(self.data_path, f) )
                
                # models are included if:
                # a) explicitly indicated in '--inc' option (or '--inc' option is empty) 
                # b) not indicated in '--exc' option.
                #if( (model.full_name() in self.exclude) or ((len(self.include) > 0) and (not model.full_name() in self.include)) ):
                if( ((model.full_name() in self.include) or (len(self.include) == 0)) and (not model.full_name() in self.exclude) ):
                    self.models.append( model )
                    self.evalues.append( eparser.parse( "%s/%s" %(self.data_path, f.replace(self.MODEL_SUFFIX, self.EVALUES_SUFFIX) ) ) )

                    included.append( model.full_name() )
                else:
                    excluded.append( model.full_name() )
                    
        if( len(included) == 0 ):
            included = ["none"]

        if( len(excluded) == 0 ):
            excluded = ["none"]

        sys.stderr.write( "\nIncluded models:  %s\n" %", ".join( included ) )
        sys.stderr.write(   "Excluded models:  %s\n" %", ".join( excluded ) )

        if( len(self.models) == 0 ):
            print("WARNING: No models read from '%s'\n" %(self.data_path))
