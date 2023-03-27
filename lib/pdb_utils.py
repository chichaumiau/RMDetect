import copy
import os

from Bio.PDB import *

import lib.mcannotate as mcannotate

#        
# get the sequence list from a pdb either raw or indexed
#
class Residue:
    def __init__(self, chain, pos, nt, res):
        self.chain = chain
        self.pos = pos
        self.nt = nt
        self.res = res
        
    def key(self):
        return "%s:%s" %(self.chain, self.pos)
    
    def __str__(self):
        return "%s:%s:%s > %s" %(self.chain, self.pos, self.nt, self.res)
    
class Index:
    def __init__(self, ndx, chain, pos):
        self.ndx = ndx      # index in the full list of residues
        
        self.chain = chain  # chain of the residue
        self.pos = pos      # position of the residue in the chain (as in the PDB)
        
        self.id = None      # position of the residue in the full ref_seq

    def key(chain, pos):
        return( "%s:%s" %(chain, pos) )
    
    def __repr__(self):
        return( "<psb_coords: %s:%s, pos: %s, index: %s>" %(self.chain, self.pos, self.ndx, self.id) )
    
    key = staticmethod( key )
    
class Interaction:
    def __init__(self, id1, id2, wc, stack):
        self.id1 = id1
        self.id2 = id2
        self.wc = wc
        self.stack = stack

    def __repr__(self):
        return( "<Interaction: %d-%d, wc: %s, stack: %s>" %(self.id1, self.id2, self.wc, self.stack) )

class PDBStruct(object):
    def __init__(self):
        self._pdb_file = None
        self._struct = None
        self._res_list = []
        self._res_seqs = []
        self._res_index = {}
        self._interactions = []
        #self._brackets = []
        #self._wcpairs = []
    
    def load(self, pdb_file, coords ):
        self._pdb_file = pdb_file
        
        self._load_struct()
        self._load_coords( coords )
        self._load_annotations_3D()
                
    def get_sequences(self):
        sequences = []
        
        for res_seq in self._res_seqs:
            seq = ""
            for i in res_seq:
                seq += self._res_list[i].nt
            
            sequences.append( seq )
    
        return sequences

    def get_interactions(self, type="ALL"):
        return self._interactions
       
    def _load_struct(self):
        parser = PDBParser()
        self._struct = parser.get_structure( "struct", self._pdb_file )
        
        if( len(self._struct) > 1 ):
            msgs.show( "WARNING", "%d models found. Only the first will be used!" %(len(self._struct)) )

        self._res_list = []
        
        # gets only the first model
        model = self._struct[0]
        self._res_index = {}
        ndx = 0
        
        for chain in model.child_list:
            for res in chain.child_list:
                self._res_list.append( Residue(chain.id, res.id[1], res.resname.strip(), res) )
                self._res_index[Index.key(chain.id, res.id[1])] = Index( ndx, chain.id, res.id[1] )
                ndx += 1

        return( True )
            
    def _load_coords(self, coords):
        self._res_seqs = []

        entries = map( lambda row: row.split( ":" ), coords.split( "," ) )
        id = 0
        
        for entry in entries:
            if( len(entry) != 3 ):
                msgs.show( "ERROR", "Bad coord entry: '%s'" %entry)
                return( False )
            
            chain = entry[0]
            pos = int(entry[1])
            count = int(entry[2])
            
            res_seq = []

            # get the index position
            index = self._get_index( chain, pos )
            
            # get the positions
            for i in xrange( index.ndx, index.ndx + count ):
                if( i >= len(self._res_list) ):
                    print("ERROR", "Bad count %d in coords entry: '%s'" %(count, entry))
                    quit()
                
                if( self._res_list[i].chain != chain ):
                    print("ERROR", "Position %d in coords entry: '%s' is outside the chain" %(i, entry))
                    quit()
                
                res_seq.append( i )
                self._res_index[Index.key(self._res_list[i].chain, self._res_list[i].pos)].id = id
                id += 1

            self._res_seqs.append( res_seq )

    def _load_annotations_3D(self):
        self._interactions = []
        
        mca = mcannotate.MCAnnotate()
        mca.load( self._pdb_file, os.path.dirname( self._pdb_file ) )
        
        for (type, chain_a, pos_a, nt_a, chain_b, pos_b, nt_b, extra1, extra2, extra3) in mca.interactions:
            # get the rank of the first position of the pair
            index_a = self._get_index( chain_a, pos_a )
            index_b = self._get_index( chain_b, pos_b )
            
            wc = ((extra2+extra1) == "cisWW")
            stack =  (type=="STACK")
            
            if( (index_a.id is not None) and (index_b.id is not None) ):
                # discard consecutive stacking
                if( (not stack) or (abs(index_a.id-index_b.id)>1)):
                    self._interactions.append( Interaction( min( index_a.id, index_b.id ), max( index_a.id, index_b.id ), wc, stack ) )
         
    def _get_index(self, chain, pos):
        index = self._res_index.get( Index.key( chain, pos ), None )
        
        if( index is None ):
            print("ERROR", "Bad index key: '%s'" %Index.key( chain, pos ))
            quit()
        
        return( index )
