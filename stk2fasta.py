import lib.file_io as file_io

fi = open( "/media/Elements/3-PHD/work/20091126-project-m_catalog/tests/dummy/lys.stk" )
stk = file_io.Stk( fi )
fi.close()

fo = open( "/media/Elements/3-PHD/work/20091126-project-m_catalog/tests/dummy/lys_shuffle.fasta", "w" )

for (i, (key, seq)) in enumerate(stk.seqs):
    if( not "SHUFFLE" in key):
        continue
    
    (name, ref, coords) = key.split( "/" )
    (name1, name2) = name.split( "_" )
    
    if( name2 == "sp." ):
        name = name1[:7] + "_sp"
    else:
        name = (name1[0] + "_" + name2)[:10]
        
    name = name + "_%04d" %(i)
    
    fo.write( ">%s\n%s\n" %(name, seq.replace( ".", "" ) ) )

fo.close()
    
