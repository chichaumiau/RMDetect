
python rmdetect.py --data-path models examples/lysine/lys.seq > result.txt
python rmout.py --data-path models < result.txt
python rmout.py --data-path models --out=lys.pdf --cands=all fold < result.txt
python rmdetect.py --data-path models examples/lysine/lys.stk > mresult.txt
python rmout.py --data-path models < mresult.txt
python rmcluster.py --data-path models --min-occur=0.1 --min-score=9.0 --min-bpp=0.01 --min-mi=0.001 < mresult.txt
python rmcluster.py --data-path models --min-occur=0.1 --min-score=9.0 --min-bpp=0.01 --min-mi=0.001 --out-dir=. < mresult.txt
python rmout.py --data-path models < cluster_001_KT_1.0.res
python rmout.py --data-path models < cluster_002_GB_1.0.res

