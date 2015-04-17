K=20
INFILE=$1
INFILE_BASE=${INFILE%.f?}
./khmer/scripts/load-into-counting.py -x 1e9 -k 20 $INFILE_BASE.kh $INFILE
# ./khmer/scripts/load-into-counting.py -x 1e8 -k 20 $INFILE_BASE.kh $INFILE # smaller tables
./khmer/scripts/abundance-dist.py $INFILE_BASE.kh $INFILE $INFILE_BASE.read_cov.dist
./khmer-recipes/005-estimate-total-genome-size/plot-coverage-dist.py $INFILE_BASE.read_cov.dist $INFILE_BASE.read_cov2.png --xmax=100 --ymax=100
./khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py $INFILE $INFILE_BASE.read_cov.dist 10
