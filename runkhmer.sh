K=20
# INFILE=data/chr14_illumina4_errFree.fa
# INFILE=data/chr14_illumina4.fq
# INFILE=data/contigs.fa
# INFILE=data/chrXV_illumina2.fq
# INFILE=data/chrXV_illumina2_errFree.fa
INFILE=data/experiment1.fa
INFILE_BASE=${INFILE%.f?}
# ./khmer/scripts/load-into-counting.py -x 1e9 -k 20 $INFILE_BASE.kh $INFILE
# ./khmer/scripts/load-into-counting.py -x 1e8 -k 20 $INFILE_BASE.kh $INFILE # smaller tables
# ./khmer/sandbox/calc-median-distribution.py $INFILE_BASE.kh $INFILE $INFILE_BASE.read_cov.dist
# ./khmer/scripts/abundance-dist.py $INFILE_BASE.kh $INFILE $INFILE_BASE.read_cov.dist
./khmer-recipes/005-estimate-total-genome-size/plot-coverage-dist.py $INFILE_BASE.read_cov.dist $INFILE_BASE.read_cov2.png --xmax=20 --ymax=4000000
./khmer-recipes/005-estimate-total-genome-size/estimate-total-genome-size.py $INFILE $INFILE_BASE.read_cov.dist 10
