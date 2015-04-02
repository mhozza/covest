#!/bin/bash
K=31
P=16

ln -s ../../Data/quakeCor/*fastq .
shuffleSequences_fastq.pl frag_1.fastq frag_2.fastq         frag.fastq
shuffleSequences_fastq.pl shortjump_1.fastq shortjump_2.fastq - | revFastq.pl > shortjump.fastq
shuffleSequences_fastq.pl longjump_1.fastq longjump_2.fastq longjump.fastq

velveth . $K -fastq -shortPaired frag.fastq -shortPaired2 shortjump.fastq  -shortPaired3 longjump.fastq
velvetg . -exp_cov auto -ins_length 180 -ins_length_sd 20 -ins_length2 3000 -ins_length2_sd 300 ins_length3 35000 -ins_length3_sd 3500 -scaffolding yes 

ln -s contigs.fa genome.scf.fasta
scf2ctg.pl < genome.scf.fasta >  genome.ctg.fasta
