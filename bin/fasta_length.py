#!/usr/bin/env python
from pyfasta import Fasta


fasta_file = Fasta('contigs.fa')
print(sum(len(seq) for seq in fasta_file.values()))
