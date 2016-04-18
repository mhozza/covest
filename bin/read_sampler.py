#! /usr/bin/env python3
import random
import argparse
from os import path
from Bio import SeqIO


BASES = ['A', 'C', 'G', 'T', ]

DEFAULT_FACTOR = 2


def load_reads(fname):
    _, ext = path.splitext(fname)
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    with open(fname, "rU") as f:
        for read in SeqIO.parse(f, fmt):
            yield read.id, read.seq


def main(fname, reads_file, factor):
    prob = 1.0 / factor
    with open(fname, 'w') as f:
        for id, read in load_reads(reads_file):
            if random.random() > prob:
                f.write('>{}\n'.format(id))
                f.write('{}\n'.format(read))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Subsample reads randomly form other reads')
    parser.add_argument('-f', '--factor', type=float,
                        default=DEFAULT_FACTOR, help='Factor')
    parser.add_argument('reads', help='Reads file')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()

    main(args.output, args.reads, args.factor)
