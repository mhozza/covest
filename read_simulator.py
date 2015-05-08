#! /usr/bin/env python3
import random
import argparse
from os import path
from Bio import SeqIO


BASES = ['A', 'C', 'G', 'T', ]

DEFAULT_COVERAGE = 10
DEFAULT_ERROR_RATE = .03
DEFAULT_READ_LENGTH = 100


def other(base):
    b = random.randrange(len(BASES) - 1)
    if BASES[b] == base:
        return BASES[-1]
    return BASES[b]


def reverse_complement(seq):
    complement = {
        'G': 'C',
        'C': 'G',
        'A': 'T',
        'T': 'A',
    }
    return ''.join(complement[x] for x in reversed(seq))


def load_reads(fname):
    _, ext = path.splitext(fname)
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    with open(fname, "rU") as f:
        for read in SeqIO.parse(f, fmt):
            yield read.seq


def main(fname, genome_file, read_length, error_rate, coverage, error_free_fname=None):
    genome = list(load_reads(genome_file))[0]
    genome_size = len(genome)
    read_count = int(round((coverage * genome_size) / float(read_length)))
    with open(fname, 'w') as f:
        if error_free_fname:
            eff = open(error_free_fname, 'w')

        for i in range(read_count):
            pos = random.randrange(genome_size - read_length)
            original_read = genome[pos:pos + read_length]
            original_read = original_read if random.randrange(2)\
                else reverse_complement(original_read)
            read = ''.join(b if random.random() >= error_rate else other(b)
                           for b in original_read)
            f.write('>read_{}-{}\n'.format(i, pos))
            f.write('{}\n'.format(read))
            if error_free_fname:
                eff.write('>read_{}-{}\n'.format(i, pos))
                eff.write('{}\n'.format(original_read))

        if error_free_fname:
            eff.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('-c', '--coverage', type=float,
                        default=DEFAULT_COVERAGE, help='Coverage')
    parser.add_argument('-e', '--error-rate', type=float,
                        default=DEFAULT_ERROR_RATE, help='Error rate')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-f', '--error-free', help='Error free output file')
    parser.add_argument('genome', help='Genome file')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()

    main(args.output, args.genome, args.read_length, args.error_rate,
         args.coverage, args.error_free)
