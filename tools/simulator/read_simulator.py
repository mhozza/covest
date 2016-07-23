#! /usr/bin/env python3
import random
import argparse
from os import path
from Bio import SeqIO

from covest.data import load_reads

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
        'N': 'N',
    }
    return ''.join(complement[x] for x in reversed(seq))


def substitute(sequence):
    subst_table = {
        'W': ('A', 'T'),
        'S': ('C', 'G'),
        'M': ('A', 'C'),
        'K': ('G', 'T'),
        'R': ('A', 'G'),
        'Y': ('C', 'T'),
        'B': ('C', 'G', 'T'),
        'D': ('A', 'G', 'T'),
        'H': ('A', 'C', 'T'),
        'V': ('A', 'C', 'G'),
        'N': ('A', 'C', 'G', 'T'),
    }

    def subst(ch):
        if ch in BASES:
            return ch
        else:
            return random.choice(subst_table[ch])

    return ''.join(
        subst(ch) for ch in sequence
    )


def main(fname, genome_file, read_length, error_rate, coverage, error_free_fname=None,
         substitude_nucleotides=False):
    total_genome_size = 0
    with open(fname, 'w') as f:
        if error_free_fname:
            eff = open(error_free_fname, 'w')
        for genome_id, genome in load_reads(genome_file):
            genome = genome.upper()
            if substitude_nucleotides:
                genome = substitute(genome)
            genome_size = len(genome)
            total_genome_size += genome_size
            read_count = int(round((coverage * genome_size) / float(read_length)))

            for i in range(read_count):
                pos = random.randrange(genome_size - read_length)
                original_read = genome[pos:pos + read_length]
                original_read = original_read if random.randrange(2)\
                    else reverse_complement(original_read)
                read = ''.join(b if random.random() >= error_rate else other(b)
                               for b in original_read)
                f.write('>read_{}_{}-{}\n'.format(genome_id, i, pos))
                f.write('{}\n'.format(read))
                if error_free_fname:
                    eff.write('>read_{}_{}-{}\n'.format(genome_id, i, pos))
                    eff.write('{}\n'.format(original_read))

        if error_free_fname:
            eff.close()

    print('Total genome size:', total_genome_size)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('-c', '--coverage', type=float,
                        default=DEFAULT_COVERAGE, help='Coverage')
    parser.add_argument('-e', '--error-rate', type=float,
                        default=DEFAULT_ERROR_RATE, help='Error rate')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-s', '--substitude', action='store_true',
                        help='Substitude extended nucleotides')
    parser.add_argument('-f', '--error-free', help='Error free output file')
    parser.add_argument('genome', help='Genome file')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()

    main(args.output, args.genome, args.read_length, args.error_rate,
         args.coverage, args.error_free, args.substitude)
