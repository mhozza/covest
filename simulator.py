#! /usr/bin/env python
import random
import argparse

BASES = ['A', 'C', 'G', 'T', ]

DEFAULT_GENOME_SIZE = 10 ** 6
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


def main(fname, genome_size, read_length, error_rate, read_count, error_free_fname=None):
    genome = ''.join(random.choice(BASES) for _ in range(genome_size))
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
    parser.add_argument('-g', '--genome-size', type=int,
                        default=DEFAULT_GENOME_SIZE, help='Genome size')
    parser.add_argument('-c', '--coverage', type=float,
                        default=DEFAULT_COVERAGE, help='Coverage')
    parser.add_argument('-e', '--error-rate', type=float,
                        default=DEFAULT_ERROR_RATE, help='Error rate')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-f', '--error-free', help='Error free output file')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()

    read_count = int(round((args.coverage * args.genome_size) / float(args.read_length)))
    main(args.output, args.genome_size, args.read_length, args.error_rate,
         read_count, args.error_free)
