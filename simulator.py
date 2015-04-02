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


def main(genome_size, read_length, error_rate, read_count):
    genome = ''.join(random.choice(BASES) for _ in range(genome_size))
    for i in range(read_count):
        pos = random.randrange(genome_size - read_length)
        read = genome[pos:pos + read_length]
        read = ''.join(b if random.random() >= error_rate else other(b) for b in read)
        print('>read_{}-'.format(i, pos))
        print(read)


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
    args = parser.parse_args()

    read_count = int(round((args.coverage * args.genome_size) / float(args.read_length)))
    main(args.genome_size, args.read_length, args.error_rate, read_count)
