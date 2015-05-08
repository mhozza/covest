#! /usr/bin/env python3
import random
import argparse

BASES = ['A', 'C', 'G', 'T', ]

DEFAULT_GENOME_SIZE = 10 ** 6


def main(fname, genome_size):
    genome = ''.join(random.choice(BASES) for _ in range(genome_size))
    with open(fname, 'w') as f:
            f.write('>simulated\n')
            f.write(genome)
            f.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('-g', '--genome-size', type=int,
                        default=DEFAULT_GENOME_SIZE, help='Genome size')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()
    main(args.output, args.genome_size)
