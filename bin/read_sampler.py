#! /usr/bin/env python3
import argparse
from covest.data import sample_reads

BASES = ['A', 'C', 'G', 'T', ]

DEFAULT_FACTOR = 2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Subsample reads randomly form other reads')
    parser.add_argument('-f', '--factor', type=float,
                        default=DEFAULT_FACTOR, help='Factor')
    parser.add_argument('reads', help='Reads file')
    parser.add_argument('output', help='Output file')

    args = parser.parse_args()

    sample_reads(args.output, args.reads, args.factor)
