#! /usr/bin/env python3
import argparse
from covest.covest import count_reads_size


def main(args):
    rc = count_reads_size(args.reads)
    if args.coverage:
        print('Genome size:', int(round(rc / args.coverage)))
    if args.genome_size:
        print('Coverage:', rc / args.genome_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('reads', help='Input histogram')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    parser.add_argument('-g', '--genome-size', type=float, help='Genome size')

    args = parser.parse_args()
    main(args)
