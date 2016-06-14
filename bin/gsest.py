#! /usr/bin/env python3
import argparse


def main(args):
    rc = args.reads_size
    if args.coverage:
        print('Genome size:', int(round(rc / args.coverage)))
    if args.genome_size:
        print('Coverage:', rc / args.genome_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute genome size from reads-size and coverage')
    parser.add_argument('-rs', '--reads-size', required=True, help='Reads size')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    parser.add_argument('-g', '--genome-size', type=float, help='Genome size')

    args = parser.parse_args()
    main(args)
