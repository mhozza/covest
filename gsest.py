#! /usr/bin/env python3
import argparse
from covest import count_reads


def main(args):
    rc = count_reads(args.reads)
    print('Genome size:', int(round(rc / args.coverage)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('reads', help='Input histogram')
    parser.add_argument('coverage', type=float, help='Coverage')

    args = parser.parse_args()
    main(args)
