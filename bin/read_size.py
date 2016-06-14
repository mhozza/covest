#! /usr/bin/env python3
import argparse
from covest.covest import count_reads_size


def main(args):
    rc = count_reads_size(args.reads)
    print(rc)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('reads', help='Input histogram')

    args = parser.parse_args()
    main(args)
