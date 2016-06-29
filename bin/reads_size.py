#! /usr/bin/env python3
import argparse

from covest.data import count_reads_stats


def main(args):
    print(count_reads_stats(args.reads)[1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('reads', help='Input histogram')

    args = parser.parse_args()
    main(args)
