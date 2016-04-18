#!/usr/bin/env python3
import argparse

from covest import load_hist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', help='Input histogram')

    parser.add_argument('-at', '--autotrim', type=int, nargs='?', const=0,
                        help='Trim histogram automatically with this treshold')

    args = parser.parse_args()
    print('Hist size:', len(load_hist(args.input_histogram, auto_trim=args.autotrim)))
