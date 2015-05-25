#! /usr/bin/env python3
import argparse


def correct_c(c, r=100, k=21):
    return c * (r - k + 1) / r


def correct_ck(c, r=151, k=21):
    return c * r / (r - k + 1)


def main(args):
    if args.coverage:
        print('Genome size:', correct_ck(correct_c(args.coverage)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')

    args = parser.parse_args()
    main(args)
