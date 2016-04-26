#! /usr/bin/env python
import argparse

from covest import config
from covest.data import load_histogram, process_histogram
from covest.covest import compute_coverage_apx


def relative_difference(a, b):
    res = abs(a - b) / a
    print(a, b, res)
    return res


def main(args):
    k, r = args.kmer_size, args.read_length
    hist_orig = load_histogram(args.input_histogram)
    sample_factor = 1
    hist, sample_factor = process_histogram(hist_orig, sample_factor=sample_factor)
    c_start, _ = compute_coverage_apx(hist, k, r)
    c = c_start
    target_coverage = 12
    while c > target_coverage:
        sample_factor += 1
        hist, sample_factor = process_histogram(hist_orig, sample_factor=sample_factor)
        c, _ = compute_coverage_apx(hist, k, r)
        print(c)
    print('Suggested sf:', sample_factor)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze histogram')
    parser.add_argument('input_histogram', help='Input histogram')
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=config.DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=config.DEFAULT_READ_LENGTH, help='Read length')

    args = parser.parse_args()
    main(args)
