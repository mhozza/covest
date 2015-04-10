#! /usr/bin/env python
from math import exp
import argparse
from inverse import inverse

DEFAULT_K = 20
DEFAULT_READ_LENGTH = 100
ERROR_RATE = 0.03


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def load_dist(fname):
    unique_kmers = 0
    all_kmers = 0.0
    observed_ones = 0

    with open(fname, 'r') as f:
        for line in f:
            i, cnt, sofar, percent = [fn(x) for fn, x in zip([int, int, int, float], line.split())]
            if i == 1:
                observed_ones = cnt
            if i >= 2:
                unique_kmers += cnt
                all_kmers += i * cnt
    return all_kmers, unique_kmers, observed_ones


def compute_coverage(all_kmers, unique_kmers, observed_ones, k, r, error_rate=None):
    total_unique_kmers = unique_kmers + observed_ones
    # compute coverage from hist >=2
    cov = all_kmers / unique_kmers
    print('Coverage:', cov)
    # fix coverage
    fn = lambda cov: (cov - cov * exp(-cov)) / (1 - exp(-cov) - cov * exp(-cov))
    fix_coverage = inverse(fn)
    cov = fix_coverage(cov)
    print('Fixed coverage:', cov)
    # fix unique kmers
    unique_kmers /= (1.0 - exp(-cov) - cov * exp(-cov))
    # compute alpha (error read ratio)
    estimated_ones = unique_kmers * cov * exp(-cov)
    error_ones = max(0, observed_ones - estimated_ones)
    alpha = error_ones / total_unique_kmers
    print('Alpha:', alpha)
    # estimate probability of correct kmer
    estimated_p = estimate_p(cov, alpha)

    # function for conversion between kmer and base coverage
    kmer_to_read_coverage = lambda c: c * r / (r - k + 1)

    if error_rate is not None:
        kmer_correct_prob = (1 - error_rate) ** k
        print('Correct kmer prob (est, given):', estimated_p, kmer_correct_prob)
        print('Error rate (est, given):', 1 - estimated_p ** (1.0 / k), error_rate)
        print('Coverage estimated from given p:', cov / kmer_correct_prob,
              kmer_to_read_coverage(cov / kmer_correct_prob))
    else:
        print('Correct kmer prob:', estimated_p)
        print('Error rate:', 1 - estimated_p ** (1.0 / k))

    print('Coverage estimated from estimated p:', cov / estimated_p,
          kmer_to_read_coverage(cov / estimated_p))


def main(args):
    error_rate = args.error_rate if 'error_rate' in args else None
    all_kmers, unique_kmers, observed_ones = load_dist(args.input_histogram)
    compute_coverage(all_kmers, unique_kmers, observed_ones, args.kmer_size,
                     args.read_length, error_rate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', help='Input histogram')
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-e', '--error-rate', type=float, help='Error rate')
    args = parser.parse_args()
    main(args)
