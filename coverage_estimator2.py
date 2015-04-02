#! /usr/bin/env python
from math import exp
import argparse

# ABUNDANCE_DIST_FILE = 'data/experiment1.read_cov.dist'
DEFAULT_K = 20
DEFAULT_READ_LENGTH = 100
ERROR_RATE = 0.03


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def load_dist(fname):
    # counts = []
    unique_kmers = 0
    all_kmers = 0.0
    # last_nonzero = 0
    observed_ones = 1

    with open(fname, 'r') as f:
        for line in f:
            i, cnt, sofar, percent = [fn(x) for fn, x in zip([int, int, int, float], line.split())]
            # counts.append(cnt)
            # if cnt > 0:
            #     last_nonzero = i

            if i == 1:
                observed_ones = cnt

            if i >= 2:
                unique_kmers += cnt
                all_kmers += i * cnt

    # counts = counts[:last_nonzero + 1]

    # print(len(counts))
    # print(counts)
    # print(unique_kmers)
    return all_kmers, unique_kmers, observed_ones


def compute_coverage(all_kmers, unique_kmers, observed_ones, k, r, error_rate=None):
    def kmer_to_read_coverage(c):
        return c * r / (r - k + 1)

    if error_rate:
        kmer_correct_prob = (1 - error_rate) ** k

    cov = all_kmers / unique_kmers
    print('Coverage:', cov)

    # one iteration
    # for i in range(1):
    estimated_ones = unique_kmers * cov * exp(-cov)
    estimated_zeros = unique_kmers * exp(-cov)
    # print estimated_ones, estimated_zeros, estimated_zeros + estimated_ones
    cov = (all_kmers + estimated_ones) / (unique_kmers + estimated_ones + estimated_zeros)

    print('Coverage:', cov)

    # old_unique_kmers = unique_kmers  # from hist >2
    # unique_kmers = unique_kmers / (1 - exp(-cov) - cov * exp(-cov))
    # # print(old_unique_kmers, unique_kmers, unique_kmers - old_unique_kmers)
    # all_kmers += unique_kmers - old_unique_kmers  # blbost, bo zaratava nie len 1tkove ale aj 0ve a tie nechceme
    # cov = all_kmers / unique_kmers

    print('Coverage:', cov)

    # unique_kmers - su od 2ky vyssie
    estimated_ones = unique_kmers * cov * exp(-cov)
    error_ones = observed_ones - estimated_ones

    # ???!!! este neviem ci tu nema byt skor all_kmers, ale asi nie
    alpha = error_ones / (unique_kmers + observed_ones)
    print('Alpha:', alpha)
    estimated_p = estimate_p(cov, alpha)

    if error_rate:
        print('Correct kmer prob (est, given):', estimated_p, kmer_correct_prob)
        print('Error rate (est, given):', 1 - estimated_p ** (1.0 / k), error_rate)
        print('Coverage estimated from given p:', cov / kmer_correct_prob,
              kmer_to_read_coverage(cov / kmer_correct_prob))
    else:
        print('Correct kmer prob:', estimated_p, kmer_correct_prob)
        print('Error rate:', 1 - estimated_p ** (1.0 / k), error_rate)

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
