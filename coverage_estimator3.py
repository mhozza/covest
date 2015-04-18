#! /usr/bin/env python
from math import exp, log, fsum
from scipy.misc import comb, factorial
from scipy.optimize import minimize
import numpy
import argparse
from inverse import inverse

DEFAULT_K = 20
DEFAULT_READ_LENGTH = 100
ERROR_RATE = 0.03

numpy.seterr('raise')


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def load_dist(fname):
    unique_kmers = 0
    all_kmers = 0.0
    observed_ones = 0
    hist = list()

    with open(fname, 'r') as f:
        for line in f:
            i, cnt, sofar, percent = [fn(x) for fn, x in zip([int, int, int, float], line.split())]
            hist.append(cnt)
            if i == 1:
                observed_ones = cnt
            if i >= 2:
                unique_kmers += cnt
                all_kmers += i * cnt
    return all_kmers, unique_kmers, observed_ones, hist


def compute_coverage_apx(all_kmers, unique_kmers, observed_ones, k, r, error_rate=None):
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
    estimated_zeros = unique_kmers * exp(-cov)
    error_ones = max(0.0, observed_ones - estimated_ones)
    alpha = error_ones / (total_unique_kmers + estimated_zeros)
    print('Alpha:', alpha)
    # estimate probability of correct kmer
    estimated_p = estimate_p(cov, alpha)

    # function for conversion between kmer and base coverage
    kmer_to_read_coverage = lambda c: c * r / (r - k + 1)

    e = 1 - estimated_p ** (1.0 / k)

    if error_rate is not None:
        kmer_correct_prob = (1 - error_rate) ** k
        print('Correct kmer prob (est, given):', estimated_p, kmer_correct_prob)
        print('Error rate (est, given):', error_rate, error_rate)
        print('Coverage estimated from given p:', cov / kmer_correct_prob,
              kmer_to_read_coverage(cov / kmer_correct_prob))
    else:
        print('Correct kmer prob:', estimated_p)
        print('Error rate:', error_rate)

    print('Coverage estimated from estimated p:', cov / estimated_p,
          kmer_to_read_coverage(cov / estimated_p))
    return cov / estimated_p, e


def tr_poisson(l, j):
    if exp(l) == 1.0:  # precision fix
        return 1.0
    if factorial(j) == 'inf':
        return 0.0
    return min(1.0, (l ** j) / (factorial(j) * (exp(l) - 1.0)))


def compute_loglikelihood(hist, r, k, c, err):
    if err < 0 or err > 1:
        return float('-inf')
    # lambda for kmers with s errors
    l_s = lambda s: c * (3 ** -s) * (1.0 - err) ** (k - s) * err ** s
    # expected probability of kmers with s errors and coverage >= 1
    n_s = lambda s: comb(k, s) * (3 ** s) * (1.0 - exp(-l_s(s)))
    sum_n_s = fsum(n_s(t) for t in range(k + 1))
    # portion of kmers with s errors
    a_s = lambda s: n_s(s) / sum_n_s
    # probability that unique kmer has coverage j (j > 0)
    p_j = lambda j: fsum(a_s(s) * tr_poisson(l_s(s), j) for s in range(k + 1))

    return sum(hist[j] * log(p_j(j)) for j in range(1, len(hist)))


def compute_coverage(hist, r, k, guessed_c=10, guessed_e=0.05, error_rate=None):
    likelihood_f = lambda x: -compute_loglikelihood(
        hist, args.read_length, args.kmer_size, x[0], x[1]
    )
    x0 = [guessed_c, guessed_e]
    res = minimize(likelihood_f, x0, bounds=((0.0, None), (0.0, 1.0)), options={'disp': False})
    cov, e = res.x

    kmer_to_read_coverage = lambda c: c * r / (r - k + 1)

    print('Final coverage:', kmer_to_read_coverage(cov))
    print('Final error_rate:', e)
    return cov, e


def main(args):
    error_rate = args.error_rate if 'error_rate' in args else None
    all_kmers, unique_kmers, observed_ones, hist = load_dist(args.input_histogram)
    cov, e = compute_coverage_apx(
        all_kmers, unique_kmers, observed_ones,
        args.kmer_size, args.read_length, error_rate
    )
    compute_coverage(hist, args.read_length, args.kmer_size, cov, e, error_rate)


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
