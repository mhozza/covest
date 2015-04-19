#! /usr/bin/env python
from math import exp, log, fsum
from scipy.misc import comb, factorial
from scipy.optimize import minimize
import argparse
from inverse import inverse

# from utils import print_wrap as pw


DEFAULT_K = 20
DEFAULT_READ_LENGTH = 100
ERROR_RATE = 0.03

# INF = 1e100
INF = float('inf')


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


def compute_coverage_apx(all_kmers, unique_kmers, observed_ones, k, r):
    if unique_kmers == 0:
        return 0.0, 1.0
    total_unique_kmers = unique_kmers + observed_ones
    # compute coverage from hist >=2
    cov = all_kmers / unique_kmers
    # fix coverage
    fn = lambda cov: (cov - cov * exp(-cov)) / (1 - exp(-cov) - cov * exp(-cov))
    fix_coverage = inverse(fn)
    cov = fix_coverage(cov)
    # fix unique kmers
    unique_kmers /= (1.0 - exp(-cov) - cov * exp(-cov))
    # compute alpha (error read ratio)
    estimated_ones = unique_kmers * cov * exp(-cov)
    estimated_zeros = unique_kmers * exp(-cov)
    error_ones = max(0.0, observed_ones - estimated_ones)
    alpha = error_ones / (total_unique_kmers + estimated_zeros)
    # estimate probability of correct kmer and error rate
    estimated_p = max(0.0, estimate_p(cov, alpha))
    e = 1 - estimated_p ** (1.0 / k)
    # return corrected coverage and error estimate
    if estimated_p > 0:
        return cov / estimated_p, e
    else:
        return 0, e


def tr_poisson(l, j):
    try:
        if exp(l) == 1.0:  # precision fix
            return 1.0
    except OverflowError:
        print l, j
        return 0.0

    if factorial(j) == 'inf':
        return 0.0
    return min(1.0, (l ** j) / (factorial(j) * (exp(l) - 1.0)))


def safe_log(x):
    if x == 0:
        return -INF
    return log(x)


def compute_probabilities(r, k, c, err):
    # lambda for kmers with s errors
    l_s = lambda s: c * (3 ** -s) * (1.0 - err) ** (k - s) * err ** s
    # expected probability of kmers with s errors and coverage >= 1
    n_s = lambda s: comb(k, s) * (3 ** s) * (1.0 - exp(-l_s(s)))
    sum_n_s = fsum(n_s(t) for t in range(k + 1))

    if sum_n_s == 0:  # division by zero fix
        sum_n_s = 1
    # portion of kmers with s errors
    a_s = lambda s: n_s(s) / sum_n_s
    # probability that unique kmer has coverage j (j > 0)
    p_j = lambda j: fsum(a_s(s) * tr_poisson(l_s(s), j) for s in range(k + 1))
    return p_j


def compute_loglikelihood(hist, r, k, c, err):
    if err < 0 or err >= 1 or c <= 0:
        return -INF
    p_j = compute_probabilities(r, k, c, err)
    return sum(hist[j] * safe_log(p_j(j)) for j in range(1, len(hist)))


def compute_probabilities_with_repeats(r, k, c, err, q1, q):
    p_oj_d = dict()

    def p_oj(o, j):
        if (o, j) not in p_oj_d:
            if o == 1:
                res = compute_probabilities(r, k, c, err)(j)
            else:
                res = fsum(p_oj(1, i) * p_oj(o - 1, j - i) for i in range(1, j + 1))
            p_oj_d[(o, j)] = res
        return p_oj_d[(o, j)]

    b_o = lambda o: q1 if o == 1 else (1 - q1) * q * (1 - q) ** (o - 2)
    p_j = lambda j: fsum(
        b_o(o) * p_oj(o, j) for o in range(1, j + 1)
    )
    return p_j


def compute_loglikelihood_with_repeats(hist, r, k, c, err, q1, q):
    if err < 0 or err >= 1 or c <= 0:
        return -INF
    p_j = compute_probabilities_with_repeats(r, k, c, err, q1, q)
    return sum(hist[j] * safe_log(p_j(j)) for j in range(1, len(hist)))


def compute_coverage(hist, r, k, guessed_c=10, guessed_e=0.05,
                     error_rate=None, orig_coverage=None):
    likelihood_f = lambda x: -compute_loglikelihood(
        hist, args.read_length, args.kmer_size, x[0], x[1]
    )
    x0 = [guessed_c, guessed_e]
    res = minimize(
        likelihood_f, x0,
        bounds=((0.0, None), (0.0, 1.0)),
        options={'disp': False}
    )
    cov, e = res.x

    # function for conversion between kmer and base coverage
    kmer_to_read_coverage = lambda c: c * r / (r - k + 1)

    print('Guessed coverage:', kmer_to_read_coverage(guessed_c))
    print('Final coverage:', kmer_to_read_coverage(cov))
    if error_rate is not None:
        print('Given error rate:', error_rate)
    print('Guessed error rate:', guessed_e)
    print('Final error rate:', e)

    if error_rate is not None and orig_coverage is not None:
        print('Original loglikelihood:', -likelihood_f([orig_coverage, error_rate]))
    print('Guessed loglikelihood:', -likelihood_f(x0))
    print('Estimated loglikelihood:', -likelihood_f([cov, e]))

    return cov, e


def compute_coverage_repeats(hist, r, k, guessed_c=10, guessed_e=0.05,
                             guessed_q1=1.0, guessed_q=0.0,
                             error_rate=None, orig_coverage=None):
    likelihood_f = lambda x: -compute_loglikelihood_with_repeats(
        hist, args.read_length, args.kmer_size, x[0], x[1], x[2], x[3],
    )
    x0 = [guessed_c, guessed_e, 1, 0]

    res = minimize(
        likelihood_f, x0,
        bounds=((0.0, None), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        options={'disp': False}
    )
    # res = minimize(likelihood_f, x0, options={'disp': False})
    cov, e, q1, q = res.x

    # function for conversion between kmer and base coverage
    kmer_to_read_coverage = lambda c: c * r / (r - k + 1)

    print('Guessed coverage:', kmer_to_read_coverage(guessed_c))
    print('Final coverage:', kmer_to_read_coverage(cov))
    if error_rate is not None:
        print('Given error rate:', error_rate)
    print('Guessed error rate:', guessed_e)
    print('Final error rate:', e)

    if error_rate is not None and orig_coverage is not None:
        print('Original loglikelihood:', -likelihood_f([orig_coverage, error_rate, q1, q]))
    print('Guessed loglikelihood:', -likelihood_f(x0))
    print('Estimated loglikelihood:', -likelihood_f([cov, e, q1, q]))

    print('Estimated q1 and q:', q1, q)

    return cov, e


def main(args):
    orig_error_rate = args.error_rate if 'error_rate' in args else None
    orig_coverage = args.coverage if 'coverage' in args else None
    all_kmers, unique_kmers, observed_ones, hist = load_dist(args.input_histogram)
    cov, e = compute_coverage_apx(
        all_kmers, unique_kmers, observed_ones,
        args.kmer_size, args.read_length
    )

    compute_coverage_repeats(
        hist, args.read_length, args.kmer_size, cov, e,
        error_rate=orig_error_rate, orig_coverage=orig_coverage,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', help='Input histogram')
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-e', '--error-rate', type=float, help='Error rate')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    args = parser.parse_args()
    main(args)
