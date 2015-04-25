#! /usr/bin/env python
from math import exp, log, fsum
import itertools
from collections import defaultdict
from scipy.misc import comb, factorial
from scipy.optimize import minimize
import argparse
from inverse import inverse
import matplotlib.pyplot as plt
import numpy
# from utils import print_wrap as pw


DEFAULT_K = 20
DEFAULT_READ_LENGTH = 100
ERROR_RATE = 0.03

# INF = 1e100
INF = float('inf')
trim_hist = True


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def hist_trim(hist):
    z = 0
    last = -1
    trim = len(hist)
    for i, v in enumerate(hist):
        if v == 0:
            z += 1
        else:
            if last == 0:
                if i > 100 and z > len(hist) // 10:
                    trim = i - z
                    break
            z = 0
        last = v
    return hist[:trim]


def load_dist(fname):
    unique_kmers = 0
    all_kmers = 0.0
    observed_ones = 0
    hist = defaultdict(int)
    max_hist = 0

    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            cnt = int(l[1])
            hist[i] = cnt
            max_hist = max(max_hist, i)
            if i == 1:
                observed_ones = cnt
            if i >= 2:
                unique_kmers += cnt
                all_kmers += i * cnt

    hist_l = [hist[i] for i in range(max_hist)]
    if trim_hist:
        hist_l = hist_trim(hist_l)
    return all_kmers, unique_kmers, observed_ones, hist_l


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
        # function for conversion between kmer and base coverage
        kmer_to_read_coverage = lambda c: c * r / (r - k + 1)
        return kmer_to_read_coverage(cov / estimated_p), e
    else:
        return 0, e


def tr_poisson(l, j):
    with numpy.errstate(over='raise'):
        try:
            if exp(l) == 1.0:  # precision fix
                return 0.0
            return min(1.0, (l ** j) / (factorial(j) * (exp(l) - 1.0)))
        except (OverflowError, FloatingPointError):
            return 0.0


def safe_log(x):
    if x == 0:
        return -INF
    return log(x)


def compute_probabilities(r, k, c, err):
    # read to kmer coverage
    c = c * (r - k + 1) / r
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
    return sum(hist[j] * safe_log(p_j(j)) for j in range(1, len(hist)) if hist[j])


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
    return sum(hist[j] * safe_log(p_j(j)) for j in range(1, len(hist)) if hist[j])


def compute_coverage(hist, r, k, guessed_c=10, guessed_e=0.05,
                     error_rate=None, orig_coverage=None):
    likelihood_f = lambda x: -compute_loglikelihood(
        hist, args.read_length, args.kmer_size, x[0], x[1]
    )
    x0 = [guessed_c, guessed_e]
    res = minimize(
        likelihood_f, x0,
        bounds=((0.0, None), (0.0, 1.0)),
        options={'disp': True}
    )
    # res = minimize(likelihood_f, x0, options={'disp': True})
    cov, e = minimize_grid(likelihood_f, res.x, bounds=((0.0, None), (0.0, 1.0)))

    print('Guessed coverage:', guessed_c)
    print('Final coverage:', cov)
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
                             guessed_q1=0.5, guessed_q=0.5,
                             error_rate=None, orig_coverage=None):

    likelihood_f = lambda x: -compute_loglikelihood_with_repeats(
        hist, args.read_length, args.kmer_size, x[0], x[1], x[2], x[3],
    )
    x0 = [guessed_c, guessed_e, guessed_q1, guessed_q]

    res = minimize(
        likelihood_f, x0,
        bounds=((0.0, None), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        options={'disp': False}
    )
    # res = minimize(likelihood_f, x0, options={'disp': False})
    # cov, e, q1, q = res.x
    cov, e, q1, q = minimize_grid(
        likelihood_f, res.x,
        bounds=((0.0, None), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
    )

    print('Guessed coverage:', guessed_c)
    print('Final coverage:', cov)
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


def plot_probs(r, k, hist, est_c, est_e, guess_c, guess_e, orig_c=None, orig_e=None):
    hs = float(sum(hist))
    hp = [f / hs for f in hist]
    ep_f = compute_probabilities(r, k, est_c, est_e)
    gp_f = compute_probabilities(r, k, guess_c, guess_e)
    if orig_c is not None and orig_e is not None:
        op_f = compute_probabilities(r, k, orig_c, orig_e)
    else:
        op_f = lambda _: 0
    ep = [0] + [ep_f(j) for j in range(1, len(hist))]
    gp = [0] + [gp_f(j) for j in range(1, len(hist))]
    op = [0] + [op_f(j) for j in range(1, len(hist))]
    plt.plot(
        range(len(hp)), hp, 'ko',
        label='hist',
        ms=8,
    )
    plt.plot(
        range(len(ep)), ep, 'ro',
        label='est: C:{:.3f} E:{:.3f}'.format(est_c, est_e),
        ms=6,
    )
    plt.plot(
        range(len(gp)), gp, 'go',
        label='guess: C:{:.3f} E:{:.3f}'.format(guess_c, guess_e),
        ms=5,
    )
    plt.plot(
        range(len(op)), op, 'co',
        label='orig: C:{:.3f} E:{:.3f}'.format(orig_c, orig_e),
        ms=4,
    )
    plt.legend()
    plt.show()


def minimize_grid(fn, initial_guess, bounds=None, oprions=None):
    def generate_grid(args, step, max_depth):
        def generate_grid_single(var):
            return (
                var * step ** d
                for d in range(-max_depth, max_depth + 1) if d != 0
            )

        def filter_bounds(var_grid, i):
            if bounds is None or len(bounds) >= i or len(bounds[1]) != 2:
                return var_grid
            low, high = bounds[i]
            return (
                var for var in var_grid
                if (low is None or var >= low) and (high is None or var <= high)
            )

        var_grids = [
            list(filter_bounds(generate_grid_single(var), i))
            for i, var in enumerate(args)
        ]
        return itertools.product(*var_grids)

    min_val = fn(initial_guess)
    min_args = initial_guess
    step = 1.1
    grid_depth = 3
    diff = 1
    while (diff > 0.1 or step > 1.001):
        diff = 0
        for args in generate_grid(min_args, step, grid_depth):
            val = fn(args)
            if val < min_val:
                min_val = val
                min_args = args
                diff += min_val - val
        if diff < 1:
            step = 1 + (step - 1) * 0.9

    return min_args


def main(args):
    orig_error_rate = args.error_rate if 'error_rate' in args else None
    orig_coverage = args.coverage if 'coverage' in args else None
    all_kmers, unique_kmers, observed_ones, hist = load_dist(args.input_histogram)
    cov, e = compute_coverage_apx(
        all_kmers, unique_kmers, observed_ones,
        args.kmer_size, args.read_length
    )

    # We were unable to guess cov and e, try to estimate from some fixed valid data instead
    if cov == 0 and e == 1:
        cov = 1
        e = 0.5

    cov_est = compute_coverage_repeats if args.repeats else compute_coverage
    cov2, e2 = cov_est(
        hist, args.read_length, args.kmer_size, cov, e,
        error_rate=orig_error_rate, orig_coverage=orig_coverage,
    )
    if args.plot:
        plot_probs(
            args.read_length, args.kmer_size, hist,
            cov2, e2, cov, e, orig_coverage, orig_error_rate,
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
    parser.add_argument('--plot', action='store_true', help='Plot probabilities')
    parser.add_argument('--repeats', action='store_true', help='Estimate vith repeats')

    args = parser.parse_args()
    main(args)
