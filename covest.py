#!/usr/bin/env python3
import sys
import itertools
from collections import defaultdict
from scipy.misc import comb
from scipy.optimize import minimize
import argparse
from inverse import inverse
import matplotlib.pyplot as plt
from perf import running_time_decorator, running_time
from functools import lru_cache
import json
from math import exp, log
from multiprocessing import Pool  # pylint: disable=E0611
import pickle
from os import path
# from utils import print_wrap as pw

# defaults
DEFAULT_K = 21
DEFAULT_READ_LENGTH = 100
DEFAULT_THREAD_COUNT = 4
DEFAULT_REPEAT_MODEL = 1

# config
VERBOSE = True
# INF = 1e100
INF = float('inf')
USE_BIGFLOAT = False


def verbose_print(message):
    if not VERBOSE:
        return
    sys.stderr.write(message + "\n")


try:
    if not USE_BIGFLOAT:
        raise ImportError("USE_BIGFLOAT is false")
    from bigfloat import exp, log
except ImportError:
    verbose_print('BigFloats are not used!\nPrecision issues may occur.')


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def get_trim(hist, precision=0):
    ss = float(sum(hist))
    s = 0.0
    trim = None
    for i, h in enumerate(hist):
        s += h
        r = s / ss
        if precision:
            r = round(r, precision)
        if r == 1 and trim is None:
            trim = i
    return trim


def load_hist(fname, autotrim=None, trim=None):
    hist = defaultdict(int)
    max_hist = 0

    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            cnt = int(l[1])
            hist[i] = cnt
            max_hist = max(max_hist, i)

    hist_l = [hist[b] for b in range(max_hist)]
    hist_trimed = hist_l[:]
    if autotrim is not None:
        trim = get_trim(hist_l, autotrim)
        verbose_print('Trimming at: {}'.format(trim))
        hist_trimed = hist_l[:trim]
    elif trim is not None:
        hist_trimed = hist_l[:trim]
    return hist_l, hist_trimed


def count_reads_size(fname):
    from Bio import SeqIO
    _, ext = path.splitext(fname)
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    try:
        with open(fname, "rU") as f:
            return sum(len(read) for read in SeqIO.parse(f, fmt))
    except FileNotFoundError as e:
        verbose_print(e)


def compute_coverage_apx(hist, k, r):
    observed_ones = hist[1]
    all_kmers = sum(i * h for i, h in enumerate(hist))
    total_unique_kmers = sum(h for h in hist)

    if total_unique_kmers == 0:
        return 0.0, 1.0

    # discard first column
    all_kmers -= observed_ones
    unique_kmers = total_unique_kmers - observed_ones
    # compute coverage from hist >=2
    try:
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
            return float(kmer_to_read_coverage(cov / estimated_p)), float(e)
        else:
            return 0.0, float(e)
    except ZeroDivisionError:
        return 0.0, 1.0


def safe_log(x):
    if x <= 0:
        return -INF
    return log(x)


def fix_zero(x, val=1):
    if x == 0:
        return val
    else:
        return x


class BasicModel:
    def __init__(self, k, r, hist, max_error=None, max_cov=None):
        self.repeats = False
        self.k = k
        self.r = r
        self.bounds = ((0.0, max_cov), (0.0, 1.0))
        self.comb = [comb(k, s) for s in range(k + 1)]
        self.hist = hist
        if max_error is None:
            self.max_error = self.k + 1
        else:
            self.max_error = min(self.k + 1, max_error)
        self.log = [INF] + [log(j) for j in range(1, len(hist) + 1)]
        self.sum_log = [None] + [sum(self.log[1:j + 1]) for j in range(1, len(hist))]

    def correct_c(self, c):
        return c * (self.r - self.k + 1) / self.r

    def log_tr_poisson(self, l, j):
        # p_{l, j} = (l^j/j!)/(e^l - 1) = (p1 / p2) / p3
        # log(p_{l, j}) = log((l^j/j!)/(e^l - 1)) = log((p1 / p2) / p3)
        #               = log(p1) - log(p2) - log(p3)
        # compute p1 and p2 in log scale
        p1 = j * log(l)
        p2 = self.sum_log[j]
        # approxiate p3 if l is too low
        try:
            if exp(l) < 1:
                p3 = log(exp(l) - 1.0)
            else:
                p3 = log(l)
        except OverflowError as e:
            verbose_print('Owerflov at l: {}, j: {}\n{}'.format(l, j, e))
            return 0
        # compute final result
        return p1 - p2 - p3

    @lru_cache(maxsize=None)
    def get_lambda_s(self, c, err):
        return [
            c * (3 ** -s) * (1.0 - err) ** (self.k - s) * err ** s
            for s in range(self.max_error)
        ]

    def compute_probabilities(self, c, err, *_):
        # read to kmer coverage
        c = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self.get_lambda_s(c, err)
        # expected probability of kmers with s errors and coverage >= 1
        n_s = [self.comb[s] * (3 ** s) * (1.0 - exp(-l_s[s])) for s in range(self.max_error)]
        sum_n_s = fix_zero(sum(n_s[t] for t in range(self.max_error)))
        # portion of kmers with s errors
        a_s = [n_s[s] / sum_n_s for s in range(self.max_error)]
        # probability that unique kmer has coverage j (j > 0)
        max_hist = len(self.hist)
        p_j = [None] + [
            sum(
                a_s[s] * exp(self.log_tr_poisson(l_s[s], j))
                for s in range(self.max_error)
            )
            for j in range(1, max_hist)
        ]
        return p_j

    def compute_loglikelihood(self, *args):
        c, err = args[:2]
        if err < 0 or err >= 1 or c <= 0:
            return -INF
        p_j = self.compute_probabilities(*args)
        return float(sum(
            self.hist[j] * safe_log(p_j[j])
            for j in range(1, len(self.hist))
            if self.hist[j]
        ))

    def plot_probs(self, est, guess, orig):
        def fmt(p):
            return ['{:.3f}'.format(x) for x in p[:20]]

        hs = float(sum(self.hist))
        hp = [f / hs for f in self.hist]
        ep = self.compute_probabilities(*est)
        gp = self.compute_probabilities(*guess)
        if orig is not None:
            op = self.compute_probabilities(*orig)
        else:
            op = [0 for j in range(len(self.hist))]
        plt.plot(
            range(len(hp)), hp, 'ko',
            label='hist',
            ms=8,
        )
        plt.plot(
            range(len(ep)), ep, 'ro',
            label='est: {}'.format(fmt(est)),
            ms=6,
        )
        plt.plot(
            range(len(gp)), gp, 'go',
            label='guess: {}'.format(fmt(guess)),
            ms=5,
        )
        plt.plot(
            range(len(op)), op, 'co',
            label='orig: {}'.format(fmt(orig)),
            ms=4,
        )
        plt.legend()
        plt.show()


class RepeatsModel2(BasicModel):
    def __init__(self, k, r, hist, max_error=None, max_cov=None, treshold=1e-8):
        super(RepeatsModel2, self).__init__(k, r, hist, max_error)
        self.repeats = False
        self.bounds = ((0.0, max_cov), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
        self.treshold = treshold

    def get_hist_treshold(self, b_o, treshold):
        hist_size = len(self.hist)
        if treshold is not None:
            for o in range(1, hist_size):
                if b_o(o) <= treshold:
                    return o
        return hist_size

    def get_b_o(self, q1, q2, q):
        o_2 = (1 - q1) * q2
        o_n = (1 - q1) * (1 - q2) * q

        def b_o(o):
            if o == 0:
                return 0
            elif o == 1:
                return q1
            elif o == 2:
                return o_2
            else:
                return o_n * (1 - q) ** (o - 3)
        return b_o

    def compute_probabilities(self, c, err, q1, q2, q):  # pylint: disable=W0221
        b_o = self.get_b_o(q1, q2, q)
        treshold_o = self.get_hist_treshold(b_o, self.treshold)
        # read to kmer coverage
        c = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self.get_lambda_s(c, err)
        # expected probability of kmers with s errors and coverage >= 1
        n_os = [None] + [
            [self.comb[s] * (3 ** s) * (1.0 - exp(o * -l_s[s])) for s in range(self.max_error)]
            for o in range(1, treshold_o)
        ]
        sum_n_os = [None] + [
            fix_zero(sum(n_os[o][t] for t in range(self.max_error))) for o in range(1, treshold_o)
        ]

        # portion of kmers wit1h s errors
        a_os = [None] + [
            [n_os[o][s] / (sum_n_os[o] if sum_n_os[o] != 0 else 1) for s in range(self.max_error)]
            for o in range(1, treshold_o)
        ]
        # probability that unique kmer has coverage j (j > 0)
        p_j = [None] + [
            sum(
                b_o(o) * sum(
                    a_os[o][s] * exp(self.log_tr_poisson(o * l_s[s], j))
                    for s in range(self.max_error)
                )
                for o in range(1, treshold_o)
            )
            for j in range(1, len(self.hist))
        ]
        return p_j


class RepeatsModel(RepeatsModel2):
    @lru_cache(maxsize=None)
    def compute_repeat_table(self, c, err, treshold_o):
        p_j = super(RepeatsModel2, self).compute_probabilities(c, err)  # pylint: disable=E1003
        p_oj = [
            [None], p_j
        ]

        hist_size = len(self.hist)
        if treshold_o is None:
            treshold_o = hist_size

        for o in range(2, treshold_o):
            p = [[None]]
            for j in range(1, hist_size):
                res = 0.0
                for i in range(1, j):
                    t = p_oj[1][i] * p_oj[o - 1][j - i]
                    res += t
                p.append(res)
            p_oj.append(p)

        return p_oj

    def compute_probabilities(self, c, err, q1, q2, q):
        b_o = self.get_b_o(q1, q2, q)
        treshold_o = self.get_hist_treshold(b_o, self.treshold)
        p_oj = self.compute_repeat_table(c, err, treshold_o)

        p_j = [0] + [
            sum(
                b_o(o) * p_oj[o][j]
                for o in range(1, min(j + 1, treshold_o))
            )
            for j in range(1, len(self.hist))
        ]
        return p_j


class CoverageEstimator:
    def __init__(self, model, fix=None):
        self.model = model
        self.fix = fix

    def likelihood_f(self, x):
        if self.fix is not None:
            x = [j if self.fix[i] is None else self.fix[i] for i, j in enumerate(x)]
        return -self.model.compute_loglikelihood(*x)

    def compute_coverage(self, guess, use_grid=True, n_threads=1):
        r = guess
        with running_time('BFGS'):
            res = minimize(
                self.likelihood_f, r,
                bounds=self.model.bounds,
            )
            r = res.x

        if use_grid:
            verbose_print('Starting grid search with guess: {}'.format(r))
            r = optimize_grid(
                self.likelihood_f, r, bounds=self.model.bounds,
                fix=self.fix, n_threads=n_threads,
            )
        return r

    def print_output(self, estimated=None, guess=None,
                     orig_coverage=None, orig_error_rate=None,
                     orig_q1=None, orig_q2=None, orig_q=None,
                     repeats=False, silent=False, reads_size=None):
        output_data = dict()
        output_data['model'] = self.model.__class__.__name__

        if guess is not None:
            output_data['guessed_loglikelihood'] = -self.likelihood_f(guess)
            output_data['guessed_coverage'] = guess[0]
            output_data['guessed_error_rate'] = guess[1]
            if repeats:
                output_data['guessed_q1'] = guess[2]
                output_data['guessed_q2'] = guess[3]
                output_data['guessed_q'] = guess[4]

        if estimated is not None:
            output_data['estimated_loglikelihood'] = -self.likelihood_f(estimated)
            output_data['estimated_coverage'] = estimated[0]
            output_data['estimated_error_rate'] = estimated[1]
            if repeats:
                output_data['estimated_q1'] = estimated[2]
                output_data['estimated_q2'] = estimated[3]
                output_data['estimated_q'] = estimated[4]
                if orig_q1 is None:
                    orig_q1 = estimated[2]
                if orig_q2 is None:
                    orig_q2 = estimated[3]
                if orig_q is None:
                    orig_q = estimated[4]

            safe_int = lambda x: int(x) if x != float('inf') else None
            output_data['estimated_genome_size'] = safe_int(round(
                sum(
                    i * h for i, h in enumerate(self.model.hist)
                ) / self.model.correct_c(estimated[0])
            ))

            if reads_size is not None:
                output_data['estimated_genome_size_r'] = safe_int(
                    round(reads_size / estimated[0])
                )

        if orig_error_rate is not None:
            output_data['original_error_rate'] = orig_error_rate
        elif estimated is not None:
            orig_error_rate = estimated[1]

        if orig_coverage is not None:
            if repeats:
                output_data['original_loglikelihood'] = -self.likelihood_f(
                    [orig_coverage, orig_error_rate, orig_q1, orig_q2, orig_q]
                )
            else:
                output_data['original_loglikelihood'] = -self.likelihood_f(
                    [orig_coverage, orig_error_rate]
                )

        if not silent:
            print(json.dumps(
                output_data, sort_keys=True, indent=4, separators=(',', ': ')
            ))

        return output_data


def unpack_call(args):
    f, data = args
    f = pickle.loads(f)
    return f(data)


@running_time_decorator
def optimize_grid(fn, initial_guess, bounds=None, maximize=False, fix=None,
                  n_threads=DEFAULT_THREAD_COUNT):
    def generate_grid(args, step, max_depth):
        def generate_grid_single(var, fix=None):
            if fix is None:
                return (
                    var * step ** d
                    for d in range(-max_depth, max_depth + 1) if d != 0
                )
            else:
                return [fix]

        def filter_bounds(var_grid, i):
            if bounds is None or len(bounds) <= i or len(bounds[i]) != 2:
                return var_grid
            low, high = bounds[i]
            return (
                var for var in var_grid
                if (low is None or var >= low) and (high is None or var <= high)
            )

        var_grids = [
            list(filter_bounds(generate_grid_single(var, fix[i]), i))
            for i, var in enumerate(args)
        ]
        return itertools.product(*var_grids)

    if fix is None:
        fix = [None] * len(initial_guess)
    sgn = -1 if maximize else 1
    f = pickle.dumps(fn, pickle.HIGHEST_PROTOCOL)
    min_val = sgn * unpack_call([f, initial_guess])
    min_args = initial_guess
    step = 1.1
    grid_depth = 3
    diff = 1
    n_iter = 0
    try:
        while (diff > 0.1 or step > 1.001):
            n_iter += 1
            diff = 0.0
            grid = list(generate_grid(min_args, step, grid_depth))
            verbose_print('Iter : {}, Grid size: {}'.format(n_iter, len(grid)))
            fn_grid = zip([f] * len(grid), grid)
            with running_time('grid iteration'):
                with Pool(n_threads) as pool:
                    res = pool.map(unpack_call, fn_grid)
            for args, val in zip(grid, res):
                # val = fn(args)
                if sgn * val < min_val:
                    diff += min_val - val
                    min_val = sgn * val
                    min_args = args
            if diff < 1.0:
                step = 1 + (step - 1) * 0.75
            verbose_print('d:{} s:{}'.format(diff, step))
            verbose_print('New args: {}, ll: {}'.format(min_args, min_val))
    except KeyboardInterrupt:
        verbose_print('Grid search interrupted')

    verbose_print('Number of iterations in grid search:{}'.format(n_iter))
    return min_args


@running_time_decorator
def main(args):
    hist_orig, hist = load_hist(
        args.input_histogram, autotrim=args.autotrim, trim=args.trim
    )

    models = [RepeatsModel, RepeatsModel2]
    if args.repeats:
        model_class = models[args.model]
    else:
        model_class = BasicModel
    model = model_class(
        args.kmer_size, args.read_length, hist,
        max_error=8, max_cov=args.max_coverage,
    )

    orig = [args.coverage, args.error_rate, args.q1, args.q2, args.q]
    if not args.repeats:
        orig = orig[:2]

    if args.ll_only:
        ll = model.compute_loglikelihood(*orig)
        print('Loglikelihood:', ll)
    else:
        verbose_print('Estimating coverage for {}'.format(args.input_histogram))
        if args.start_original:
            cov, e, q1, q2, q = orig
        else:
            # compute guess
            cov, e = compute_coverage_apx(hist_orig, args.kmer_size, args.read_length)
            # We were unable to guess cov and e.
            # Try to estimate from some fixed valid data instead.
            if cov == 0 and e == 1:
                cov, e = 1, 0.5
            if args.fix:
                if args.coverage is not None:
                    cov = args.coverage
                if args.error_rate is not None:
                    e = args.error_rate
            q1, q2, q = [0.5 if i is None else i for i in [args.q1, args.q2, args.q]]
        if args.repeats:
            guess = [cov, e, q1, q2, q]
        else:
            guess = [cov, e]

        verbose_print('Initial guess: {} ll:{}'.format(
            guess, model.compute_loglikelihood(*guess)
        ))

        fix = [args.coverage, args.error_rate, args.q1, args.q2, args.q] if args.fix else None
        estimator = CoverageEstimator(model, fix)
        res = estimator.compute_coverage(
            guess, use_grid=args.grid, n_threads=args.thread_count
        )
        reads_size = None
        if args.genome_size is not None:
            # genome_size contains read file name
            reads_size = count_reads_size(args.genome_size)

        estimator.print_output(
            res, guess, args.coverage, args.error_rate, args.q1, args.q2, args.q,
            repeats=args.repeats, reads_size=reads_size,
        )

        if args.plot:
            model.plot_probs(
                res, guess, orig,
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', help='Input histogram')
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('--plot', action='store_true', help='Plot probabilities')
    parser.add_argument('-rp', '--repeats', action='store_true', help='Estimate with repeats')
    parser.add_argument('-ll', '--ll-only', action='store_true',
                        help='Only compute log likelihood')
    parser.add_argument('-t', '--trim', type=int, help='Trim histogram at this value')
    parser.add_argument('-M', '--max-coverage', type=int, help='Upper coverage limit')
    parser.add_argument('-at', '--autotrim', type=int, nargs='?', const=0,
                        help='Trim histogram automatically with this treshold')
    parser.add_argument('-g', '--grid', action='store_true', default=False,
                        help='Use grid search')
    parser.add_argument('-e', '--error-rate', type=float, help='Error rate')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    parser.add_argument('-q1', type=float, help='q1')
    parser.add_argument('-q2', type=float, help='q2')
    parser.add_argument('-q', type=float, help='q')
    parser.add_argument('-so', '--start-original', action='store_true',
                        help='Start form given values')
    parser.add_argument('-f', '--fix', action='store_true',
                        help='Fix some vars, optimize others')
    parser.add_argument('-m', '--model', default=DEFAULT_REPEAT_MODEL, type=int,
                        help='Model to use for estimation')
    parser.add_argument('-T', '--thread-count', default=DEFAULT_THREAD_COUNT, type=int,
                        help='Thread count')
    parser.add_argument('-s', '--genome-size', help='Calculate genome size from reads')

    main(parser.parse_args())
