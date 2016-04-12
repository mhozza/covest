#!/usr/bin/env python3
import argparse
import itertools
import json
import pickle
import random
from collections import defaultdict, namedtuple
from functools import lru_cache
from math import exp
from multiprocessing import Pool
from os import path

from scipy.optimize import minimize

import config
from inverse import inverse
from models import BasicModel, RepeatsModel
from perf import running_time, running_time_decorator
from utils import verbose_print


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


def load_hist(fname, auto_trim=None, trim=None):
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
    hist_trimmed = hist_l[:]
    if auto_trim is not None:
        trim = get_trim(hist_l, auto_trim)
        verbose_print('Trimming at: {}'.format(trim))
        hist_trimmed = hist_l[:trim]
    elif trim is not None:
        hist_trimmed = hist_l[:trim]
    return hist_l, hist_trimmed


@lru_cache(maxsize=None)
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
    def fix_coverage(coverage):
        return inverse(lambda c: (c - c * exp(-c)) / (1 - exp(-c) - c * exp(-c)))(coverage)

    def kmer_to_read_coverage(coverage):
        return coverage * r / (r - k + 1)

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
            return float(kmer_to_read_coverage(cov / estimated_p)), float(e)
        else:
            return 0.0, float(e)
    except ZeroDivisionError:
        return 0.0, 1.0


def safe_int(x):
    return int(x) if x != float('inf') else None


class CoverageEstimator:
    GRID_SEARCH_TYPE_NONE = 0
    GRID_SEARCH_TYPE_PRE = 1
    GRID_SEARCH_TYPE_POST = 2

    def __init__(self, model, err_scale=1, fix=None):
        self.model = model
        self.fix = fix
        self.err_scale = err_scale
        self.bounds = list(self.model.bounds)
        self.bounds[1] = self.bounds[1][0], self.bounds[1][1] * self.err_scale

    def likelihood_f(self, x):
        args = list(x)
        if self.fix is not None:
            args = [j if self.fix[i] is None else self.fix[i] for i, j in enumerate(args)]
        args[1] /= self.err_scale
        return -self.model.compute_loglikelihood(*args)

    def _optimize(self, r):

        return minimize(
            self.likelihood_f, r,
            method=config.OPTIMIZATION_METHOD,
            bounds=self.bounds,
            options={'disp': False}
        )

    def compute_coverage(
            self,
            guess,
            grid_search_type=GRID_SEARCH_TYPE_PRE,
            n_threads=config.DEFAULT_THREAD_COUNT,
    ):
        r = list(guess)
        r[1] *= self.err_scale

        success = True
        try:
            verbose_print('Bounds: {}'.format(self.bounds))
            if grid_search_type == self.GRID_SEARCH_TYPE_NONE:
                with running_time('First optimization'):
                    res = self._optimize(r)
                    success = res.success
                    if not success:
                        verbose_print(
                            'Optimization unsuccessful.\n'
                            'Initial params:{}\nResult{}'.format(r, res)
                        )
                    r = res.x
            else:
                params = initial_grid(r, bounds=self.bounds)
                with running_time('Initial grid optimization'):
                    min_r = None
                    with Pool(n_threads) as pool:
                        results = list(pool.map(self._optimize, params))
                    for res in results:
                        if min_r is None or min_r > res.fun:
                            min_r = res.fun
                            success = res.success
                            if not success:
                                verbose_print(
                                    'Optimization unsuccessful.\n'
                                    'Initial params:{}\nResult{}'.format(r, res)
                                )
                            r = res.x

            if grid_search_type == self.GRID_SEARCH_TYPE_POST and not success:
                verbose_print('Starting grid search with guess: {}'.format(r))
                r = optimize_grid(
                    self.likelihood_f, r, bounds=self.bounds,
                    fix=self.fix, n_threads=n_threads,
                )
        except KeyboardInterrupt:
            pass
        r[1] /= self.err_scale
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

        output_data['hist_size'] = len(self.model.hist)

        if not silent:
            print(json.dumps(
                output_data, sort_keys=True, indent=4, separators=(',', ': ')
            ))

        return output_data


def parse_data(data):
    model_class_name = data['model']

    guess = list()
    guess.append(data.get('guessed_coverage', None))
    guess.append(data.get('guessed_error_rate', None))
    if model_class_name == 'RepeatsModel':
        guess.append(data.get('guessed_q1', None))
        guess.append(data.get('guessed_q2', None))
        guess.append(data.get('guessed_q', None))

    estimated = list()
    estimated.append(data.get('estimated_coverage', None))
    estimated.append(data.get('estimated_error_rate', None))
    if model_class_name == 'RepeatsModel':
        estimated.append(data.get('estimated_q1', None))
        estimated.append(data.get('estimated_q2', None))
        estimated.append(data.get('estimated_q', None))

    return namedtuple('ParsedData', ('estimated', 'guess'))(estimated=estimated, guess=guess)


def unpack_call(args):
    f, data = args
    f = pickle.loads(f)
    return f(data)


@running_time_decorator
def optimize_grid(fn, initial_guess, bounds=None, maximize=False, fix=None,
                  n_threads=config.DEFAULT_THREAD_COUNT):
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
    step = config.STEP
    grid_depth = config.GRID_DEPTH
    diff = 1
    n_iter = 0
    try:
        while diff > 0.1 or step > 1.001:
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


def initial_grid(initial_guess, count=config.INITIAL_GRID_COUNT, bounds=None, fix=None, ):
    def generate_grid(step):
        def apply_bounds(interval, i):
            if bounds is None or len(bounds) <= i or len(bounds[i]) != 2:
                return interval
            lb, rb = bounds[i]
            li, ri = interval
            if lb is not None:
                li = max(li, lb)
            if rb is not None:
                ri = min(ri, rb)
            return li, ri

        def generate_random_params():
            bounds = [
                apply_bounds((var / step, var * step), i) for i, var in enumerate(initial_guess)
            ]
            return [
                random.uniform(*interval) if fix is None or fix[i] is None else fix[i]
                for i, interval in enumerate(bounds)
            ]

        if count < 1:
            grid = []
        else:
            grid = [initial_guess]
            for _ in range(count - 1):
                grid.append(generate_random_params())
        return grid

    if fix is None:
        fix = [None] * len(initial_guess)
    return list(generate_grid(config.INITIAL_GRID_STEP))


@running_time_decorator
def main(args):
    hist_orig, hist = load_hist(
        args.input_histogram, auto_trim=args.auto_trim, trim=args.trim
    )
    err_scale = args.error_scale

    # Model selection
    if args.repeats:
        model_class = RepeatsModel
    else:
        model_class = BasicModel
    model = model_class(
        args.kmer_size, args.read_length, hist,
        max_error=8, max_cov=args.max_coverage,
        min_single_copy_ratio=args.min_q1,
    )

    orig = [args.coverage, args.error_rate, args.q1, args.q2, args.q]
    if not args.repeats:
        orig = orig[:2]

    if args.ll_only:
        ll = model.compute_loglikelihood(*orig)
        print('Loglikelihood:', ll)
    else:
        if args.load:
            with open(args.load) as f:
                data = json.load(f)
                parsed_data = parse_data(data)
                guess = parsed_data.guess
                res = parsed_data.estimated
        else:
            verbose_print('Estimating coverage for {}'.format(args.input_histogram))
            if args.start_original:
                if args.repeats:
                    cov, e, q1, q2, q = orig
                else:
                    cov, e = orig
            else:
                # Compute initial guess
                cov, e = compute_coverage_apx(hist_orig, args.kmer_size, args.read_length)
                if cov == 0 and e == 1:
                    # We were unable to guess cov and e.
                    # Try to estimate from some fixed valid data instead.
                    cov, e = 1, 0.5
                if args.fix:
                    if args.coverage is not None:
                        cov = args.coverage
                    if args.error_rate is not None:
                        e = args.error_rate
                q1, q2, q = [0.5 if i is None else i for i in [args.q1, args.q2, args.q]]

            # Initial guess
            if args.repeats:
                guess = [cov, e, q1, q2, q]
            else:
                guess = [cov, e]

            verbose_print('Initial guess: {} ll:{}'.format(
                guess, model.compute_loglikelihood(*guess)
            ))

            fix = [args.coverage, args.error_rate, args.q1, args.q2, args.q] if args.fix else None
            estimator = CoverageEstimator(model, err_scale=err_scale, fix=fix)
            res = estimator.compute_coverage(
                guess,
                grid_search_type=args.grid,
                n_threads=args.thread_count,
            )
            reads_size = None
            if args.genome_size is not None:
                # genome_size contains read file name
                reads_size = count_reads_size(args.genome_size)

            estimator.print_output(
                res, guess, args.coverage, args.error_rate, args.q1, args.q2, args.q,
                repeats=args.repeats, reads_size=reads_size,
            )

        if args.plot is not None:
            model.plot_probs(
                res, guess, orig, cumulative=args.plot,
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', type=str, help='Input histogram')
    parser.add_argument('--load', type=str, help='Load json')
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=config.DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=config.DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('--plot', type=bool, nargs='?', const=False,
                        help='Plot probabilities (use --plot 1 to plot probs * j)')
    parser.add_argument('-rp', '--repeats', action='store_true', help='Estimate with repeats')
    parser.add_argument('-ll', '--ll-only', action='store_true',
                        help='Only compute log likelihood')
    parser.add_argument('-t', '--trim', type=int, help='Trim histogram at this value')
    parser.add_argument('-M', '--max-coverage', type=int, help='Upper coverage limit')
    parser.add_argument('-at', '--auto-trim', type=int, nargs='?', const=0,
                        help='Trim histogram automatically with this threshold')
    parser.add_argument('-g', '--grid', type=int, default=0,
                        help='Grid search type: 0 - None, 1 - Pre-grid, 2 - Post-grid')
    parser.add_argument('-e', '--error-rate', type=float, help='Error rate')
    parser.add_argument('-es', '--error-scale', type=float, default=config.DEFAULT_ERR_SCALE,
                        help='Error scale')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    parser.add_argument('-q1', type=float, help='q1')
    parser.add_argument('-mq1', '--min-q1', type=float, default=config.DEFAULT_MIN_SINGLECOPY_RATIO,
                        help='minimum single copy ratio')
    parser.add_argument('-q2', type=float, help='q2')
    parser.add_argument('-q', type=float, help='q')
    parser.add_argument('-so', '--start-original', action='store_true',
                        help='Start form given values')
    parser.add_argument('-f', '--fix', action='store_true',
                        help='Fix some vars, optimize others')
    parser.add_argument('-T', '--thread-count', default=config.DEFAULT_THREAD_COUNT, type=int,
                        help='Thread count')
    parser.add_argument('-s', '--genome-size', help='Calculate genome size from reads')

    main(parser.parse_args())
