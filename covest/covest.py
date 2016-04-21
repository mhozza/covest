import argparse
import json
from math import exp
from multiprocessing import Pool

from scipy.optimize import minimize

from . import config
from .grid import initial_grid, optimize_grid
from .inverse import inverse
from .data import count_reads_size, load_hist, parse_data, print_output
from .models import BasicModel, RepeatsModel
from .perf import running_time, running_time_decorator
from .utils import verbose_print


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


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
        verbose_print('Status: %s' % ('success' if success else 'failure',))
        r[1] /= self.err_scale
        return r


@running_time_decorator
def main(args):
    # Load histogram
    hist_orig, hist, sample_factor = load_hist(
        args.input_histogram, tail_sum=config.ESTIMATE_TAIL, auto_trim=args.auto_trim,
        trim=args.trim, auto_sample=args.auto_sample, sample_factor=args.sample_factor,
    )
    reads_size = None
    if args.read_file is not None:
        reads_size = count_reads_size(args.read_file)
    err_scale = args.error_scale
    if sample_factor and args.coverage:
        args.coverage /= sample_factor

    # Model selection
    if args.repeats:
        model_class = RepeatsModel
    else:
        model_class = BasicModel
    # Model initialisation
    model = model_class(
        args.kmer_size, args.read_length, hist,
        max_error=8, max_cov=args.max_coverage,
        min_single_copy_ratio=args.min_q1,
    )

    orig = (args.coverage, args.error_rate, args.q1, args.q2, args.q)[:model.param_count]
    fix = orig if args.fix else None

    if args.ll_only:
        ll = model.compute_loglikelihood(*orig)
        print('Loglikelihood:', ll)
    else:
        if args.load:  # Load saved data, don't estimate anything
            with open(args.load) as f:
                parsed_data = parse_data(json.load(f))
                guess = parsed_data.guess
                res = parsed_data.estimated
        else:
            verbose_print('Estimating coverage for {}'.format(args.input_histogram))
            # Compute initial guess
            if args.start_original:
                guess = list(orig)
            else:
                guess = list(model.defaults)
                cov, e = compute_coverage_apx(hist_orig, args.kmer_size, args.read_length)
                if not (cov == 0 and e == 1):  # We were able to guess cov and e
                    guess[:2] = cov, e
                if fix:
                    for i, v in fix:
                        if v is not None:
                            guess[i] = v
            verbose_print('Initial guess: {} ll:{}'.format(
                guess, model.compute_loglikelihood(*guess)
            ))

            # Estimate coverage
            estimator = CoverageEstimator(model, err_scale=err_scale, fix=fix)
            res = estimator.compute_coverage(
                guess,
                grid_search_type=args.grid,
                n_threads=args.thread_count,
            )

            print_output(
                model,
                res, guess, args.coverage, args.error_rate, args.q1, args.q2, args.q,
                sample_factor=sample_factor, repeats=args.repeats, reads_size=reads_size,
            )

        # Draw plot
        if args.plot is not None:
            model.plot_probs(
                res, guess, orig, cumulative=args.plot, log_scale=config.PLOT_LOG_SCALE
            )


def run():
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
    parser.add_argument('-M', '--max-coverage', type=int, help='Upper coverage limit')
    parser.add_argument('-t', '--trim', type=int, help='Trim histogram at this value')
    parser.add_argument('-at', '--auto-trim', type=int, nargs='?', const=0,
                        help='Trim histogram automatically with this threshold')
    parser.add_argument('-sf', '--sample-factor', type=int,
                        help='Sample histogram with this factor')
    parser.add_argument('-as', '--auto-sample', type=int, nargs='?',
                        const=config.DEFAULT_SAMPLE_TARGET,
                        help='Sample histogram automatically to this target size')
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
    parser.add_argument('-s', '--genome-size', dest='read_file',
                        help='Calculate genome size from reads')

    main(parser.parse_args())


if __name__ == '__main__':
    run()
