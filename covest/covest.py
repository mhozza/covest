import argparse
import json
from multiprocessing import Pool

from scipy.optimize import minimize

from . import config
from .data import count_reads_size, parse_data, print_output, load_histogram
from .grid import initial_grid, optimize_grid
from .histogram import process_histogram
from .models import BasicModel, RepeatsModel
from .perf import running_time, running_time_decorator
from .utils import verbose_print


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
    hist_orig = load_histogram(args.input_histogram)
    # Process histogram and obtain first guess for c and e
    hist, tail, sample_factor, guess_c, guess_e = process_histogram(
        hist_orig, args.kmer_size, args.read_length, auto_trim=args.auto_trim,
        trim=args.trim, auto_sample=args.auto_sample, sample_factor=args.sample_factor,
    )
    reads_size = None
    if args.read_file is not None:
        reads_size = count_reads_size(args.read_file)
    err_scale = args.error_scale
    if sample_factor is None:
        sample_factor = 1
    if args.coverage:
        args.coverage /= sample_factor

    # Model selection
    if args.repeats:
        model_class = RepeatsModel
    else:
        model_class = BasicModel
    # Model initialisation
    model = model_class(
        args.kmer_size, args.read_length, hist, tail,
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
                if not (guess_c == 0 and guess_e == 1):  # We were able to guess cov and e
                    guess[:2] = guess_c, guess_e
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
                hist_orig,
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
    parser.add_argument('-at', '--auto-trim', action='store_true',
                        help='Trim histogram automatically with this threshold')
    parser.add_argument('-sf', '--sample-factor', type=int,
                        help='Sample histogram with this factor')
    parser.add_argument('-as', '--auto-sample', action='store_true',
                        help='Sample histogram automatically')
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
