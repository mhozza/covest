import argparse
from multiprocessing import Pool
from pathlib import Path

from scipy.optimize import minimize

from covest import version_string
from . import constants
from .data import load_histogram, parse_data, print_output, save_histogram
from .grid import initial_grid, optimize_grid
from .histogram import process_histogram
from .models import models, select_model
from .perf import running_time, running_time_decorator
from .utils import nonefloat, verbose_print


class CoverageEstimator:
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
            method=constants.OPTIMIZATION_METHOD,
            bounds=self.bounds,
            options={'disp': False}
        )

    def compute_coverage(
        self,
        guess,
        starting_points=1,
        use_grid_search=False,
        n_threads=constants.DEFAULT_THREAD_COUNT,
    ):
        r = list(guess)
        r[1] *= self.err_scale
        success = True
        try:
            verbose_print('Bounds: {}'.format(self.bounds))
            if starting_points == 1:
                with running_time('First optimization'):
                    res = self._optimize(r)
                    success = res.success
                    if not success:
                        verbose_print(
                            'Optimization unsuccessful.\n'
                            'Initial params:{}\nResult{}'.format(r, res)
                        )
                    r = res.x
            elif starting_points > 1:
                params = initial_grid(r, count=starting_points, bounds=self.bounds)
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

            # If use_grid_search is none, run grid search only on failure
            if use_grid_search is None and not success:
                use_grid_search = True

            if use_grid_search:
                verbose_print('Starting grid search with guess: {}'.format(r))
                r = optimize_grid(
                    self.likelihood_f, r, bounds=self.bounds,
                    fix=self.fix, n_threads=n_threads,
                )
        except KeyboardInterrupt:
            pass
        verbose_print('Estimation finished with status: %s.' % (
            'success' if success else 'failure',
        ))
        r[1] /= self.err_scale
        return r, success


@running_time_decorator
def main(args):
    # Load histogram
    verbose_print('Loading histogram {} with parameters k={} r={}.'.format(
        args.input_histogram, args.kmer_size, args.read_length,
    ))
    hist_orig, meta = load_histogram(args.input_histogram)
    # Process histogram and obtain first guess for c and e
    hist, tail, sample_factor, guess_c, guess_e = process_histogram(
        hist_orig, args.kmer_size, args.read_length,
        trim=args.trim, sample_factor=args.sample_factor,
    )

    orig_sample_factor = 1
    if 'sample_factor' in meta:
        try:
            orig_sample_factor = int(meta['sample_factor'])
        except ValueError as e:
            print(e)
    if sample_factor > 1:
        fname = '%s.covest.sampled_x%d.hist' % (Path(args.input_histogram).stem, sample_factor)
        save_histogram(hist, fname, {
            'tool': version_string,
            'sample_factor': sample_factor * orig_sample_factor,
        })
    err_scale = args.error_scale
    if sample_factor is None:
        sample_factor = 1
    if args.coverage:
        args.coverage /= sample_factor

    # Model initialisation
    model = select_model(args.model)(
        args.kmer_size, args.read_length, hist, tail,
        max_error=constants.MAX_ERRORS, max_cov=args.max_coverage,
        min_single_copy_ratio=args.min_q1,
    )

    orig = [None] * model.param_count
    for i, v in zip(
        range(model.param_count), (args.coverage, args.error_rate) + tuple(args.params)
    ):
        orig[i] = v
    fix = orig if args.fix else None

    if args.ll_only:
        ll = model.compute_loglikelihood(*orig)
        print('Loglikelihood:', ll)
    else:
        if args.load:  # Load saved data, don't estimate anything
            with open(args.load) as f:
                parsed_data = parse_data(f)
                guess = parsed_data.guess
                res = parsed_data.estimated
        else:
            verbose_print('Estimating coverage...')
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
            guess_ll = model.compute_loglikelihood(*guess)
            if guess_ll == -constants.INF:
                verbose_print(
                    'Unable to compute likelihood. Please, try to trim the histogram, or use more complex model'
                )
                exit(1)
            verbose_print('Initial guess: {} ll:{}'.format(
                guess, guess_ll
            ))

            # Estimate coverage
            estimator = CoverageEstimator(model, err_scale=err_scale, fix=fix)
            res, success = estimator.compute_coverage(
                guess,
                starting_points=args.starting_points,
                use_grid_search=args.grid,
                n_threads=args.thread_count,
            )

            print_output(
                hist_orig, model, success, sample_factor,
                res, guess, orig,
                reads_size=args.reads_size,
                orig_sample_factor=orig_sample_factor,
            )

        # Draw plot
        if args.plot is not None:
            model.plot_probs(
                res, guess, orig, cumulative=args.plot, log_scale=constants.PLOT_LOG_SCALE
            )


def run():
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', type=str, help='Input histogram')
    parser.add_argument('-v', '--version', action='version', version=version_string,
                        help='Print version and exit.')
    parser.add_argument('-m', '--model', type=str, default='basic',
                        help='Select models for estimation. Options: {}'.format(
                            list(models.keys()))
                        )
    parser.add_argument('-k', '--kmer-size', type=int,
                        default=constants.DEFAULT_K, help='Kmer size')
    parser.add_argument('-r', '--read-length', type=int,
                        default=constants.DEFAULT_READ_LENGTH, help='Read length')
    parser.add_argument('-rs', '--reads-size', type=int,
                        help='Calculate genome size from reads size')
    parser.add_argument('-sp', '--starting-points', type=int, default=1,
                        help='Number of point to start optimization from.')
    parser.add_argument('-T', '--thread-count', default=constants.DEFAULT_THREAD_COUNT, type=int,
                        help='Thread count')
    parser.add_argument('--plot', type=bool, nargs='?', const=False,
                        help='Plot probabilities (use --plot 1 to plot probs * j)')
    parser.add_argument('--load', type=str, help='Load json')
    parser.add_argument('-t', '--trim', type=int, default=None,
                        help='Trim histogram at this value. '
                             'Set to 0 to disable automatic trimming.')
    parser.add_argument('-sf', '--sample-factor', type=int, default=None,
                        help='Use fixed sample factor for histogram sampling instead of automatic.'
                             ' Set to 1 to not sample at all.')
    parser.add_argument('-g', '--grid', action='store_true', default=False,
                        help='Use grid search for fine-tuning.')
    parser.add_argument('-f', '--fix', action='store_true',
                        help='Fix some params, optimize others')
    parser.add_argument('-c', '--coverage', type=float, help='Coverage')
    parser.add_argument('-M', '--max-coverage', type=int, help='Upper coverage limit')
    parser.add_argument('-e', '--error-rate', type=float, help='Error rate')
    parser.add_argument('-es', '--error-scale', type=float, default=constants.DEFAULT_ERR_SCALE,
                        help='Error scale')
    parser.add_argument('-mq1', '--min-q1', type=float,
                        default=constants.DEFAULT_MIN_SINGLECOPY_RATIO,
                        help='minimum single copy ratio')
    parser.add_argument('-p', '--params', type=nonefloat, nargs='*', default=tuple(),
                        help='Additional model parameters.')
    parser.add_argument('-ll', '--ll-only', action='store_true',
                        help='Only compute log likelihood from provided values')
    parser.add_argument('-so', '--start-original', action='store_true',
                        help='Start optimization form provided values')

    main(parser.parse_args())


if __name__ == '__main__':
    run()
