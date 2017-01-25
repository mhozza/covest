import itertools
import pickle
import random
from multiprocessing import Pool

from . import constants
from .perf import running_time, running_time_decorator
from .utils import verbose_print


def unpack_call(args):
    f, data = args
    f = pickle.loads(f)
    return f(data)


@running_time_decorator
def optimize_grid(fn, initial_guess, bounds=None, maximize=False, fix=None,
                  n_threads=constants.DEFAULT_THREAD_COUNT):
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
    step = constants.STEP
    grid_depth = constants.GRID_DEPTH
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


def initial_grid(initial_guess, count=constants.INITIAL_GRID_COUNT, bounds=None, fix=None):
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
    return list(generate_grid(constants.INITIAL_GRID_STEP))
