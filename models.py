import itertools
import multiprocessing
from functools import lru_cache
from math import exp, log

import matplotlib.pyplot as plt
import numpy
from scipy.misc import comb

import config
from utils import verbose_print


def fix_zero(x, val=1):
    if x == 0:
        return val
    else:
        return x


def safe_log(x):
    if x is None or x <= 0:
        return -config.INF
    return log(x)


def tr_poisson(l, j):
    if l == 0:
        return 0.0
    with numpy.errstate(over='raise'):
        try:
            p1 = 1.0
            for i in range(1, j + 1):
                p1 *= l / i
            while l > config.MAX_EXP and p1 > 0:
                p1 /= exp(config.MAX_EXP)
                l -= config.MAX_EXP
            if l > 1e-8 and p1 > 0:
                p3 = exp(l) - 1.0
            else:
                p3 = l
            res = p1 / p3
            return res
        except (OverflowError, FloatingPointError) as e:
            verbose_print(
                'Precision error at l:{}, j:{}\n Consider histogram trimming\n{}'.format(
                    l, j, e
                )
            )
            return 0.0


def tail_prob(l, j):
    last_p = 1.0
    sp = 0
    for i in range(1, j + 1):
        p = last_p * (l / i)
        sp += p
        last_p = p
    return 1.0 - exp(-l) * sp


class BasicModel:
    def __init__(self, k, r, hist, max_error=None, max_cov=None, *args, **kwargs):
        self.repeats = False
        self.k = k
        self.r = r
        self.bounds = ((0.01, max_cov), (0.0, 0.5))
        self.comb = [comb(k, s) * (3 ** s) for s in range(k + 1)]
        self.hist = hist
        if max_error is None:
            self.max_error = self.k + 1
        else:
            self.max_error = min(self.k + 1, max_error)

    def check_bounds(self, args):
        for arg, (l, r) in zip(args, self.bounds):
            if arg is None:
                continue
            if arg == float('NaN'):
                return False
            if (l is not None and arg < l) or (r is not None and arg > r):
                return False
        return True

    def correct_c(self, c):
        return c * (self.r - self.k + 1) / self.r

    @lru_cache(maxsize=None)
    def _get_lambda_s(self, c, err):
        return [
            c * (3 ** -s) * (1.0 - err) ** (self.k - s) * err ** s
            for s in range(self.max_error)
        ]

    def compute_probabilities(self, c, err, *_):
        # read to kmer coverage
        ck = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self._get_lambda_s(ck, err)
        # expected probability of kmers with s errors and coverage >= 1
        n_s = [self.comb[s] * (1.0 - exp(-l_s[s])) for s in range(self.max_error)]
        sum_n_s = fix_zero(sum(n_s[t] for t in range(self.max_error)))
        # portion of kmers with s errors
        a_s = [n_s[s] / sum_n_s for s in range(self.max_error)]
        # probability that unique kmer has coverage j (j > 0)
        max_hist = len(self.hist)
        p_j = [None] + [
            sum(
                a_s[s] * tr_poisson(l_s[s], j) for s in range(self.max_error)
            )
            for j in range(1, max_hist)
        ]
        return p_j

    def compute_loglikelihood(self, *args):
        if not self.check_bounds(args):
            return -config.INF
        p_j = self.compute_probabilities(*args)
        if config.ESTIMATE_TAIL:
            sp_j = sum(p_j[1:-1])
            p_j[-1] = 1 - sp_j
        return float(sum(
            self.hist[j] * safe_log(p_j[j]) for j in range(1, len(self.hist)) if self.hist[j]
        ))

    def compute_loglikelihood_multi(self, args_list, thread_count=config.DEFAULT_THREAD_COUNT):
        if thread_count is None:  # do not use multiprocessing
            likelihoods = itertools.starmap(self.compute_loglikelihood, args_list)
        else:
            pool = multiprocessing.Pool(processes=thread_count)
            likelihoods = pool.starmap(self.compute_loglikelihood, args_list)
        return {
            tuple(args): likelihood for args, likelihood in zip(args_list, likelihoods)
        }

    def plot_probs(self, est, guess, orig, cumulative=False):
        def fmt(p):
            return ['{:.3f}'.format(x) if x is not None else 'None' for x in p[:20]]

        def adjust_probs(probs):
            if cumulative:
                return [0 if p is None else i * p for i, p in enumerate(probs)]
            else:
                return probs

        hs = float(sum(self.hist))
        hp = adjust_probs([f / hs for f in self.hist])
        ep = adjust_probs(self.compute_probabilities(*est))
        gp = adjust_probs(self.compute_probabilities(*guess))
        if orig is not None and None not in orig:
            op = adjust_probs(self.compute_probabilities(*orig))
        else:
            op = adjust_probs([0 for _ in range(len(self.hist))])

        plt.yscale('log')
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
        try:
            plt.show()
        except KeyboardInterrupt:
            pass


class RepeatsModel(BasicModel):
    def __init__(self, k, r, hist, max_error=None, max_cov=None, threshold=1e-8,
                 min_single_copy_ratio=0.3, *args, **kwargs):
        super(RepeatsModel, self).__init__(k, r, hist, max_error=max_error)
        self.repeats = True
        self.bounds = (
            (0.01, max_cov), (0.0, 0.5), (min_single_copy_ratio, 0.9999), (0.0, 0.99), (0.0, 0.99))
        self.threshold = threshold

    def get_hist_threshold(self, b_o, threshold):
        hist_size = len(self.hist)
        if threshold is not None:
            for o in range(1, hist_size):
                if b_o(o) <= threshold:
                    return o
        return hist_size

    @staticmethod
    def get_b_o(q1, q2, q):
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

    # noinspection PyMethodOverriding
    def compute_probabilities(self, c, err, q1, q2, q, *_):
        b_o = self.get_b_o(q1, q2, q)
        threshold_o = self.get_hist_threshold(b_o, self.threshold)
        # read to kmer coverage
        c = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self._get_lambda_s(c, err)
        # expected probability of kmers with s errors and coverage >= 1
        # noinspection PyTypeChecker
        n_os = [None] + [
            [self.comb[s] * (1.0 - exp(o * -l_s[s])) for s in range(self.max_error)]
            for o in range(1, threshold_o)
        ]
        sum_n_os = [None] + [
            fix_zero(sum(n_os[o][t] for t in range(self.max_error))) for o in range(1, threshold_o)
        ]

        # portion of kmers wit1h s errors
        # noinspection PyTypeChecker
        a_os = [None] + [
            [n_os[o][s] / (sum_n_os[o] if sum_n_os[o] != 0 else 1) for s in range(self.max_error)]
            for o in range(1, threshold_o)
        ]
        # probability that unique kmer has coverage j (j > 0)
        p_j = [None] + [
            sum(
                b_o(o) * sum(
                    a_os[o][s] * tr_poisson(o * l_s[s], j) for s in range(self.max_error)
                ) for o in range(1, threshold_o)
            ) for j in range(1, len(self.hist))
        ]
        return p_j
