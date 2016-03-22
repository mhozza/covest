from scipy.misc import comb
import matplotlib.pyplot as plt
from functools import lru_cache
import numpy
from utils import verbose_print
import config

try:
    if not config.USE_BIGFLOAT:
        raise ImportError("USE_BIGFLOAT is false")
    from bigfloat import BigFloat, exp, log
except ImportError:
    from math import exp, log
    verbose_print('BigFloats are not used!\nPrecision issues may occur.')


def fix_zero(x, val=1):
    if x == 0:
        return val
    else:
        return x


def safe_log(x):
    if x <= 0:
        return -config.INF
    return log(x)


class BasicModel:
    def __init__(self, k, r, hist, max_error=None, max_cov=None):
        self.repeats = False
        self.k = k
        self.r = r
        self.bounds = ((0.01, max_cov), (0.0, 0.5))
        self.comb = [comb(k, s) for s in range(k + 1)]
        self.hist = hist
        if max_error is None:
            self.max_error = self.k + 1
        else:
            self.max_error = min(self.k + 1, max_error)
        self.log = [config.INF] + [log(j) for j in range(1, len(hist) + 1)]

        self.factorial = [1]
        self.sum_log = [0]
        for i in range(1, len(hist)):
            t = self.factorial[-1] * i
            if config.USE_BIGFLOAT:
                t = BigFloat(t)
            else:
                t = float(t)
            self.factorial.append(t)
            self.sum_log.append(self.sum_log[-1] + self.log[i])

    def correct_c(self, c):
        return c * (self.r - self.k + 1) / self.r

    def _tr_poisson(self, l, j):
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
                    'Exception at l:{}, j:{}\n Consider histogram trimming\n{}'.format(l, j, e)
                )
                return 0.0

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
        n_s = [self.comb[s] * (3 ** s) * (1.0 - exp(-l_s[s])) for s in range(self.max_error)]
        sum_n_s = fix_zero(sum(n_s[t] for t in range(self.max_error)))
        # portion of kmers with s errors
        a_s = [n_s[s] / sum_n_s for s in range(self.max_error)]
        # probability that unique kmer has coverage j (j > 0)
        max_hist = len(self.hist)
        p_j = [None] + [
            sum(
                a_s[s] * self._tr_poisson(l_s[s], j)
                for s in range(self.max_error)
            )
            for j in range(1, max_hist)
        ]
        return p_j

    def compute_loglikelihood(self, *args):
        c, err = args[:2]
        if err < 0 or err >= 1 or c <= 0:
            return -config.INF
        p_j = self.compute_probabilities(*args)
        return float(sum(
            self.hist[j] * safe_log(p_j[j])
            for j in range(1, len(self.hist))
            if self.hist[j]
        ))

    def plot_probs(self, est, guess, orig):
        def fmt(p):
            return ['{:.3f}'.format(x) if x is not None else 'None' for x in p[:20]]

        def adjust_probs(probs):
            return [0 if p is None else i * p for i, p in enumerate(probs)]

        hs = float(sum(self.hist))
        hp = adjust_probs([f / hs for f in self.hist])
        ep = adjust_probs(self.compute_probabilities(*est))
        gp = adjust_probs(self.compute_probabilities(*guess))
        if orig is not None and None not in orig:
            op = adjust_probs(self.compute_probabilities(*orig))
        else:
            op = adjust_probs([0 for j in range(len(self.hist))])

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
        plt.show()


class RepeatsModel(BasicModel):
    def __init__(self, k, r, hist, max_error=None, max_cov=None, treshold=1e-8):
        super(RepeatsModel, self).__init__(k, r, hist, max_error)
        self.repeats = False
        self.bounds = ((0.01, max_cov), (0.0, 0.5), (0.0, 0.999), (0.0, 0.999), (0.0, 0.999))
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

    def compute_probabilities(self, c, err, q1, q2, q):
        b_o = self.get_b_o(q1, q2, q)
        treshold_o = self.get_hist_treshold(b_o, self.treshold)
        # read to kmer coverage
        c = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self._get_lambda_s(c, err)
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
                    a_os[o][s] * self._tr_poisson(o * l_s[s], j)
                    for s in range(self.max_error)
                )
                for o in range(1, treshold_o)
            )
            for j in range(1, len(self.hist))
        ]
        return p_j
