from scipy.misc import comb
import matplotlib.pyplot as plt
from perf import running_time_decorator
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

        # print(self.log)
        # print(self.sum_log)

    def correct_c(self, c):
        return c * (self.r - self.k + 1) / self.r

    def tr_poisson_old(self, l, j):
        with numpy.errstate(over='raise'):  # pylint: disable=E1101
            try:
                # if exp(l) == 1.0:  # precision fix
                #     return 0.0
                p1 = pow(l, j)
                p2 = self.factorial[j]
                if l > 1e-8:
                    p3 = exp(l) - 1.0
                else:
                    p3 = l
                res = p1 / (p2 * p3)
                return float(res)
            except (OverflowError, FloatingPointError) as e:
                verbose_print(
                    'Exception at l:{}, j:{}\n Consider histogram trimming\n{}'.format(l, j, e)
                )
                return 0.0

    def tr_poisson(self, l, j):
        if l == 0:
            return 0.0
        with numpy.errstate(over='raise'):  # pylint: disable=E1101
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
    def get_lambda_s(self, c, err):
        return [
            c * (3 ** -s) * (1.0 - err) ** (self.k - s) * err ** s
            for s in range(self.max_error)
        ]

    def compute_probabilities(self, c, err, *_):
        # read to kmer coverage
        ck = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self.get_lambda_s(ck, err)
        # expected probability of kmers with s errors and coverage >= 1
        n_s = [self.comb[s] * (3 ** s) * (1.0 - exp(-l_s[s])) for s in range(self.max_error)]
        sum_n_s = fix_zero(sum(n_s[t] for t in range(self.max_error)))
        # portion of kmers with s errors
        a_s = [n_s[s] / sum_n_s for s in range(self.max_error)]
        # probability that unique kmer has coverage j (j > 0)
        max_hist = len(self.hist)
        p_j = [None] + [
            sum(
                a_s[s] * self.tr_poisson(l_s[s], j)
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

    def der_correct_c(self, var, *args):
        if var == 0:
            return (self.r - self.k + 1) / self.r
        else:
            return 0

    def der_l(self, var, s, *args):
        c, err = args[:2]
        if var == 0:
            nzs = err ** s * (3 ** -s) if err > 0 else 0
            res = self.der_correct_c(
                var, *args
            ) * (1.0 - err) ** (self.k - s) * nzs
        elif var == 1:
            if s > 0:
                res = (
                    - self.correct_c(c) * (3 ** - s) * err ** (s - 1) *
                    (1 - err) ** (self.k - s - 1) * (self.k * err - s)
                )
            else:
                res = - self.correct_c(c) * self.k * (1 - err) ** (self.k - 1)
        else:
            res = 0
        return res

    def der_tr_poisson(self, var, s, l_s, j, *args):
        with numpy.errstate(over='raise', divide='raise'):  # pylint: disable=E1101
            try:
                p1 = self.der_l(var, s, *args)
                # verbose_print('p1_0: {}'.format(p1))
                l = l_s[s]
                if l == 0:
                    return 0
                for i in range(1, j):
                    p1 *= l / (i + 1)
                # verbose_print('p1_1: {}'.format(p1))
                while l > config.MAX_EXP and p1 > 0:
                    p1 /= exp(config.MAX_EXP)
                    l -= config.MAX_EXP
                # verbose_print('p1_2: {}'.format(p1))
                if l > 1e-8 and p1 > 0:
                    p3 = exp(l) - 1.0
                else:
                    p3 = l
                # verbose_print('p3: {}'.format(p3))
                if l_s[s] <= config.MAX_EXP:
                    p2 = j - l * exp(l) / p3
                else:
                    p2 = j - l_s[s]
                # verbose_print('p2: {}'.format(p2))
                res = (p1 / p3) * p2
                # verbose_print('res: {}'.format(res))
                return res
            except (OverflowError, FloatingPointError) as e:
                verbose_print(
                    'Exception at l_orig:{}, l:{}, j:{}\n'
                    'Consider histogram trimming\n{} L:312'.format(
                        l_s[s], l, j, e
                    )
                )
                return 0.0

    def der_tr_poisson_old(self, var, s, l_s, j, *args):
        with numpy.errstate(over='raise', divide='raise'):  # pylint: disable=E1101
            try:
                # if exp(l_s[s]) == 1.0:  # precision fix
                #     return 0.0
                p2 = self.factorial[j]
                if l_s[s] > 1e-8:
                    p3 = exp(l_s[s]) - 1.0
                else:
                    p3 = l_s[s]
                p1 = (
                    j * self.der_l(var, s, *args) * p3 * pow(l_s[s], j - 1) -
                    pow(l_s[s], j) * exp(l_s[s]) * self.der_l(var, s, *args)
                )
                res = p1 / (p2 * p3 * p3)
                return float(res)
            except (OverflowError, FloatingPointError) as e:
                verbose_print(
                    'Exception at l:{}, j:{}\n Consider histogram trimming\n{}'.format(
                        l_s[s], j, e
                    )
                )
                return 0.0

    def der_n(self, var, s, l_s, *args):
        res = self.comb[s] * (3 ** s) * exp(-l_s[s]) * self.der_l(var, s, *args)
        return res

    def der_a(self, var, s, n_s, sum_n_s, l_s, *args):
        res = (
            self.der_n(var, s, l_s, *args) * sum_n_s - n_s[s] * sum(
                self.der_n(var, s, l_s, *args)
                for s in range(self.max_error)
            )
        ) / (sum_n_s * sum_n_s)
        return res

    def der_p(self, var, *args):
        c, err = args[:2]
        # read to kmer coverage
        ck = self.correct_c(c)
        # lambda for kmers with s errors
        l_s = self.get_lambda_s(ck, err)
        # expected probability of kmers with s errors and coverage >= 1
        n_s = [self.comb[s] * (3 ** s) * (1.0 - exp(-l_s[s])) for s in range(self.max_error)]
        sum_n_s = fix_zero(sum(n_s[t] for t in range(self.max_error)))
        # portion of kmers with s errors
        a_s = [n_s[s] / sum_n_s for s in range(self.max_error)]
        # probability that unique kmer has coverage j (j > 0)
        res = [None] + [
            sum(
                self.der_a(
                    var, s, n_s, sum_n_s, l_s, *args
                ) * self.tr_poisson(l_s[s], j) + a_s[s] * self.der_tr_poisson(
                    var, s, l_s, j, *args
                )
                for s in range(self.max_error)
            )
            for j in range(1, len(self.hist))
        ]
        return res

    def der_likelihood(self, var, *args):
        p_j = self.compute_probabilities(*args)
        der_p_j = self.der_p(var, *args)

        def col(j):
            # verbose_print('{} {} {}'.format(self.hist[j], der_p_j[j], p_j[j]))
            if p_j[j] == 0:
                return config.INF
            res = self.hist[j] * der_p_j[j] / p_j[j]
            return res

        res = sum(
            col(j)
            for j in range(1, len(self.hist))
            if self.hist[j]
        )
        return res

    def gradient(self, *args):
        res = [self.der_likelihood(var, *args) for var in range(len(args))]
        return res

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


class RepeatsModel(BasicModel):
    def __init__(self, k, r, hist, max_error=None, max_cov=None, treshold=1e-8):
        super(RepeatsModel, self).__init__(k, r, hist, max_error)
        self.repeats = False
        self.bounds = ((0.01, max_cov), (0.0, 0.5), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
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

    def compute_probabilities_new(self, c, err, q1, q2, q):  # pylint: disable=W0221
        b_o = self.get_b_o(q1, q2, q)
        treshold_o = self.get_hist_treshold(b_o, self.treshold)
        p_o_j = [None] + [
            super(RepeatsModel, self).compute_probabilities(c * o, err)
            for o in range(1, treshold_o)
        ]
        p_j = [None] + [
            sum(
                b_o(o) * p_o_j[o][j]
                for o in range(1, treshold_o)
            )
            for j in range(1, len(self.hist))
        ]
        p_j_old = self.compute_probabilities_old(c, err, q1, q1, q)
        print(p_j_old)
        print(p_j)
        # assert(p_j_old == p_j)
        return p_j

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
                    a_os[o][s] * self.tr_poisson(o * l_s[s], j)
                    for s in range(self.max_error)
                )
                for o in range(1, treshold_o)
            )
            for j in range(1, len(self.hist))
        ]
        return p_j

    def der_b_o(self, var, *args):
        q1, q2, q = args[2:5]

        def b_o(o):
            if var < 2 or o == 0:
                return 0
            elif o == 1:
                return 1 if var == 2 else 0
            elif o == 2:
                return -q2 if var == 2 else (1 - q1) if var == 3 else 0
            else:
                if var == 4:
                    return (
                        (1 - q1) * (1 - q2) * (
                            (1 - q) ** (o - 3) - q * (o - 3) * (1 - q) ** (o - 4)
                        )
                    )
                else:
                    t = q2 if var == 2 else q1
                    return -(1 - t) * q * (1 - q) ** (o - 3)
        return b_o

    def der_p_o(self, var, o, *args):
        args = list(args)
        args[0] = o * args[0]
        return super(RepeatsModel, self).der_p(var, *args)

    def der_p(self, var, *args):
        c, err, q1, q2, q = args[:5]
        b_o = self.get_b_o(q1, q2, q)
        treshold_o = self.get_hist_treshold(b_o, self.treshold)
        if var < 2:
            der_p_o_j = [None] + [self.der_p_o(var, o, *args) for o in range(1, treshold_o)]
            res = [None] + [
                sum(
                    b_o(o) * der_p_o_j[o][j]
                    for o in range(1, treshold_o)
                )
                for j in range(1, len(self.hist))
            ]
        else:
            der_b_o = self.der_b_o(var, *args)
            p_o_j = [None] + [
                super(RepeatsModel, self).compute_probabilities(c * o, err)
                for o in range(1, treshold_o)
            ]
            res = [None] + [
                sum(
                    der_b_o(o) * p_o_j[o][j]
                    for o in range(1, treshold_o)
                )
                for j in range(1, len(self.hist))
            ]
        return res

    @running_time_decorator
    def gradient(self, *args):
        verbose_print('Gradient for {}'.format(args))
        res = [self.der_likelihood(var, *args) for var in range(len(args))]
        verbose_print('{}'.format(res))
        return res
