import math
from copy import copy
from functools import lru_cache

from scipy.stats import poisson, binom

INF = float('inf')
KMER_SIZE = 20


@lru_cache(maxsize=None)
def prod(iterable):
    res = 1
    for i in iterable:
        res *= i
    return res


@lru_cache(maxsize=None)
def factorial(n):
    math.factorial(n)


@lru_cache(maxsize=None)
def log_factorial(n):
    return math.fsum(math.log(i) for i in range(1, n + 1))


@lru_cache(maxsize=None)
def single_base_cov_prob(genome_size, read_length, read_count, exact=False):
    prob_base = read_count * (read_length - KMER_SIZE) / genome_size
    if exact:
        return binom(read_length, prob_base)
    else:
        return poisson(read_length * prob_base)


def check_log_multinomial_params(n, p_list, k_list):
    k_sum = sum(k_list)
    p_sum = math.fsum(math.exp(i) for i in p_list)
    print(n, '==', k_sum)
    assert n == k_sum
    print(1, '==', p_sum)
    assert abs(p_sum - 1) < 10**-5


def multinomial(n, p_list, k_list):
    return (
        prod(pow(p, k) for p, k in zip(p_list, k_list))
        * factorial(n)
        / prod(factorial(k) for k in k_list)
    )


def log_multinomial(n, log_p_list, k_list):
    nfact = log_factorial(n)
    kfact = math.fsum(log_factorial(k) for k in k_list)
    p_sum = math.fsum(k * p for p, k in zip(log_p_list, k_list))
    print(nfact, kfact, p_sum, nfact - kfact + p_sum)
    return nfact - kfact + p_sum


def prob_hist(genome_size, read_length, read_count, hist, log=True):
    genome_size -= KMER_SIZE
    modified_hist = copy(hist)
    covered = sum(hist[1:])
    modified_hist[0] = genome_size - covered
    if covered > genome_size:
        return -INF

    assert all((i >= 0) for i in modified_hist)

    def pmf(i):
        dist = single_base_cov_prob(genome_size, read_length, read_count)
        if log:
            return dist.logpmf(i)
        else:
            return dist.pmf(i)

    print(read_length // genome_size, read_length, genome_size)

    p_list = [
        pmf(i) for i in range(len(modified_hist))
    ]
    print(p_list)
    if log:
        # check_log_multinomial_params(genome_size, p_list, modified_hist)
        return log_multinomial(genome_size, p_list, modified_hist)
    else:
        return multinomial(genome_size, p_list, modified_hist)


def all_cov_prob_dist(total_read_length, read_count, hist, log=True):
    return (
        [prob_hist(round(total_read_length / c), total_read_length / read_count, read_count, hist, log=log) for c in (0.2, 0.4, 0.6, 0.8)] +
        [prob_hist(total_read_length // c, total_read_length / read_count, read_count, hist, log=log) for c in range(1, len(hist))]
    )
