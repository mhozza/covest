import math
from functools import lru_cache

from scipy.stats import poisson


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
def single_base_cov_prob(genome_size, read_length):
    prob_base = 1 / genome_size
    # return poisson(read_length * prob_base)
    return prob_base


def multinomial(n, p_list, k_list):
    return (
        prod(pow(p, k) for p, k in zip(p_list, k_list))
        * factorial(n)
        / prod(factorial(k) for k in k_list)
    )


def log_multinomial(n, log_p_list, k_list):
    return (
        math.fsum(k * p for p, k in zip(log_p_list, k_list))
        + log_factorial(n)
        - math.fsum(log_factorial(k) for k in k_list)
    )


def prob_hist(genome_size, read_length, hist, log=True):
    def pmf(i):
        if log:
            return math.log(single_base_cov_prob(genome_size, read_length))#.logpmf(i)
        else:
            return single_base_cov_prob(genome_size, read_length)#.pmf(i)

    print (read_length // genome_size, read_length, genome_size)

    p_list = [
        pmf(i) for i in range(len(hist))
    ]
    print(p_list)
    if log:
        return log_multinomial(read_length, p_list, hist)
    else:
        return multinomial(read_length, p_list, hist)


def all_cov_prob_dist(read_length, hist, log=True):
    return [prob_hist(read_length // c, read_length, hist, log=log) for c in range(1, len(hist))]
