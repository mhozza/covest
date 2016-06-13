import random
from collections import defaultdict
from covest_poisson import poisson_dist
from functools import partial
from math import exp, floor, ceil
from multiprocessing.pool import Pool, ThreadPool

from scipy.stats import binom

from covest import constants
from .utils import estimate_p, kmer_to_read_coverage, fix_coverage, verbose_print


def compute_coverage_apx(hist, k, r):
    observed_ones = hist.get(1, 0)
    all_kmers = sum(i * h for i, h in hist.items())
    total_unique_kmers = sum(h for h in hist.values())

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
            return float(kmer_to_read_coverage(cov / estimated_p, k, r)), float(e)
        else:
            return 0.0, float(e)
    except ZeroDivisionError:
        return 0.0, 1.0


def compute_binomial_probs(prob, i):
    if i < 100:
        b = binom(i, prob)
        return [b.pmf(j) for j in range(1, i + 1)]
    else:
        return poisson_dist(i * prob, i)


def compute_hist_item_from_probs(hist, prob_list, j):
    return sum(prob_list[i][j] * v for i, v in enumerate(hist.values()) if j < len(prob_list[i]))


def sample_histogram(hist, factor=2, trim=None, n_threads=constants.DEFAULT_THREAD_COUNT):
    if trim is None:
        if len(hist) > 300:
            trim = get_trim(hist)
        else:
            trim = max(hist)
    else:
        trim = min(max(hist), trim * factor)

    hist = {k: v for k, v in hist.items() if k < trim}
    compute_probs = partial(compute_binomial_probs, 1.0 / factor)
    # probs_list = map(compute_probs, hist.keys())
    with ThreadPool(n_threads) as pool:
        probs_list = pool.map(compute_probs, hist.keys())

        h = list(pool.map(partial(compute_hist_item_from_probs, hist, probs_list), range(trim)))

    for i, v in enumerate(h):
        d = v - round(v)
        h[i] = ceil(v) if random.random() < d else floor(v)
    return {k + 1: v for k, v in enumerate(h) if v > 0}   # remove 0 elements


def auto_sample_hist(hist, k, r, trim=None):
    h = dict(hist)
    f = 1
    s = 1
    c, e = compute_coverage_apx(hist, k, r)

    while c > constants.AUTO_SAMPLE_TARGET_COVERAGE:
        f += s
        s *= 2
        h = sample_histogram(hist, factor=f, trim=trim)
        c, e = compute_coverage_apx(h, k, r)

    s //= 4
    f2 = f - s
    while s >= 1:
        h2 = sample_histogram(hist, factor=f2, trim=trim)
        c, e = compute_coverage_apx(h2, k, r)
        if c > constants.AUTO_SAMPLE_TARGET_COVERAGE:
            f2 += s
        else:
            h = h2
            f = f2
            f2 -= s
        s //= 2

    return h, f, c, e


def remove_noise(hist):
    total = sum(hist.values())
    hist_denoised = {k: v for k, v in hist.items() if v / total > constants.NOISE_THRESHOLD}
    return hist_denoised


def get_trim(hist):
    hist = remove_noise(hist)
    ss = float(sum(hist.values()))
    s = 0.0
    trim = max(hist)
    for i, h in sorted(hist.items()):
        s += h
        r = s / ss
        r = round(r, constants.AUTO_TRIM_PRECISION)
        if r == 1:
            trim = i
            break
    return trim


def trim_hist(hist, threshold):
    if threshold >= max(hist):
        return hist, 0
    h = {k: v for k, v in hist.items() if k < threshold}
    tail = sum(v for k, v in hist.items() if k >= threshold)
    # remove 0 elements
    return {k: v for k, v in h.items() if v > 0}, tail


def process_histogram(hist, k, r, trim=None, sample_factor=None):
    hist = dict(hist)
    tail = 0
    if sample_factor is not None and sample_factor > 1:
        verbose_print('Sampling histogram...')
        hist = sample_histogram(hist, sample_factor, trim)
    if sample_factor is None:
        verbose_print('Sampling histogram...')
        hist, sample_factor, c, e = auto_sample_hist(hist, k, r, trim=trim)
        verbose_print('Histogram sampled with factor {}.'.format(sample_factor))
    else:
        c, e = compute_coverage_apx(hist, k, r)
    if trim is None:
        trim = get_trim(hist)
        verbose_print('Trimming at: {}'.format(trim))
        hist, tail = trim_hist(hist, trim)
    elif trim > 0:
        verbose_print('Trimming at: {}'.format(trim))
        hist, tail = trim_hist(hist, trim)
    return hist, tail, sample_factor, c, e
