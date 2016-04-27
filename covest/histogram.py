import random
from collections import defaultdict
from covest_poisson import poisson
from math import exp, floor, ceil

from scipy.stats import binom

from covest import config
from .utils import estimate_p, kmer_to_read_coverage, fix_coverage, verbose_print


def compute_coverage_apx(hist, k, r):
    hist = dict(hist)
    tail = hist.pop('tail', 0)
    hist[max(hist) + 1] = tail
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


def sample_histogram(hist, factor=2, trim=None):
    if trim is None:
        if len(hist) > 300:
            trim = get_trim(hist)
        else:
            trim = max(hist)
    else:
        trim = min(max(hist), trim * factor)
    hist = {k: v for k, v in hist.items() if k < trim}
    h = defaultdict(int)
    prob = 1.0 / factor
    for i, v in hist.items():
        if i < 100:
            b = binom(i, prob)
            probs = [b.pmf(j) for j in range(1, i + 1)]
        else:
            probs = [poisson(i * prob, j) for j in range(1, i + 1)]

        for j, p in enumerate(probs):
            h[j + 1] += v * p

    h = dict(h)

    for i, v in h.items():
        d = v - round(v)
        if random.random() < d:
            h[i] = ceil(v)
        else:
            h[i] = floor(v)
    return {k: v for k, v in h.items() if v > 0}   # remove 0 elements


# @Todo: return c and e
def auto_sample_hist(hist, k, r, trim=None):
    h = dict(hist)
    f = 1
    c, e = compute_coverage_apx(hist, k, r)
    while c > config.AUTO_SAMPLE_TARGET_COVERAGE:
        f += 1
        h = sample_histogram(hist, factor=f, trim=trim)
        c, e = compute_coverage_apx(h, k, r)
    return h, f


def remove_noise(hist):
    total = sum(hist.values())
    hist_denoised = {k: v for k, v in hist.items() if v/total > config.NOISE_THRESHOLD}
    tail = sum(v for v in hist.values() if v/total <= config.NOISE_THRESHOLD)
    return hist_denoised, tail


def get_trim(hist):
    hist, _ = remove_noise(hist)
    ss = float(sum(hist.values()))
    s = 0.0
    trim = max(hist)
    for i, h in sorted(hist.items()):
        s += h
        r = s / ss
        r = round(r, config.AUTO_TRIM_PRECISION)
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


def process_histogram(hist, k, r, auto_trim=True, trim=None, auto_sample=True, sample_factor=None):
    hist = dict(hist)
    tail = 0
    if sample_factor is not None:
        verbose_print('Sampling histogram...')
        hist = sample_histogram(hist, sample_factor, trim)
    if auto_sample:
        verbose_print('Sampling histogram...')
        hist, sample_factor = auto_sample_hist(hist, k, r, trim=trim)
        verbose_print('Histogram sampled with factor {}.'.format(sample_factor))
    if auto_trim:
        trim = get_trim(hist)
        verbose_print('Trimming at: {}'.format(trim))
        hist, tail = trim_hist(hist, trim)
    elif trim is not None:
        hist, tail = trim_hist(hist, trim)
    hist['tail'] = tail
    if sample_factor is None:
        sample_factor = 1
    return hist, sample_factor
