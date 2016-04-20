import math
import random
from collections import defaultdict
from covest_poisson import poisson
from scipy.stats import binom
from functools import lru_cache
from os import path

from .utils import verbose_print


def get_trim(hist, precision=0):
    ss = float(sum(hist))
    s = 0.0
    trim = None
    for i, h in enumerate(hist):
        s += h
        r = s / ss
        if precision:
            r = round(r, precision)
        if r == 1 and trim is None:
            trim = i
    return trim


def sample_hist(hist, factor=2, trim=None):
    if trim is None:
        if len(hist) > 200:
            trim = get_trim(hist, 5)
        else:
            trim = len(hist)
    else:
        trim = min(len(hist), trim*factor)

    h = [0 for _ in range(trim)]
    prob = 1.0 / factor
    max_h = 0
    for i, v in enumerate(hist):
        if i >= trim:
            break
        if i < 100:
            b = binom(i, prob)
            probs = [b.pmf(j) for j in range(1, i+1)]
        else:
            probs = [poisson(i*prob, j) for j in range(1, i + 1)]

        for j, p in enumerate(probs):
            h[j+1] += v*p

    for i, v in enumerate(h):
        d = v - round(v)
        if random.random() < d:
            h[i] = math.ceil(v)
        else:
            h[i] = math.floor(v)
        if h[i]:
            max_h = i
    return h[:max_h+1]


def auto_sample_hist(hist, target_size=50, sample_factor=2):
    h = list(hist)
    f = 1
    while len(h) > target_size:
        h = sample_hist(h, sample_factor)
        f *= sample_factor
    return h, f


def load_hist(fname, tail_sum=False, auto_trim=None, trim=None, auto_sample=None, sample_factor=None):
    hist = defaultdict(int)
    max_hist = 0

    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            cnt = int(l[1])
            hist[i] = cnt
            max_hist = max(max_hist, i)

    hist_l = [hist[b] for b in range(max_hist+1)]
    if auto_sample is None and sample_factor is not None:
        hist_l = sample_hist(hist_l, sample_factor, trim)
    if auto_sample is not None:
        if sample_factor is None:
            sample_factor = config.DEFAULT_SAMPLE_FACTOR
        hist_l, sample_factor = auto_sample_hist(hist_l, auto_sample, sample_factor)
        verbose_print('Histogram sampled with factor {}.'.format(sample_factor))
    hist_trimmed = list(hist_l)
    if auto_trim is not None:
        trim = get_trim(hist_l, auto_trim)
        verbose_print('Trimming at: {}'.format(trim))
        hist_trimmed = hist_l[:trim]
    elif trim is not None:
        hist_trimmed = hist_l[:trim]
    if tail_sum:
        tail = sum(hist_l[trim:]) if trim is not None else 0
        hist_trimmed.append(tail)
    return hist_l, hist_trimmed, sample_factor


@lru_cache(maxsize=None)
def count_reads_size(fname):
    from Bio import SeqIO
    _, ext = path.splitext(fname)
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    try:
        with open(fname, "rU") as f:
            return sum(len(read) for read in SeqIO.parse(f, fmt))
    except FileNotFoundError as e:
        verbose_print(e)
