import json
import math
import random
from collections import namedtuple, defaultdict
from functools import lru_cache
from os import path

from scipy.stats import binom

from covest_poisson import poisson

from . import config
from .utils import safe_int, verbose_print


def get_trim(hist, precision=0):
    ss = float(sum(hist.values()))
    s = 0.0
    trim = None
    for i, h in sorted(hist.items()):
        s += h
        r = s / ss
        if precision:
            r = round(r, precision)
        if r == 1 and trim is None:
            trim = i
    return trim


def sample_hist(hist, factor=2, trim=None):
    verbose_print('Sampling histogram...')
    if trim is None:
        if len(hist) > 300:
            trim = get_trim(hist, 5)
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
            h[i] = math.ceil(v)
        else:
            h[i] = math.floor(v)
    return {k: v for k, v in h.items() if v > 0}   # remove 0 elements


def auto_sample_hist(hist, target_size=50, sample_factor=2):
    h = dict(hist)
    f = 1
    while max(h) > target_size:
        h = sample_hist(h, sample_factor)
        f *= sample_factor
    return h, f


class InvalidFormatException(Exception):
    def __init__(self, fname):
        self.fname = fname
    def __str__(self):
        return 'Unable to parse %s. Unsupported format.' % self.fname


def load_histogram(fname):
    hist = dict()
    with open(fname, 'r') as f:
        for line in f:
            try:
                l = line.split()
                i = int(l[0])
                cnt = int(l[1])
                hist[i] = cnt
            except (ValueError, KeyError):
                raise InvalidFormatException(fname)
    return hist


def trim_hist(hist, threshold):
    if threshold >= max(hist):
        return hist, 0
    h = {k: v for k, v in hist.items() if k < threshold}
    tail = sum(v for k, v in hist.items() if k >= threshold)
    # remove 0 elements
    return {k: v for k, v in h.items() if v > 0}, tail


def process_histogram(hist, auto_trim=None, trim=None, auto_sample=None, sample_factor=None):
    hist = dict(hist)
    tail = 0
    if auto_sample is None and sample_factor is not None:
        hist = sample_hist(hist, sample_factor, trim)
    if auto_sample is not None:
        if sample_factor is None:
            sample_factor = config.DEFAULT_SAMPLE_FACTOR
        hist, sample_factor = auto_sample_hist(hist, auto_sample, sample_factor)
        verbose_print('Histogram sampled with factor {}.'.format(sample_factor))
    if auto_trim is not None:
        trim = get_trim(hist, auto_trim)
        verbose_print('Trimming at: {}'.format(trim))
        hist, tail = trim_hist(hist, trim)
    elif trim is not None:
        hist, tail = trim_hist(hist, trim)
    hist['tail'] = tail
    return hist, sample_factor


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


def parse_data(data):
    model_class_name = data['model']

    guess = list()
    guess.append(data.get('guessed_coverage', None))
    guess.append(data.get('guessed_error_rate', None))
    if model_class_name == 'RepeatsModel':
        guess.append(data.get('guessed_q1', None))
        guess.append(data.get('guessed_q2', None))
        guess.append(data.get('guessed_q', None))

    estimated = list()
    estimated.append(data.get('estimated_coverage', None))
    estimated.append(data.get('estimated_error_rate', None))
    if model_class_name == 'RepeatsModel':
        estimated.append(data.get('estimated_q1', None))
        estimated.append(data.get('estimated_q2', None))
        estimated.append(data.get('estimated_q', None))

    return namedtuple('ParsedData', ('estimated', 'guess'))(estimated=estimated, guess=guess)


def print_output(
    hist_orig,
    model,
    estimated=None,
    guess=None,
    orig_coverage=None,
    orig_error_rate=None,
    orig_q1=None,
    orig_q2=None,
    orig_q=None,
    sample_factor=None,
    repeats=False,
    silent=False,
    reads_size=None,
):
    output_data = dict()
    output_data['model'] = model.__class__.__name__

    if guess is not None:
        output_data['guessed_loglikelihood'] = model.compute_loglikelihood(*guess)
        output_data['guessed_coverage'] = guess[0] * sample_factor
        output_data['guessed_error_rate'] = guess[1]
        if repeats:
            output_data['guessed_q1'] = guess[2]
            output_data['guessed_q2'] = guess[3]
            output_data['guessed_q'] = guess[4]

    if estimated is not None:
        output_data['estimated_loglikelihood'] = model.compute_loglikelihood(*estimated)
        output_data['estimated_coverage'] = estimated[0] * sample_factor
        output_data['estimated_error_rate'] = estimated[1]
        if repeats:
            output_data['estimated_q1'] = estimated[2]
            output_data['estimated_q2'] = estimated[3]
            output_data['estimated_q'] = estimated[4]
            if orig_q1 is None:
                orig_q1 = estimated[2]
            if orig_q2 is None:
                orig_q2 = estimated[3]
            if orig_q is None:
                orig_q = estimated[4]

        output_data['estimated_genome_size'] = safe_int(round(
            sum(
                i * h for i, h in hist_orig.items()
            ) / model.correct_c(estimated[0] * sample_factor)
        ))

        if reads_size is not None:
            output_data['estimated_genome_size_r'] = safe_int(
                round(reads_size / (estimated[0] * sample_factor))
            )

    if orig_error_rate is not None:
        output_data['original_error_rate'] = orig_error_rate
    elif estimated is not None:
        orig_error_rate = estimated[1]

    if orig_coverage is not None:
        if repeats:
            output_data['original_loglikelihood'] = model.compute_loglikelihood(
                orig_coverage, orig_error_rate, orig_q1, orig_q2, orig_q
            )
        else:
            output_data['original_loglikelihood'] = model.compute_loglikelihood(
                orig_coverage, orig_error_rate
            )

    output_data['hist_size'] = len(model.hist)
    output_data['tail_included'] = config.ESTIMATE_TAIL
    output_data['sample_factor'] = sample_factor

    if not silent:
        print(json.dumps(
            output_data, sort_keys=True, indent=4, separators=(',', ': ')
        ))

    return output_data
