import json
from collections import namedtuple
from functools import lru_cache
from os import path

from .utils import safe_int, verbose_print


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
    silent=False,
    reads_size=None,
):
    output_data = dict()
    output_data['model'] = model.__class__.__name__

    if guess is not None:
        output_data['guessed_loglikelihood'] = model.compute_loglikelihood(*guess)
        output_data['guessed_coverage'] = guess[0] * sample_factor
        output_data['guessed_error_rate'] = guess[1]
        if model.repeats:
            output_data['guessed_q1'] = guess[2]
            output_data['guessed_q2'] = guess[3]
            output_data['guessed_q'] = guess[4]

    if estimated is not None:
        output_data['estimated_loglikelihood'] = model.compute_loglikelihood(*estimated)
        output_data['estimated_coverage'] = estimated[0] * sample_factor
        output_data['estimated_error_rate'] = estimated[1]
        if model.repeats:
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
        if model.repeats:
            output_data['original_loglikelihood'] = model.compute_loglikelihood(
                orig_coverage, orig_error_rate, orig_q1, orig_q2, orig_q
            )
        else:
            output_data['original_loglikelihood'] = model.compute_loglikelihood(
                orig_coverage, orig_error_rate
            )

    output_data['hist_size'] = len(model.hist)
    output_data['sample_factor'] = sample_factor

    if not silent:
        print(json.dumps(
            output_data, sort_keys=True, indent=4, separators=(',', ': ')
        ))

    return output_data
