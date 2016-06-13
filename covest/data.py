from collections import namedtuple
from functools import lru_cache
from os import path

import yaml
from first import first

from covest import __version__
from .models import models, BasicModel
from .utils import safe_int, verbose_print


class InvalidFormatException(Exception):
    def __init__(self, fname):
        self.fname = fname

    def __str__(self):
        return 'Unable to parse %s. Unsupported format.' % self.fname


def load_histogram(fname):
    hist = dict()
    meta = dict()
    with open(fname, 'r') as f:
        for line in f:
            if line[0] == '#':
                try:
                    k, v = line[1:].strip().split(':')
                    meta[k] = v
                except IndexError:
                    pass
            else:
                try:
                    l = line.split()
                    i = int(l[0])
                    cnt = int(l[1])
                    hist[i] = cnt
                except (ValueError, KeyError):
                    raise InvalidFormatException(fname)
    return hist, meta


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


def parse_data(f):
    data = yaml.load(f)
    try:
        model_class_name = data['model']
        model = first(models.values(), key=lambda x: x.__name__ == model_class_name)
    except KeyError:
        model = BasicModel
    guess = [
        data.get('guessed_coverage', None),
        data.get('guessed_error_rate', None),
    ]
    estimated = [data.get(k, None) for k in model.params]
    return namedtuple('ParsedData', ('estimated', 'guess', 'model'))(
        estimated=estimated, guess=guess, model=model
    )


def replace_none(dest, src):
    if dest is None or src is None:
        raise ValueError('Invalid arguments.')
    dest = list(dest)
    if len(dest) != len(src):
        raise ValueError('Length of arguments should be equal.')
    for i in range(len(dest)):
        if dest[i] is None:
            dest[i] = src[i]
    return dest


def print_output(
    hist_orig,
    model,
    success,
    sample_factor,
    estimated=None,
    guess=None,
    orig=None,
    reads_size=None,
    silent=False,
):
    def params_to_dict(names, values):
        nonlocal sample_factor
        if values is None or names is None:
            return dict()
        values = [float(v) if v is not None else v for v in values]
        # apply sample factor to coverage
        if values[0] is not None and sample_factor is not None:
            values[0] *= sample_factor
        data = dict()
        for k, v in zip(names, values):
            if v is not None:
                data[k] = v
        return data

    output_data = {
        'model': model.short_name(),
        'hist_size': max(model.hist),
        'sample_factor': sample_factor,
        'success': success,
        'version': __version__,
    }

    if guess is not None:
        output_data.update(params_to_dict(('guessed_coverage', 'guessed_error_rate'), guess))
        output_data['guessed_loglikelihood'] = model.compute_loglikelihood(*guess)
    if estimated is not None:
        output_data.update(params_to_dict(model.params, estimated))
        output_data['loglikelihood'] = model.compute_loglikelihood(*estimated)
        output_data['genome_size'] = safe_int(round(
            sum(
                i * h for i, h in hist_orig.items()
            ) / model.correct_c(estimated[0] * sample_factor)
        ))
        if reads_size is not None:
            output_data['genome_size_reads'] = safe_int(
                round(reads_size / (estimated[0] * sample_factor))
            )
    if orig is not None and any(orig):
        output_data.update(params_to_dict(('provided_%s' % name for name in model.params), orig))
        try:
            output_data['provided_loglikelihood'] = model.compute_loglikelihood(
                *replace_none(orig, estimated)
            )
        except ValueError:
            pass

    if not silent:
        print(yaml.dump(output_data, indent=4, default_flow_style=False))

    return output_data


def save_histogram(hist, fname, meta=None):
    with open(fname, 'w') as f:
        if meta:
            for k, v in meta.items():
                f.write('#{}:{}\n'.format(k, v))
        for k, v in hist.items():
            f.write('%d %d\n' % (k, v))
