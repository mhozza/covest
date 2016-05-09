#! /usr/bin/env python
import glob
import yaml
import os
import sys
from collections import defaultdict

SEPARATE_EF = True


def parse_fname(fname, error=True):
    base, ext = os.path.splitext(os.path.splitext(fname)[0])
    base = os.path.basename(base)
    parts = base.split('_')
    seq_name = parts[0]
    if error:
        ef = parts[2][0] == 'f'
        cov = float(parts[1][1:])
        e = float(parts[2][1:])
        k = int(parts[3][1:])
    else:
        ef = False
        e = None
        cov = parts[1][1:]
        if cov[-1] == 'f':
            ef = True
            cov = cov[:-1]
        k = int(parts[2][1:])
        cov = float(cov)

    if ef and SEPARATE_EF:
        e = 0.0
        return seq_name, cov, e, k, ext, False
    else:
        return seq_name, cov, e, k, ext, ef


def parse_khmer(fname):
    s = 'Estimated single-genome coverage is: '
    l = 2
    with open(fname) as f:
        for i, line in enumerate(f):
            if i == l:
                return float(line[len(s):])


def parse_williams(fname):
    s = 'gsize'
    with open(fname) as f:
        for i, line in enumerate(f):
            parts = line.split()
            if parts[0] == s:
                return float(parts[1])


def parse_estimate(fname):
    with open(fname, 'rU') as f:
        return yaml.load(f)


def kmer_to_read_coverage(c, k, read_length=100):
    if c is not None:
        return c * read_length / (read_length - k + 1)


def compute_average(table_lines, std_key_suffix='_std'):
    table_cnt = defaultdict(lambda: defaultdict(int))
    table_sum = defaultdict(lambda: defaultdict(float))
    table_avg = defaultdict(lambda: defaultdict(float))
    table_std_sum = defaultdict(lambda: defaultdict(float))
    for key, val in table_lines.items():
        for k, v in val.items():
            try:
                table_sum[key[1:]][k] += v
                table_cnt[key[1:]][k] += 1.0
            except TypeError:
                pass

    for key, val in table_sum.items():
        for k, v in val.items():
            if table_cnt[key][k] == 0:
                table_avg[key][k] = None
            else:
                table_avg[key][k] = v / table_cnt[key][k]

    for key, val in table_lines.items():
        for k, v in val.items():
            try:
                table_std_sum[key[1:]][k] += (v - table_avg[key[1:]][k]) ** 2
            except TypeError:
                pass

    for key, val in table_std_sum.items():
        for k, v in val.items():
            if table_cnt[key][k] <= 1:
                table_avg[key][k + std_key_suffix] = 0
            else:
                table_avg[key][k + std_key_suffix] = (v / (table_cnt[key][k] - 1)) ** 0.5

    return table_avg


def parse_all(path, file_filter, err=False, legacy=False):
    # path = args.path
    # files = sorted(glob.glob(os.path.join(path, args.filter)))
    files = sorted(glob.glob(os.path.join(path, file_filter)))
    # err = not args.no_error

    table_lines = defaultdict(dict)
    sequences = set()

    for fname in files:
        try:
            seq_name, cov, error, k, ext, ef = parse_fname(fname, err)
            repeats = ext[-1] == 'r'
            sequences.add(seq_name)
            key = (seq_name, cov, error, k, repeats)

            table_lines[key]['provided_coverage'] = cov
            table_lines[key]['provided_error_rate'] = error
            table_lines[key]['provided_k'] = k
            table_lines[key]['repeats'] = repeats
            table_lines[key]['fname'] = fname
            table_lines[key]['seq_name'] = seq_name

            if ext == '.est' or ext == '.est_r':
                d = parse_estimate(fname)
                table_lines[key]['coverage'] = d.get(
                    'coverage', None)
                table_lines[key]['error_rate'] = d.get(
                    'error_rate', None)
                table_lines[key]['loglikelihood'] = d.get(
                    'loglikelihood', None)
                table_lines[key]['q1'] = d.get('q1', None)
                table_lines[key]['q2'] = d.get('q2', None)
                table_lines[key]['q'] = d.get('q', None)
                table_lines[key]['guessed_coverage'] = d.get('guessed_coverage', None)
                table_lines[key]['guessed_error_rate'] = d.get('guessed_error_rate', None)
                table_lines[key]['guessed_loglikelihood'] = d.get(
                    'guessed_loglikelihood', None)
                table_lines[key]['original_loglikelihood'] = d.get(
                    'original_loglikelihood', None)
                table_lines[key]['genome_size'] = d.get(
                    'genome_size', None)
            elif ext == '.fit':
                table_lines[key]['williams_genome_size'] = parse_williams(fname)
            else:
                table_lines[key]['khmer_coverage'] = kmer_to_read_coverage(
                    parse_khmer(fname), k)
        except Exception as e:
            print('Unable to process {}\n{}'.format(fname, e), file=sys.stderr)
    return table_lines
