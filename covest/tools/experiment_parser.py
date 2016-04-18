#! /usr/bin/env python
import sys
import glob
import os
import json
from collections import defaultdict
from .table_generator import format_table

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
    err = "'Error rate (est"
    err2 = "'Error rate:'"
    err3 = "'Final error rate:'"
    cov = "'Coverage estimated from estimated p:'"
    cov2 = "'Final coverage:'"
    gcov = "'Guessed coverage:'"
    gerr = "'Guessed error rate:'"
    orig_l = "'Original loglikelihood:'"
    guessed_l = "'Guessed loglikelihood:'"
    est_l = "'Estimated loglikelihood:'"
    est_q = "'Estimated q1 and q:'"
    # ', 3.1537225216529357, 4.44186270655343)
    coverage = None
    error = None
    guessed_c = None
    guessed_e = None
    orig_likelihood = None
    guessed_likelihood = None
    est_likelihood = None
    d = dict()
    with open(fname) as f:
        for _, line in enumerate(f):
            line = line[1:len(line) - 2]
            parts = line.split(',')
            if parts[0] == err:
                error = float(parts[2])
            if parts[0] == err2 or parts[0] == err3:
                error = float(parts[1])
            if parts[0] == cov:
                coverage = float(parts[2])
            if parts[0] == cov2:
                coverage = float(parts[1])
            if parts[0] == gcov:
                guessed_c = float(parts[1].strip(" '\""))
            if parts[0] == gerr:
                guessed_e = float(parts[1].strip(" '\""))
            if parts[0] == orig_l:
                orig_likelihood = float(parts[1])
            if parts[0] == guessed_l:
                guessed_likelihood = float(parts[1])
            if parts[0] == est_l:
                est_likelihood = float(parts[1])
            if parts[0] == est_q:
                d['estimated_q1'] = float(parts[1])
                d['estimated_q'] = float(parts[2])

    d['estimated_coverage'] = coverage
    d['estimated_error_rate'] = error
    d['estimated_loglikelihood'] = est_likelihood
    d['guessed_coverage'] = guessed_c
    d['guessed_error_rate'] = guessed_e
    d['guessed_loglikelihood'] = guessed_likelihood
    d['original_loglikelihood'] = orig_likelihood
    return d


def parse_estimate2(fname):
    with open(fname, 'rU') as f:
        return json.load(f)


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


def parse_all(path, filter, err=False, legacy=False):
    # path = args.path
    # files = sorted(glob.glob(os.path.join(path, args.filter)))
    files = sorted(glob.glob(os.path.join(path, filter)))
    # err = not args.no_error

    table_lines = defaultdict(dict)
    sequences = set()

    for fname in files:
        try:
            seq_name, cov, error, k, ext, ef = parse_fname(fname, err)
            repeats = ext[-1] == 'r'
            sequences.add(seq_name)
            key = (seq_name, cov, error, k, repeats)

            table_lines[key]['original_coverage'] = cov
            table_lines[key]['original_error_rate'] = error
            table_lines[key]['original_k'] = k
            table_lines[key]['repeats'] = repeats
            table_lines[key]['fname'] = fname
            table_lines[key]['seq_name'] = seq_name

            if ext == '.est' or ext == '.est_r':
                if legacy:
                    d = parse_estimate(fname)
                else:
                    d = parse_estimate2(fname)
                if ef:
                    table_lines[key]['estimated_ef_coverage'] = d.get('estimated_coverage', None)
                    table_lines[key]['guessed_ef_coverage'] = d.get('guessed_coverage', None)
                else:
                    table_lines[key]['estimated_coverage'] = d.get(
                        'estimated_coverage', None)
                    table_lines[key]['estimated_error_rate'] = d.get(
                        'estimated_error_rate', None)
                    table_lines[key]['estimated_loglikelihood'] = d.get(
                        'estimated_loglikelihood', None)
                    table_lines[key]['estimated_q1'] = d.get('estimated_q1', None)
                    table_lines[key]['estimated_q2'] = d.get('estimated_q2', None)
                    table_lines[key]['estimated_q'] = d.get('estimated_q', None)
                    table_lines[key]['guessed_coverage'] = d.get('guessed_coverage', None)
                    table_lines[key]['guessed_error_rate'] = d.get('guessed_error_rate', None)
                    table_lines[key]['guessed_loglikelihood'] = d.get(
                        'guessed_loglikelihood', None)
                    table_lines[key]['original_loglikelihood'] = d.get(
                        'original_loglikelihood', None)
                    table_lines[key]['estimated_genome_size'] = d.get(
                        'estimated_genome_size', None)
            elif ext == '.fit':
                table_lines[key]['williams_genome_size'] = parse_williams(fname)
            else:
                if ef:
                    table_lines[key]['khmer_ef_coverage'] = kmer_to_read_coverage(
                        parse_khmer(fname), k)
                else:
                    table_lines[key]['khmer_coverage'] = kmer_to_read_coverage(
                        parse_khmer(fname), k)
        except Exception as e:
            print('Unable to process {}\n{}'.format(fname, e), file=sys.stderr)
    return table_lines
