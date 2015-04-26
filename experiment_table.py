#! /usr/bin/env python
import argparse
import glob
import os
import json
from collections import defaultdict


def parse_fname(fname, error=True):
    base, ext = os.path.splitext(os.path.splitext(fname)[0])
    base = os.path.basename(base)
    parts = base.split('_')
    if error:
        cov = parts[1][1:]
        error = parts[2][1:]
        k = parts[3][1:]
        return float(cov), float(error), int(k), ext, parts[2][0] == 'f'
    else:
        ef = False
        cov = parts[1][1:]
        if cov[-1] == 'f':
            ef = True
            cov = cov[:-1]
        k = parts[2][1:]
        return float(cov), None, int(k), ext, ef


def parse_khmer(fname):
    s = 'Estimated single-genome coverage is: '
    l = 2
    with open(fname) as f:
        for i, line in enumerate(f):
            if i == l:
                return float(line[len(s):])


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
        for i, line in enumerate(f):
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


def format_table(header, lines, template, line_template, header_cell_template,
                 cell_template, round_floats=4):
    lbody = ''
    for c in header:
        lbody += header_cell_template.format(cbody=c)
    thead = line_template.format(lbody=lbody)
    tbody = ''
    for l in lines:
        lbody = ''
        for c in header:
            val = l.get(c, None)
            if round_floats and type(val) is float:
                val = round(val, round_floats)
            lbody += cell_template.format(cbody=val)
        tbody += line_template.format(lbody=lbody)
    return template.format(thead=thead, tbody=tbody)


def format_table_html(header, lines):
    template = '''<!DOCTYPE html>
    <html>
    <head>
    <!-- Latest compiled and minified CSS -->
    <link rel="stylesheet"
        href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
    </head>
    <table class="table table-condensed table-hover">
        {thead}
        {tbody}
    </table>
    </html>
    '''
    line_template = '\t<tr>\n{lbody}\t</tr>\n'
    cell_template = '\t\t<td>{cbody}</td>\n'
    header_cell_template = '\t\t<th>{cbody}</th>\n'
    return format_table(header, lines, template, line_template,
                        header_cell_template, cell_template)


def kmer_to_read_coverage(c, k, read_length=100):
    if c is not None:
        return c * read_length / (read_length - k + 1)


def main(args):
    path = args.path
    files = sorted(glob.glob(os.path.join(path, '*.out')))
    error = not args.no_error

    table_lines = defaultdict(dict)

    for fname in files:
        cov, error, k, ext, ef = parse_fname(fname, error)
        repeats = ext[-1] == 'r'

        key = (cov, error, k, repeats)

        table_lines[key]['original_coverage'] = cov
        table_lines[key]['original_error_rate'] = error
        table_lines[key]['original_k'] = k
        table_lines[key]['repeats'] = repeats

        if ext == '.est' or ext == '.est_r':
            if args.legacy:
                d = parse_estimate(fname)
            else:
                d = parse_estimate2(fname)
            if ef:
                table_lines[key]['estimated_ef_coverage'] = d.get('estimated_coverage', None)
                table_lines[key]['guessed_ef_coverage'] = d.get('guessed_coverage', None)
            else:
                table_lines[key]['estimated_coverage'] = d.get('estimated_coverage', None)
                table_lines[key]['estimated_error_rate'] = d.get('estimated_error_rate', None)
                table_lines[key]['estimated_loglikelihood'] = d.get('estimated_loglikelihood', None)
                table_lines[key]['estimated_q1'] = d.get('estimated_q1', None)
                table_lines[key]['estimated_q'] = d.get('estimated_q', None)
                table_lines[key]['guessed_coverage'] = d.get('guessed_coverage', None)
                table_lines[key]['guessed_error_rate'] = d.get('guessed_error_rate', None)
                table_lines[key]['guessed_loglikelihood'] = d.get('guessed_loglikelihood', None)
                table_lines[key]['original_loglikelihood'] = d.get('original_loglikelihood', None)
        else:
            if ef:
                table_lines[key]['khmer_ef_coverage'] = kmer_to_read_coverage(parse_khmer(fname), k)
            else:
                table_lines[key]['khmer_coverage'] = kmer_to_read_coverage(parse_khmer(fname), k)

    # header = [
    #     'original_coverage', 'original_error_rate', 'original_k',
    #     'estimated_coverage', 'estimated_error_rate',
    #     'estimated_ef_coverage',
    #     'guessed_coverage', 'guessed_error_rate', 'guessed_ef_coverage',
    #     'original_loglikelihood', 'estimated_loglikelihood', 'guessed_loglikelihood',
    #     'khmer_coverage',
    #     'khmer_ef_coverage',
    # ]

    header = [
        'original_coverage', 'original_k', 'repeats',
        'estimated_coverage', 'estimated_error_rate',
        'estimated_q1', 'estimated_q',
        'guessed_coverage', 'guessed_error_rate',
        'estimated_loglikelihood', 'guessed_loglikelihood',
    ]

    print(format_table_html(
        header,
        sorted(
            list(table_lines.values()),
            key=lambda x: (x['original_coverage'], x['original_error_rate'], x['original_k'], x.get('repeats', False))
        )
    ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse experiment output and generate table')
    parser.add_argument('path', help='Experiment')
    parser.add_argument('-f', '--format', default='html', help='Table format')
    parser.add_argument('-ne', '--no-error', action='store_true', help='Error is unknown')
    parser.add_argument('--legacy', action='store_true', help='Run in legacy mode')
    args = parser.parse_args()
    main(args)
