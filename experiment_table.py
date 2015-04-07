#! /usr/bin/env python
import glob
import os
from collections import defaultdict

path = 'experiment1/'
files = sorted(glob.glob(os.path.join(path, '*.out')))


def parse_fname(fname):
    basename, ext = os.path.splitext(os.path.splitext(fname)[0])
    parts = basename.split('_')
    cov = parts[1][1:]
    error = parts[2][1:]
    k = parts[3][1:]
    return int(cov), float(error), int(k), ext, parts[2][0] == 'f'


def parse_khmer(fname):
    s = 'Estimated single-genome coverage is: '
    l = 2
    with open(fname) as f:
        for i, line in enumerate(f):
            if i == l:
                return float(line[len(s):])


def parse_estimate(fname):
    err = "'Error rate (est"
    cov = "'Coverage estimated from estimated p:'"
    # ', 3.1537225216529357, 4.44186270655343)
    coverage = None
    error = None
    with open(fname) as f:
        for i, line in enumerate(f):
            line = line[1:len(line) - 2]
            parts = line.split(',')
            if parts[0] == err:
                error = float(parts[2])
            if parts[0] == cov:
                coverage = float(parts[2])
    return coverage, error


def format_table(header, lines, template, line_template, header_cell_template,
                 cell_template, round_floats=2):
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
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
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

table_lines = defaultdict(dict)

for fname in files:
    cov, error, k, ext, ef = parse_fname(fname)

    key = (cov, error, k)

    table_lines[key]['original_coverage'] = cov
    table_lines[key]['original_error_rate'] = error
    table_lines[key]['original_k'] = k

    if ext == '.est':
        if ef:
            est_cov_ef, _ = parse_estimate(fname)
            table_lines[key]['estimated_ef_coverage'] = est_cov_ef
        else:
            est_cov, est_err = parse_estimate(fname)
            table_lines[key]['estimated_coverage'] = est_cov
            table_lines[key]['estimated_error_rate'] = est_err
    else:
        if ef:
            table_lines[key]['khmer_ef_coverage'] = parse_khmer(fname)
        else:
            table_lines[key]['khmer_coverage'] = parse_khmer(fname)

header = [
    'original_coverage', 'original_error_rate', 'original_k',
    'estimated_coverage', 'estimated_error_rate',
    'estimated_ef_coverage',
    'khmer_coverage',
    'khmer_ef_coverage',
]

print(format_table_html(
    header,
    sorted(
        list(table_lines.values()),
        key=lambda x: (x['original_coverage'], x['original_error_rate'], x['original_k'])
    )
))
