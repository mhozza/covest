#! /usr/bin/env python
import argparse
from collections import defaultdict

from tools import templates
from tools.experiment_parser import parse_all
from tools.table_generator import format_table

SEPARATE_EF = True


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


def main(args):
    table_lines = parse_all(args.path, args.filter, not args.no_error, legacy=args.legacy)

    header = [
        'seq_name',
        'provided_coverage', 'provided_error_rate', 'provided_k',
        'coverage', 'error_rate',
        'genome_size',
        'q1', 'q2', 'q',
        'guessed_coverage', 'guessed_error_rate',
        'provided_loglikelihood', 'loglikelihood', 'guessed_loglikelihood',
    ]

    header_avg = [
        'provided_coverage',
        'provided_error_rate',
        'provided_k',
        'coverage', 'coverage_std',
        'error_rate', 'error_rate_std',
        'genome_size', 'genome_size_std',
        'q1', 'q1_std',
        'q2', 'q2_std',
        'q', 'q_std',
        'guessed_coverage', 'guessed_coverage_std',
        'guessed_error_rate', 'guessed_error_rate_std',
        'provided_loglikelihood', 'provided_loglikelihood_std',
        'loglikelihood', 'loglikelihood_std',
        'guessed_loglikelihood', 'guessed_loglikelihood_std',
    ]

    # header = [
    #     'provided_coverage', 'provided_error_rate',
    #     'coverage', 'error_rate',
    # ]

    format_templates = {
        'html': templates.html,
        'csv': templates.csv,
        'tex': templates.tex,
    }

    format_escape = {
        'tex': lambda x: x.replace('_', '\\_'),
    }

    titles = {
        'provided_coverage': 'Coverage',
        'provided_error_rate': 'Error Rate',
        'coverage': 'Est. Coverage',
        'coverage_std': 'Est. Coverage Std',
        'error_rate': 'Est. Error Rate',
        'error_rate_std': 'Est. Error Rate Std',
        'genome_size': 'Est. Genome Size',
        'genome_size_std': 'Est. Genome Size Std',
    }

    if args.average:
        table_lines = compute_average(table_lines)
        header = header_avg

    print(format_table(
        header,
        titles,
        sorted(
            list(table_lines.values()),
            key=lambda x: (
                x['provided_coverage'],
                x['provided_error_rate'],
                x['provided_k'],
                x.get('repeats', False),
                x['seq_name'],
            )
        ),
        template_file=format_templates[args.format],
        escape=format_escape.get(args.format, None),
    ))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse experiment output and generate table')
    parser.add_argument('path', help='Experiment')
    parser.add_argument('-f', '--format', default='html', help='Table format')
    parser.add_argument('-i', '--filter', default='*.out', help='Filter files')
    parser.add_argument('-a', '--average', action='store_true',
                        help='Compute average from all sequences')
    parser.add_argument('-ne', '--no-error', action='store_true', help='Error is unknown')
    parser.add_argument('--legacy', action='store_true', help='Run in legacy mode')
    args = parser.parse_args()
    main(args)
