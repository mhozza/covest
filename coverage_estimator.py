#!/usr/bin/env python
from collections import defaultdict
import argparse
import json

from matplotlib.pyplot import bar, draw, show, axvline
import pysam

from intervaltree import IntervalTree
# from coverage_probability import all_cov_prob_dist


def compute_histogram(ref_sequences, intervals):
    segments = dict()
    read_length = 0

    for id, length in ref_sequences:
        segments[id] = IntervalTree(0, length)

    for id, start, end in intervals:
        read_length += end - start
        segments[id].add_to_interval(start, end)

    histogram = defaultdict(int)
    max_c = 0
    for s in segments.values():
        for i in range(s.original_end):
            c = s.get_one(i)
            histogram[c] += 1
            max_c = max(max_c, c)

    return [histogram[i] for i in range(max_c + 1)], read_length


def plot_hist(hist, avg=None):
    bar(range(len(hist)), hist)
    if avg is not None:
        axvline(avg, color='r')
    draw()
    show()


def load_summary(fname):
    with open(fname, 'r') as f:
        for line in f:
            prefix = 'Sequence-'
            if line.startswith(prefix):
                id, name, length = line.split('\t')
                id = int(id[len(prefix):]) - 1
                yield id, int(length)


def load_sam(fname):
    mapping = pysam.AlignmentFile(fname, "r")
    for m in mapping:
        if m.reference_start >= 0:
            yield m.reference_id, m.reference_start, m.reference_end


def main(args):
    fname = args.input[0]
    summary_fname = args.summary[0]
    if args.plot:
        with open(fname, 'r') as f:
            data = json.load(f)
            if 'hist' in data:
                hist = data['hist']
                read_length = data['length']
            else:
                hist = data
                read_length = (0, 0)
    else:
        if args.output:
            out_basename = args.output[0]
        else:
            out_basename = fname

        ref_sequences = load_summary(summary_fname)
        intervals = load_sam(fname)

        hist, read_length = compute_histogram(ref_sequences, intervals)

        hist_fname = '%s.hist.json' % out_basename
        with open(hist_fname, 'w') as f:
            json.dump({'hist': hist, 'length': read_length}, f)

    print(hist)
    print(read_length)
    if args.genome_length:
        avg = read_length / args.genome_length[0]
    else:
        avg = None
    plot_hist(hist, avg)

    if args.estimate:
        # dist = all_cov_prob_dist(read_length, hist)
        print(sum(i * k for i, k in enumerate(hist)), read_length)
        # plot_hist(dist)
        # print(dist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Computes average sequence coverage')
    parser.add_argument('input', nargs=1, type=str, help='input SAM file name')
    parser.add_argument('summary', nargs=1, type=str, help='index summary file name')
    parser.add_argument('-o', '--output', nargs=1, type=str, help='output base name')
    parser.add_argument('-p', '--plot', action='store_true', help='only plot histogram')
    parser.add_argument('-e', '--estimate', action='store_true',
                        help='estimate coverage distribution')
    parser.add_argument('--genome_length', nargs=1, type=int, help='length of genome')
    args = parser.parse_args()
    main(args)
