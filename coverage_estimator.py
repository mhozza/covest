#!/usr/bin/env python
from collections import defaultdict
import argparse
import json

from matplotlib.pyplot import bar, draw, show, axvline
import pysam

from intervaltree import IntervalTree
from coverage_probability import all_cov_prob_dist
import kmer_hist


def compute_histogram(ref_sequences, intervals):
    segments = dict()
    read_count = 0

    for id, length in ref_sequences:
        segments[id] = IntervalTree(0, length)

    for id, start, end in intervals:
        segments[id].add_to_interval(start, end)
        read_count += 1

    histogram = defaultdict(int)
    max_c = 0
    for s in segments.values():
        for i in range(s.original_end):
            c = s.get_one(i)
            histogram[c] += 1
            max_c = max(max_c, c)

    return [histogram[i] for i in range(max_c + 1)], read_count


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
            else:
                hist = data

            if 'read_count' in data:
                read_count = data['read_count']
            else:
                read_count = None
    else:
        if args.output:
            out_basename = args.output[0]
        else:
            out_basename = fname

        # ref_sequences = load_summary(summary_fname)
        # intervals = load_sam(fname)

        # hist, read_count = compute_histogram(ref_sequences, intervals)

        reads = kmer_hist.load_fastq(fname)

        hist, read_count = kmer_hist.compute_histogram(reads)

        hist_fname = '%s.hist.json' % out_basename
        with open(hist_fname, 'w') as f:
            json.dump({'hist': hist, 'read_count': read_count}, f)

    print(hist)
    total_read_length = sum(i * k for i, k in enumerate(hist))
    if args.genome_length:
        avg = total_read_length / args.genome_length[0]
    else:
        avg = None
    plot_hist(hist, avg)

    if read_count is None:
        read_count = total_read_length / 101

    if args.estimate:
        dist = all_cov_prob_dist(total_read_length, read_count, hist)
        plot_hist(dist)
        print(dist)


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
