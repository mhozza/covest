#!/usr/bin/env python
from collections import defaultdict
import argparse
import json

from matplotlib.pyplot import bar, draw, show
import pysam

from intervaltree import IntervalTree


def compute_histogram(intervals):
    segments = dict()
    for i in intervals:
        id, start, end, rlen = i
        if not id in segments:
            segments[id] = IntervalTree(0, rlen)
        segments[id].add_to_interval(start, end)

    histogram = defaultdict(int)
    max_c = 0
    for s in segments.values():
        for i in range(s.original_end):
            c = s.get_one(i)
            histogram[c] += 1
            max_c = max(max_c, c)

    return [histogram[i] for i in range(max_c + 1)]


def compute_lengths(intervals):
    read_length = 0
    ref_length = 0
    segments = dict()
    for i in intervals:
        id, start, end, rlen = i
        segments[id] = rlen
        read_length += rlen

    ref_length = sum(segments.values())
    return ref_length, read_length


def plot_hist(hist):
    bar(range(len(hist)), hist)
    draw()
    show()


def load_sam(fname):
    mapping = pysam.AlignmentFile(fname, "r")
    for m in mapping:
        if m.reference_start >= 0:
            yield m.reference_id, m.reference_start, m.reference_end, m.rlen


def main(args):
    fname = args.input[0]
    if args.plot:
        with open(fname, 'r') as f:
            data = json.load(f)
            if 'hist' in data:
                hist = data['hist']
                lengths = data['lenghts']
            else:
                hist = data
                lengths = (0, 0)
    else:
        if args.output:
            out_basename = args.output[0]
        else:
            out_basename = fname

        intervals = load_sam(fname)
        hist = compute_histogram(intervals)
        lengths = compute_lengths(intervals)

        hist_fname = '%s.hist.json' % out_basename
        with open(hist_fname, 'w') as f:
            json.dump({'hist': hist, 'lengths': lengths}, f)

    print(hist)
    print(lengths)
    plot_hist(hist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Computes average sequence coverage')
    parser.add_argument('input', nargs=1, type=str, help='input SAM file name')
    parser.add_argument('-o', '--output', nargs=1, type=str, help='output base name')
    parser.add_argument('-p', '--plot', action='store_true', help='only plot histogram')
    args = parser.parse_args()
    main(args)
