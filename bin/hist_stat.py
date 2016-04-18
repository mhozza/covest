#! /usr/bin/env python
from collections import defaultdict
import argparse


def load_dist(fname):
    unique_kmers = 0
    all_kmers = 0.0
    observed_ones = 0
    hist = defaultdict(int)
    max_hist = 0

    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            cnt = int(l[1])
            hist[i] = cnt
            max_hist = max(max_hist, i)
            if i == 1:
                observed_ones = cnt
            if i >= 2:
                unique_kmers += cnt
                all_kmers += i * cnt

    hist_l = [hist[h] for h in range(max_hist)]
    return all_kmers, unique_kmers, observed_ones, hist_l


def hist_stat(hist):
    z = 0
    last = -1
    for v in hist:
        if v == 0:
            z += 1
        else:
            if last == 0:
                print z
            z = 0
        last = v


def hist_port(hist, precision=None):
    ss = float(sum(hist))
    s = 0.0
    trim = None
    for i, h in enumerate(hist):
        s += h
        r = s / ss
        if precision:
            r = round(r, precision)
        if r == 1 and trim is None:
            trim = i
        print(i, r)
    return trim


def main(args):
    all_kmers, unique_kmers, observed_ones, hist = load_dist(args.input_histogram)
    t = hist_port(hist, 6)
    print 'Trim at: ', t

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('input_histogram', help='Input histogram')

    args = parser.parse_args()
    main(args)
