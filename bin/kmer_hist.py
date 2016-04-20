#!/usr/bin/env python
import argparse
import random
from collections import defaultdict
from os import path

from Bio import SeqIO

NS_IGNORE = 0
NS_SINGLE = 1
NS_RANDOM = 2


def single_hash(b):
    return {'a': 0, 'c': 1, 'g': 2, 't': 3}[b]


def hash_kmer(kmer):
    h = 0
    for b in kmer:
        h <<= 2
        h |= single_hash(b)
    return h


def rehash(old_hash, b, k):
    l = 2 * k
    h = old_hash & ((1 << (l - 2)) - 1)
    h <<= 2
    h |= single_hash(b)
    return h


def compute_counts(seq, prev_counts=None, k=20):
    counts = defaultdict(int) if prev_counts is None else prev_counts
    h = hash_kmer(seq[:k])
    counts[h] += 1
    for b in seq[k:]:
        h = rehash(h, b, k)
        counts[h] += 1
    return counts


def preprocess(seq, nstrategy=NS_IGNORE):
    seq = seq.lower()
    if nstrategy == NS_IGNORE:
        seq = ''.join(filter(lambda x: x != 'n', seq))
    elif nstrategy == NS_SINGLE:
        seq = ''.join(('a' if b == 'n' else b) for b in seq)
    elif nstrategy == NS_RANDOM:
        seq = ''.join(random.choice(['a', 'c', 'g', 't']) if b == 'n' else b for b in seq)
    else:
        raise 'Invalid N strategy'
    return seq


def compute_histogram(counts):
    histogram = defaultdict(int)
    max_c = 0
    for c in counts.values():
        histogram[c] += 1
        max_c = max(max_c, c)

    return [histogram[i] for i in range(max_c + 1)]


def load_reads(fname):
    _, ext = path.splitext(fname)
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    with open(fname, "rU") as f:
        for read in SeqIO.parse(f, fmt):
            yield read.seq


def main(fname, out_fname, k, n_strategy):
    counts = defaultdict(int)
    for seq in load_reads(fname):
        seq = preprocess(seq, n_strategy)
        counts = compute_counts(seq, prev_counts=counts)
    hist = compute_histogram(counts)

    if out_fname:
        with open(out_fname, 'w') as f:
            for i, v in enumerate(hist):
                f.write('{} {}\n'.format(i, v))
    else:
        print(hist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute kmer_histogram form reads')
    parser.add_argument('-n', '--n-strategy',
                        default='IGNORE', help='Strategy to handle Ns')
    parser.add_argument('-k', type=int,
                        default=20, help='Kmer size')
    parser.add_argument('input', help='input file')
    parser.add_argument('-o', '--output', default=None, help='Output file')

    args = parser.parse_args()

    try:
        n_strategy_d = {
            'IGNORE': NS_IGNORE, 'SINGLE': NS_SINGLE, 'RANDOM': NS_RANDOM
        }
        n_strategy = n_strategy_d[args.n_strategy] if args.n_strategy in n_strategy_d\
            else int(args.n_strategy)
    except:
        n_strategy = NS_IGNORE

    main(args.input, args.output, args.k, n_strategy)
