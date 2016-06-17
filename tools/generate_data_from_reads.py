#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess

from Bio import SeqIO
from covest.constants import DEFAULT_K
from pathlib import Path

jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile} -o {infile}.jf'
jellyfish_hist = 'jellyfish histo {infile}.jf -o {outfile}'


def avg_read_length_and_size(file_path):
    ext = file_path.suffix
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    s = n = 0
    with open(str(file_path), "rU") as f:
        for read in SeqIO.parse(f, fmt):
            s += len(read.seq)
            n += 1
    return round(s / n), s


def run(command, shell=False, output=None):
    f = open(output, 'w') if output else None
    if not shell:
        command = command.split()
    return subprocess.call(command, shell=shell, stdout=f)


def main(args):
    dest = Path(args.dest).resolve()

    # copy reads file
    source = Path(args.source).resolve()
    seq_file = dest / ('reads%s' % source.suffix)
    if source.exists():
        if args.copy:
            if not seq_file.exists():
                print('Copying file...')
                shutil.copy2(str(source), str(seq_file))
            else:
                print('File already exists. Not overriding.')
        else:
            os.symlink(str(source), str(seq_file))
    else:
        print('File does not exist: %s' % source)
        exit(1)

    # calculate read_length and size
    if args.read_length and args.reads_size:
        read_len, reads_size = args.read_length, args.reads_size
    else:
        print('Calculating read length and size...')
        read_len, reads_size = avg_read_length_and_size(seq_file)
        print('read length: %d, size: %d' % (read_len, reads_size))

    # generate histogram
    print('Generating histogram...')
    hist_file = dest / ('%s.hist' % seq_file.stem)
    params = {
        'k': DEFAULT_K,
        'infile': seq_file,
        'outfile': hist_file,
    }
    run(jellyfish_count.format(**params), shell=True)
    run(jellyfish_hist.format(**params), shell=True)

    # generate_config
    config = {
        'reads': str(seq_file.name),
        'reads_size': reads_size,
        'hist': str(hist_file.name),
        'k': DEFAULT_K,
        'r': read_len,
    }
    print('Generating config...', config)
    with open(str(dest / 'config.json'), 'w') as f:
        json.dump(config, f)
    print('Done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dest', help='directory')  # required param from prepare_experiment_script
    parser.add_argument('source', help='fasta file with reads')
    parser.add_argument('read_length', nargs='?', default=None, help='average read length')
    parser.add_argument('reads_size', nargs='?', default=None, help='reads size')
    parser.add_argument('--copy', action='store_true', help='copy reads instead of linking')
    main(parser.parse_args())
