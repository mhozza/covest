import argparse
import json
import shutil
import subprocess

from Bio import SeqIO
from covest.constants import DEFAULT_K
from pathlib import Path

jellyfish_count = 'jellyfish count -m {k} -s 500M -t 16 -C {infile} -o {infile}.jf'
jellyfish_hist = 'jellyfish histo {infile}.jf -o {infile_base}_k{k}.hist'


def avg_read_length(file_path):
    ext = file_path.suffix
    fmt = 'fasta'
    if ext == '.fq' or ext == '.fastq':
        fmt = 'fastq'
    s = n = 0
    with open(str(file_path), "rU") as f:
        for read in SeqIO.parse(f, fmt):
            s += len(read.seq)
            n += 1
    return round(s / n)


def run(command, output=None):
    f = open(output, 'w') if output else None
    return subprocess.call(command.split(), stdout=f)


def main(args):
    dest = Path(args.dest)

    # copy reads file
    source = Path(args.source)
    seq_file = dest / 'reads%s' % source.suffix
    if source.exists():
        shutil.copy2(str(source), str(seq_file))

    # calculate read_length
    read_len = avg_read_length(seq_file)

    # generate histogram
    params = {
        'k': DEFAULT_K,
        'infile': seq_file
    }
    run(jellyfish_count.format(**params))
    run(jellyfish_hist.format(**params))

    # generate_config
    config = {
        'reads': str(seq_file),
        'hist': '%s.hist' % seq_file.stem,
        'k': DEFAULT_K,
        'r': read_len,
    }
    with open(str(dest / 'config.json'), 'w') as f:
        json.dump(f, config)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dest', help='directory')  # required param from prepare_experiment_script
    parser.add_argument('source', help='fasta file with reads')
