#! /usr/bin/env python3
import argparse
import json
from pathlib import Path

from covest.data import count_reads_stats


def main(args):
    ar, rs = count_reads_stats(args.reads)
    print(rs, ar)
    if args.config:
        with open(str(args.config)) as f:
            config = json.load(f)
        config['reads_size'] = rs
        config['r'] = ar
        with open(str(args.config), 'w') as f:
            json.dump(config, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('reads', help='Input histogram')
    parser.add_argument('-c', '--config', type=Path, help='Add to config')

    args = parser.parse_args()
    main(args)
