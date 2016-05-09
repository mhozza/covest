#! /usr/bin/env python3
import argparse
import json


def fix_genome_size(j):
    d = json.loads(j)
    d['estimated_genome_size'] = int(round(d['estimated_genome_size'] * 100))
    return json.dumps(d, sort_keys=True, indent=4, separators=(',', ': '))


def main(args):
    with open(args.fname, 'rU') as f:
        j = f.read()
    j = fix_genome_size(j)
    with open(args.fname, 'w') as f:
        f.write(j)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate reads form random genome with errors')
    parser.add_argument('fname')

    args = parser.parse_args()
    main(args)
