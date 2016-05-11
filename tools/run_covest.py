#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
from pathlib import Path


covest = 'covest {hist} -s {reads} -k {k} -r {r} -sp 16 -m repeat > output.yml'
wd = Path(__file__).parent


def run(command, output=None):
    f = open(output, 'w') if output else None
    return subprocess.call(command.split(), stdout=f)


def main(args):
    os.chdir(wd.path)
    cfg_file = args.config_file
    with open(cfg_file) as f:
        cfg = json.load(f)
    print(cfg)
    run(covest.format(**cfg))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', nargs='?', default='config.json', help='config_file')
    main(parser.parse_args())