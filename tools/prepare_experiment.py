#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess
import sys

from pathlib import Path


def main(args):
    # create dir
    p = Path(args.name)
    if not args.force and p.exists():
        print('Path already exists. Use --force to use it.', file=sys.stderr)
        exit(1)
    try:
        p.mkdir()
    except FileExistsError:
        pass

    # gather data
    if args.datasource:
        ds = Path(args.datasource)
        if ds.exists() and ds.is_dir():
            for f in ds.iterdir():
                shutil.copy2(str(f), str(p))
        elif ds.exists() and not os.access(str(ds), os.X_OK):
            shutil.copy2(str(ds), str(p))
        else:
            print('Executing "%s":' % ds, file=sys.stderr)
            command = [str(ds.resolve()), args.name] + args.params
            if args.copy:
                command.append('--copy')
            ret = subprocess.call(command)
            if ret:
                print('"%s" failed with return code: %d.' % (ds, ret), file=sys.stderr)
                exit(2)
            else:
                print('"%s" succeeded.' % (ds,), file=sys.stderr)

    # copy runner script
    if args.runner:
        r = Path(args.runner).resolve()
        if r.exists():
            if args.copy:
                shutil.copy2(str(r), str(p / 'run'))
            else:
                os.symlink(str(r), str(p / 'run'))
        else:
            print('File does not exist: %s' % r, file=sys.stderr)
            exit(2)

    if args.config:
        cfg = p / 'config.json'
        cfg_json = json.loads(args.config)
        json.dump(cfg_json, cfg.open(mode='w'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare experiment')
    parser.add_argument('name', help='directory')
    parser.add_argument('-t', '--type', type=str, choices=['simulated', 'real'],
                        default='simulated', help='Type of experiment')
    parser.add_argument('-d', '--datasource', type=str,
                        help='Data source: either generation script or directory'
                             ' or kmer histogram file')
    parser.add_argument('-r', '--runner', type=str,
                        help='Experiment runner: script for running the experiment')
    parser.add_argument('-c', '--config', type=str,
                        help='parameters for the experiment runner')
    parser.add_argument('-f', '--force', action='store_true',
                        help='force use existing path')
    parser.add_argument('-p', '--params', nargs='*', default=tuple(),
                        help='generation script params')
    parser.add_argument('--copy', action='store_true', help='copy reads instead of linking')

    main(parser.parse_args())
