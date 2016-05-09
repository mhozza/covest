import argparse
import json
import subprocess

covest = 'covest {hist} -s {reads} -k {k} -r {r} -sp 16'


def run(command, output=None):
    f = open(output, 'w') if output else None
    return subprocess.call(command.split(), stdout=f)


def main(args):
    cfg_file = args.config_file
    with open(cfg_file) as f:
        cfg = json.load(f)
    run(covest % cfg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', default='config.json', help='config_file')
