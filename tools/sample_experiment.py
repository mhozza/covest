#!/usr/bin/env python3
import argparse


def main(args):
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', help='experiment_directory')
    parser.add_argument('dest', help='destination directory')
    parser.add_argument('-c', '--coverage', help='target coverage')
    parser.add_argument('-t', '--times', help='divide current data coverage by this value')
    main(parser.parse_args())
