import subprocess
import sys
from math import log, exp

from covest import constants
from .inverse import inverse


def print_wrap(x, label='', cond=True):
    if cond:
        print(label, x)
    return x


def verbose_print(message):
    if not constants.VERBOSE:
        return
    sys.stderr.write(message + "\n")


def safe_int(x):
    return int(x) if x != float('inf') else None


def fix_zero(x, val=1):
    if x == 0:
        return val
    else:
        return x


def safe_log(x):
    if x is None or x <= 0:
        return -constants.INF
    return log(x)


def estimate_p(cc, alpha):
    return (cc * (alpha - 1)) / (alpha * cc - alpha - cc)


def kmer_to_read_coverage(coverage, k, r):
    return coverage * r / (r - k + 1)


def fix_coverage(coverage):
    return inverse(lambda c: (c - c * exp(-c)) / (1 - exp(-c) - c * exp(-c)))(coverage)


def nonefloat(x):
    try:
        return float(x)
    except ValueError:
        return None


def run(command, shell=False, output=None, verbose=False):
    if verbose:
        print(command, file=sys.stderr)
    f = open(output, 'w') if output else None
    if not shell:
        command = command.split()
    return subprocess.call(command, shell=shell, stdout=f)
