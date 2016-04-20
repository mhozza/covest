import sys
from math import log

from covest import config


def print_wrap(x, label='', cond=True):
    if cond:
        print(label, x)
    return x


def verbose_print(message):
    if not config.VERBOSE:
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
        return -config.INF
    return log(x)
