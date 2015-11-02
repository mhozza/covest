import sys
import config


def print_wrap(x, label='', cond=True):
    if cond:
        print(label, x)
    return x


def verbose_print(message):
    if not config.VERBOSE:
        return
    sys.stderr.write(message + "\n")
