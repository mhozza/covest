import time
import sys
from contextlib import contextmanager

stack = list()
messages = list()


def indent():
    if len(stack) < 2:
        return ''
    return ('| ' * (len(stack) - 2)) + '+'


def push(cnt=1):
    for _ in range(cnt):
        stack.append(time.time())


def pop(cnt=1):
    ret = None
    for _ in range(cnt):
        ret = stack.pop()
    return ret


def replace():
    pop()
    push()


def get_time(back=0):
    slen = len(stack)
    if back >= slen or back == -1:
        back = slen - 1
    return stack[slen - back - 1]


def print_all():
    global messages
    for message in messages:
        sys.stderr.write(message + "\n")
    messages = []


def msg(message, back=0):
    timediff = time.time() - get_time(back)
    messages.append(indent() + message.format(time=timediff))
    print_all()


def running_time_decorator(fn):
    def wrapped(*p, **k):
        push(2)
        ret = fn(*p, **k)
        msg("Function " + fn.__name__ + " (" + fn.__module__ +
            ") took {time} seconds.", 1)
        pop(2)
        return ret
    return wrapped


@contextmanager
def running_time(text):
    push(2)
    yield
    msg(text + " took {time} seconds.", 1)
    pop(2)
