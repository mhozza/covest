# Binary Search
def inverse_bs(f, precision=1e-8):
    """Given a function y = f(x) that is a monotonically increasing function on
    non-negative numbers, return the function x = f_1(y) that is an approximate
    inverse, picking the closest value to the inverse, within delta."""
    def f_1(y):
        low, high = 0, float(y)
        last, mid = 0, high / 2
        while abs(mid - last) > precision:
            if f(mid) < y:
                low = mid
            else:
                high = mid
            last, mid = mid, (low + high) / 2
        return mid
    return f_1


# Newton's Method
def inverse(f, delta=1e-8):
    """Given a function y = f(x) that is a monotonically increasing function on
    non-negative numbers, return the function x = f_1(y) that is an approximate
    inverse, picking the closest value to the inverse, within delta."""
    def derivative(func):
        return lambda y: (func(y + delta) - func(y)) / delta

    def root(y):
        return lambda x: f(x) - y

    def newton(y, iters=15):
        guess = float(y) / 2
        rootfunc = root(y)
        derifunc = derivative(rootfunc)
        d = (rootfunc(guess) / derifunc(guess))
        while abs(d) > delta:
            guess -= d
            d = (rootfunc(guess) / derifunc(guess))

        return guess
    return newton
