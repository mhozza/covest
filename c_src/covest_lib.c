#include <math.h>

#define MAX_EXP 200

long double truncated_poisson(long double l, int j) {
    if (l == 0 || l != l) {
        return 0;
    }
    long double p1 = 1;
    long double p3 = l;
    for (int i = 1; i<=j; i++) {
        p1 *= l / i;
    }
    while (l > MAX_EXP && p1 > 0) {
        p1 /= expl(MAX_EXP);
        l -= MAX_EXP;
    }
    if (l > 1e-8 && p1 > 0) {
        p3 = expl(l) - 1;
    }
    long double res = p1 / p3;
    return res;
}
