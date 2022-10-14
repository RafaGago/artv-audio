#!/usr/bin/env python3

from sympy import *
from mpmath import mp

# https://mpmath.org/doc/current/calculus/approximation.html

init_printing(use_latex="mathjax")

var("x a")

mp.dps = 30 # more precision than double can handle...

def f_x(x):
    return x / sqrt(x**2 + 1)

order=7

# build the correction polynomial
# https://www.kvraudio.com/forum/viewtopic.php?t=521377

taylor = mp.taylor (f_x, 0, order)
correction = 0

for idx, tblv in enumerate (taylor):
    term = x**idx * abs(tblv)
    apwr = idx - 2
    if apwr > 0:
        term *= a**apwr
    correction += term

correction = nsimplify(correction)

rsqrt_sigmoid = x / sqrt(x**2 + 1)

f = rsqrt_sigmoid.subs(x, correction)

fdx = f.diff(x)
fda = f.diff(a)

pprint(fdx)
pprint(fda)

stationary_points = solve(f.diff(x), x, dict=True)
print(stationary_points)

#fd = f.diff(x)
#fdd = fd.diff(x)
#
#roots = solveset(f, x)
