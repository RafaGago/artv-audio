#!/usr/bin/env python3

from sympy import *

init_printing(use_latex="mathjax")

var("x a b c")

r1 = (-b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)
r2 = (-b - sqrt(b ** 2 - 4 * a * c)) / (2 * a)
f = a * (x - r1) * (x - r2)
pprint(simplify(f))
