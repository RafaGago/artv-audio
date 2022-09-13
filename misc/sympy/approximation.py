#!/usr/bin/env python3

from sympy import *
import mpmath as mp

# https://mpmath.org/doc/current/calculus/approximation.html

init_printing(use_latex="mathjax")

var("v")

def fx(x):
    return sin (x)

taylor = mp.taylor (fx, 0, 12)

print("Taylor")
# Note list of coeffs reversal
pprint (Poly(taylor[::-1], v).as_expr())

num, den = mp.pade (taylor, 6, 6)

print("Pade")
# Note list of coeffs reversal
pprint (Poly(num[::-1], v) / Poly(den[::-1], v))

print("Chebishev") # Not Remez but a one liner
chebyshev, err = mp.chebyfit(fx, [0, 1.6], 6, error=True)
# Note NO list of coeffs reversal
pprint (Poly(chebyshev, v).as_expr())
pprint (err)
