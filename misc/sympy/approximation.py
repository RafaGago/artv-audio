#!/usr/bin/env python3

from sympy import *
from mpmath import mp

# https://mpmath.org/doc/current/calculus/approximation.html

init_printing(use_latex="mathjax")

var("v w")

mp.dps = 30 # more precision than double can handle...

def fx(x):
    return x / sqrt(x**2 + 1)

taylor = mp.taylor (fx, 0, 30)

i = 0;
for lv in taylor :
    print('{},'.format(lv))

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

def pan(x):
    return sin(x*1.5707963267948966)

print("Pan") # Not Remez but a one liner. At this point probably lolremez or
             # some python Remez implementation is also worth trying
chebyshev, err = mp.chebyfit(pan, [0, 1], 3, error=True)
pprint (Poly(chebyshev, v).as_expr())

print("Thiran") # Not Remez but a one liner
def thiran1(x):
    return 1 - (x + 0.418) / (1 + (x + 0.418))

thiransp = 1 - (v + 0.418) / (1 + (v + 0.418))

chebyshev, err = mp.chebyfit(thiran1, [0, 1], 5, error=True)
pprint (Poly(chebyshev, v).as_expr())


range = (v, 0, 1)
p1 = plot(thiransp, range, show=False, line_color='b')
p2 = plot (Poly(chebyshev, v).as_expr(), range, show=False, line_color='r')
p1.extend(p2)
p1.show()

p3 = plot (thiransp - Poly(chebyshev, v).as_expr(), range, show=False, line_color='r')
p3.show()
