#!/usr/bin/env python3

from sympy import *

init_printing(use_latex="mathjax")

# var("h s a q")
## Peaking EQ H(s) RBJ
# num = s ** 2 + s * (a / q) + 1
# den = s ** 2 + s / (a * q) + 1
#
# zeros = roots(num, s, multiple=True)
# zeros[0] = simplify(zeros[0])
# zeros[1] = simplify(zeros[1])
# poles = roots(den, s, multiple=True)
# poles[0] = simplify(poles[0])
# poles[1] = simplify(poles[1])
#
# print("zeros")
# print(zeros)
# print("poles")
# print(poles)
# print("zeros verification")
# print(simplify((s + zeros[0]) * (s + zeros[1])))
# print("poles verification")
# print(simplify((s - poles[0]) * (s - poles[1])))

# From the RBJ peaking filter biquad.

var("h z a0 a1 a2 b0 b1 b2 alpha w0 A Q")

num = b0 + (b1 * z) + (b2 * z ** 2)
den = a0 + (a1 * z) + (a2 * z ** 2)

pprint(num)
pprint(den)

zeros = roots(num, z, multiple=True)
poles = roots(den, z, multiple=True)

zeros = roots(num, z, multiple=True)
poles = roots(den, z, multiple=True)

# Let the program write the quadratic equation formulas by using biquad coeffs
# name
print("poles")
pprint(poles)
print("zeros")
pprint(zeros)
# The takeaway from this verification is that b2 and a2 are the gain.
print("zeros verification")
pprint(simplify(b2 * (z - zeros[0]) * (z - zeros[1])))
print("poles verification")
pprint(simplify(a2 * (z - poles[0]) * (z - poles[1])))

# substitute
for pole in poles:
    print("pole")
    f = pole.subs(
        a0,
    )
    f = f.subs(a1, -2 * cos(w0))
    f = f.subs(a2, 1 - alpha / A)
    f = f.subs(alpha, sin(w0) / (2 * Q))
    f = simplify(f)
    pprint(f)
    print("pole non-pretty")
    print(f)

print("poles gain")
print()
print("poles gain non-pretty")


for zero in zeros:
    print("zero")
    f = zero.subs(b0, 1 + alpha * A)
    f = f.subs(b1, -2 * cos(w0))
    f = f.subs(b2, 1 - alpha * A)
    f = f.subs(alpha, sin(w0) / (2 * Q))
    f = simplify(f)
    print("zero non-pretty")
    pprint(simplify(f))
