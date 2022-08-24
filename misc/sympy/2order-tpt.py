#!/usr/bin/env python3

from sympy import *

init_printing(use_latex="mathjax")

var("yHP yBP yLP yLP2 yAP x s1 s2 R g k")

# k = 2 * R

# ART of VA filter design. 4.14

yHP = (x - (k + g)*s1 - s2) / (1 + k*g + g*g)
#y_bp = Eq (yBP, g * (x - (k*yBP + g*yBP + s2)) + s1)
#y_bp = solve(y_bp, yBP)
yBP = (g * (x - s2) + s1) / (1 + k*g + g*g)
# From Cytomic SVF
yLP = s2 + (g*(s1+g*(-s2+x)))/(g*g + k*g + 1)
yAP = 1 - 2*k*yBP

print("HP")
pprint(simplify(yHP))
print("\nBP")
pprint(yBP)
print("\nLP")
pprint(simplify(yLP))
print("\nAP")
pprint(simplify(yAP))
print("\nAP expanded")
pprint(expand(yAP))


##sqrt sig
##unsolvable
#var("x k S G u", real=True)
#fback_sigmoid = Eq ((x-k*S)-u, k*G*(u/sqrt(u*u+1)))
#pprint(solve(fback_sigmoid, u))

# tanh vaneev
# sympy doesn't like support abs(?)
var("x k S G u", real=True)
vaneev_u = (x  * ( 2.45550750702956 +  2.45550750702956 * abs(u) +
    ( 0.893229853513558 +  0.821226666969744 * abs(u)) * u*u)
    / ( 2.44506634652299 + ( 2.44506634652299 + u*u) *
    abs (u +  0.814642734961073 * u * abs(u))))
fback_sigmoid = Eq ((x-k*S)-u, k*G*vaneev_u)
pprint(solve(fback_sigmoid, u))
