#!/usr/bin/env python3

from sympy import *
from sympy.plotting import plot

#init_printing(use_latex="mathjax")

var("yHP yBP yLP yLP2 yAP x s1 s2 R g k A")

# k = 2 * R

# ART of VA filter design. 4.14
yHP = (x - (k + g)*s1 - s2) / (1 + k*g + g*g)
yBP = (g * (x - s2) + s1) / (1 + k*g + g*g)
yLP = g * yBP + s2
yBPnorm = yBP * k
yAP = x - 2 * yBPnorm
#yBell = x + yBP * 2 * k * A
yNotch = x - yBPnorm
yPeak = yLP -yHP
####

yLP1 = g * (x + s1) / (1 + g)
yHP1 = x - yLP1
yAP1 = yLP1 - yHP1
yLS1 = x + k * yLP1
yHS1 = x + k * yHP1

def print_fun(f, name):
    pprint(Eq(S('{}_y'.format(name)),simplify(f)))
    f_GS = div (f, x)
    pprint(Eq(S('{}_G'.format(name)),simplify(f_GS[0])))
    pprint(Eq(S('{}_S'.format(name)),simplify(f_GS[1])))

print_fun (yHP, "HP")
print_fun (yLP, "LP")
print_fun (yBPnorm, "BP")
print_fun (yPeak, "Peak")
print_fun (yBP, "BP_raw")
print_fun (yAP, "AP")
print_fun (yNotch, "Notch")

print_fun (yLP1, "LP1")
print_fun (yHP1, "HP1")
print_fun (yAP1, "AP1")
print_fun (yLS1, "LS1")
print_fun (yHS1, "HS1")
##sqrt sig
##unsolvable
#var("x k S G u", real=True)
#fback_sigmoid = Eq ((x-k*S)-u, k*G*(u/sqrt(u*u+1)))
#pprint(solve(fback_sigmoid, u))

## tanh vaneev
## sympy doesn't like support abs(?)
#var("x k S G u", real=True)
#vaneev_u = (x  * ( 2.45550750702956 +  2.45550750702956 * abs(u) +
#    ( 0.893229853513558 +  0.821226666969744 * abs(u)) * u*u)
#    / ( 2.44506634652299 + ( 2.44506634652299 + u*u) *
#    abs (u +  0.814642734961073 * u * abs(u))))
#fback_sigmoid = Eq ((x-k*S)-u, k*G*vaneev_u)
#pprint(solve(fback_sigmoid, u))
#
## parabola (axis reversed)
#var ("xv yv y1")
#
#yparab = solve (Eq(0, (2 * y1 - 1) / (y1**2) * y**2 - x * y + x - y), x)
#
#print (yparab[0])
#pprint (simplify (yparab[0]))
#pprint (simplify (yparab[1]))
##p1 = plot(yparab, show=False)
##p1.show()
