#!/usr/bin/env python3

from sympy import *
from sympy.plotting import plot

var("S k wc T z w f fc, argAP", real=True)

HAPnorm = (S**2 - k*S + 1) / (S**2 + k*S + 1)
# denormalization of the respose
HAP = HAPnorm.subs(S, S/wc)
# Function of jw
HAPjw = HAP.subs(S, I*w)
# Expand all frequencies
HAPfull = HAPjw.subs(w, 2*pi*f).subs(wc,2*pi*fc)

# https://en.wikipedia.org/wiki/Argument_(complex_analysis)
# Denominator is always bigger than zero so using atan(im/re)
HAP_phase_resp = atan(im(HAPfull)/re(HAPfull))
# phase response of the allpass

pprint(Eq(argAP, simplify(HAP_phase_resp)))

#HAP_BLT=HAP.subs(S, (2 * (z - 1))/( T * (z + 1)))
p1 = plot(
    HAP_phase_resp.subs(k, 0.001).subs(fc, 1000),
    (f, 1, 20000),
    x_var=f,
    xscale='log',
    show=False,
#    adaptive=False,
#    nb_of_points=20000
    )
p1.show()
