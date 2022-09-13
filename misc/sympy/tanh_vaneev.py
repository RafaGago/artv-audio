#!/usr/bin/env python3

from sympy import *

init_printing(use_latex="mathjax")

var("x h")

num = x * (2.45550750702956 + 2.45550750702956 * abs(x) + (0.893229853513558 + 0.821226666969744 * abs (x)) * x * x)
den = 2.44506634652299 + (2.44506634652299 + x * x) * abs (x + 0.814642734961073 * x * abs (x))

tvaneev = num / den

pprint (limit(tvaneev, x, oo))

tvaneev = num / (den * h)
tvaneev = tvaneev.subs(x, x * h)

pprint (limit(tvaneev, x, oo))
# print('tanh_vaneev with drive lim: {:.35f}'.format(limit(tvaneev, x, oo)))
