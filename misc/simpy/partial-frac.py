#!/usr/bin/env python3

from sympy import *

init_printing(use_latex="mathjax")


# 4 pole
var("p1_re p1_im p2_re p2_im z1 h")
p1 = p1_re + I * p1_im
p1c = p1_re - I * p1_im
p2 = p2_re + I * p2_im
p2c = p2_re - I * p2_im
h_fn = Eq(h, 1 / ((1 - p1 * z1) * (1 - p1c * z1) * (1 - p2 * z1) * (1 - p2c * z1)))
h_fn_part = apart(h_fn, z1)


# 4 pole, scipy dies

# 4 pole
# var("p1_re p1_im p2_re p2_im z h")
# z1 = pow(z, -1)
# p1 = p1_re + I * p1_im
# p2 = p2_re + I * p2_im
# h_fn = Eq(
#    h,
#    1
#    / (
#        (1 - p1 * z1)
#        * (1 - conjugate(p1) * z1)
#        * (1 - p2 * z1)
#        * (1 - conjugate(p2) * z1)
#    ),
# )
# h_fn_part = apart(h_fn, z)

# 2 pole
# var("p1_re p1_im z h")
# z1 = pow(z, -1)
# p1 = p1_re + I * p1_im
# p1c = p1_re - I * p1_im
# h_fn = Eq(h, 1 / ((1 - p1 * z1) * (1 - p1c * z1)))
# h_fn_part = apart(h_fn, z)
