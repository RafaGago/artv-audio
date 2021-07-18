#!/usr/bin/env python

# This requires pysr.

import math
import numpy as np
from pysr import pysr, best, best_row

def linterp (x0, x1, y0, y1, x):
    return (y0 * (x1 -x) + y1 * (x - x0)) / (x1 - x0)

x_max = 40
x_max_weight_lim = 1
w_min = 0.3
# If the function is symmetric this will be handled when implementing on C by.
# running fabs on the input and capturing the sign. Saving the program the
# trouble.
f_is_symmetric = True

unary_operators=[
    "x2(x) = x^2",
    "x3(x) = x^3",
    "x4(x) = x^4",
    "x5(x) = x^5"]

extra_sympy_mappings={
    "x2": lambda x : x ** 2,
    "x3": lambda x : x ** 3,
    "x4": lambda x : x ** 4,
    "x5": lambda x : x ** 5}

if not f_is_symmetric:
    unary_operators.append ("abs")

# only positive ranges
x = np.geomspace (0.0000000000000001, x_max, 50)
# weighting trapezoidal window
w = [1 if v <= 1. else linterp (x_max_weight_lim, x_max, 1, math.sqrt(w_min), v)
     for v in x]
# weighting quadratic window (assumes normalized values)
w = [v * v for v in w]
w = np.array (w)

if f_is_symmetric:
    # add zero
    x = np.concatenate ((np.array ([0]), x), axis=0)
    w = np.concatenate ((np.array ([0]), w), axis=0)
else:
    # add negative range
    x = np.concatenate ((np.flip (x) * -1, x), axis=0)
    w = np.concatenate ((np.flip (w), w), axis=0)

y = np.tanh(x)

# Learn equations
equations = pysr(
    x,
    y,
    weights=w,
    #niterations=2,
    variable_names=["x"],
    binary_operators=["+", "*", "-", "/"],
    unary_operators=unary_operators,
    extra_sympy_mappings=extra_sympy_mappings
    )

# (you can use ctl-c to exit early)
print(best(equations))
