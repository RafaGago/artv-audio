#!/usr/bin/env python

# This requires pysr.

import math
import numpy as np
from pysr import pysr, best, best_row

def linterp (x0, x1, y0, y1, x):
    return (y0 * (x1 -x) + y1 * (x - x0)) / (x1 - x0)

x_max = 40
x_max_weight_lim = 2
w_min = 0.3
# If the function is symmetric this will be handled when implementing on C by.
# running fabs on the input and capturing the sign. Saving the program the
# trouble.
f_is_symmetric = True

unary_operators=[
    'pow2(x) = x^2',
    'sqrt',
    #'sigm(x) = x / (x+1)',
    #'pow3(x) = x^3',
    #'pow4(x) = x^4',
    #'pow5(x) = x^5',
    ]

extra_sympy_mappings={
    'pow2': lambda x : x ** 2,
    #'sigm(x)': lambda x : x / (x + 1),
    #'pow3': lambda x : x ** 3,
    #'pow4': lambda x : x ** 4,
    #'pow5': lambda x : x ** 5,
    }

if not f_is_symmetric:
    unary_operators.append ('abs')

# only positive ranges
x = np.geomspace (0.00000001, x_max, 60)
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

y = x - np.log(np.cosh(x))
sigm = x / (x + 1) # adding a feature with a raedy to go sigmoid

equations = pysr(
    #np.transpose (np.array ([x, sigm], ndmin=2)),
    x,
    y,
    weights=w,
    niterations=1000, # ctrl + c
    variable_names=['x'],
    binary_operators=['+', '*', '-', '/'],
    unary_operators=unary_operators,
    procs=10,
    populations=16*10,
    #update=False,
    multithreading=True,
    maxsize=35,
    #select_k_features=4,
    parsimony=1e-4, # how complex it becomes so fast
    shouldOptimizeConstants=True,
    #warmupMaxsizeBy=0.02,
    extra_sympy_mappings=extra_sympy_mappings,
    useFrequency=False
    )

# (you can use ctl-c to exit early)
print(best(equations))

# more or less OK tanh:
#
# sigm(x) = x / (x+1)
#
# 12|8.0615973e-7|sigm(x + pow2((pow2(x) + 2.9968638) * (x * 0.3682524)))
# 19|3.3242973e-8|sigm(pow2(pow2(x * ((((1.2004293 - x) - pow2(x)) * 0.04444416) + -0.6549788)) + # x) + x)
# 21|2.1726796e-9|sigm(pow2(pow2(x * (((1.3981539 - pow2((x - -0.043486368) + 1.1953219)) * 0.# 03136203) + -0.5794995)) + x) + x)
# 34|9.061694e-10|sigm(pow2(pow2(x * ((((0.6233191 - pow2((x - -0.044878364) + 0.9664551)) + sigm# (pow2((((0.90934205 - x) + -0.07033085) - pow2(x)) * 0.20001386))) * 0.037465252) + -0.5673225)) # + x) + x)
# 35|4.664905e-10|sigm(pow2(pow2(x * ((((1.2158864 - pow2((x - -0.03984681) + 0.67147267)) + sigm# (pow2(((pow2(2.306366 - x) + -0.5386064) - pow2(x)) * -0.14637995))) * 0.04096665) + -0.62395185)# ) + x) + x)

# Learn equations
# (pow2(pow2(x * (((1.3981539 - pow2((x - -0.043486368) + 1.1953219)) * 0.03136203) + -0.5794995)) + x) + x)
#
# bool pos = x >= 0.;
# x = abs (x);
# auto v = (x - -0.043486368) + 1.1953219)
# v *= v;
# v = 1.3981539 - v * 0.03136203;
# v = x * v + -0.5794995;
# v = v * v + x;
# v = v * v + x;
# return (pos ? v : -v) / (v + 1);
