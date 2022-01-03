#!/usr/bin/env python

# fit a straight line to the economic data
import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def objective (x, a, b, c, d, e, f):
    return a * x**7 + b * x**5 + c * x**3 + d * x**2 + e * x + f

def linterp (x0, x1, y0, y1, x):
    return (y0 * (x1 -x) + y1 * (x - x0)) / (x1 - x0)

x_max = 2
x_max_weight_lim = 2
w_min = 0.3
f_is_symmetric = True

# only positive ranges
x = np.geomspace (0.00000001, x_max, 100)
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

popt, _ = curve_fit (objective, x, y)
# summarize the parameter values
a, b, c, d, e, f = popt

print (f'{a} * x**5 + {b} * x**4 + {c} * x**3 + {d} * x**2 + {e} * x + {f}')
print(y)

plotx = np.linspace (0., x_max, 20000)

plt.plot (plotx, np.tanh (plotx))
plt.plot (plotx, objective (plotx, a, b, c, d, e, f))
plt.show()
plt.show()
