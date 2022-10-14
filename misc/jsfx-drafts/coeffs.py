#!/usr/bin/env python
from math import *
import sympy as sp

def f(x):
  return 1 / (x + 1 - log(2))

def get_start_times(del_spls, n_ap):
  print (f'{del_spls}[{sp.prevprime (del_spls)},{sp.nextprime (del_spls)}]')
  for i in range (n_ap):
    spls = del_spls * f (i + 1)
    print (f'{spls}[{sp.prevprime (spls)},{sp.nextprime (spls)}]')

get_start_times(3000, 6)
