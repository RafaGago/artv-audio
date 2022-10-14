#!/usr/bin/env python

# This seems to be single threaded, A C++ version should be better.

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from geneticalgorithm import geneticalgorithm as ga

class Delay:
  def __init__(self, samples):
    self.spls = []
    self.spls.extend([0.0] * int(samples))
    self.idx = 0

  def pop (self):
    size = len (self.spls)
    idx = self.idx + size
    if idx >= size:
      idx -= size
    return self.spls[idx]

  def push (self, spl):
    self.spls[self.idx] = spl
    self.idx += 1
    if self.idx == len(self.spls):
      self.idx = 0

def allpass (x, yn, g):
  u = x + yn * g
  return (yn - u * g, u)

t_sec = 10
fs = 33600
N = t_sec * fs

# initial noise burst
noise_power = 0.001 * fs / 2
noise = np.random.normal (scale = np.sqrt (noise_power), size=fs)

# run an iteration of allpasses
def run_iter (n_spls):
  delays = []
  for spls in n_spls:
    delays.append (Delay (spls))

  # main loop
  feedback = 0

  outs = np.zeros (N) # 1 dimension, size N
  for i in range (0, N):
    x = 0
    if i < noise.size:
        x = noise[i]
    x += feedback

    for j in range (0, len (delays)):
        yn = delays[j].pop()
        ret = allpass (x, yn, 0.618)
        x = ret[0]
        delays[j].push(ret[1])

    outs[i] = x
    feedback = x

  return outs

# run for some seconds and check the standard deviation of the spectrum
def optim_fun (args):
  res = run_iter (args.tolist())
  _, y = signal.welch (res, fs, nperseg=32*1024)
  return np.std (y) # TODO: autocorrelation

def do_plot(outs):
  f1, y1 = signal.welch (outs, fs, nperseg=32*1024)
  f2, y2 = signal.welch (noise, fs, nperseg=32*1024)

  print (f'Std: {np.std(y1)}, Noise Std {np.std(y2)}')

  plt.figure()
  plt.subplot (211)

  plt.semilogy (f1, y1)
  plt.ylim ([0.1, 1])
  plt.xlabel ('frequency [Hz]')
  plt.ylabel ('AP PSD [V**2/Hz]')

  plt.subplot (212)
  plt.semilogy (f2, y2)
  plt.ylim ([0.001, 1])
  plt.xlabel ('frequency [Hz]')
  plt.ylabel ('Noise PSD [V**2/Hz]')

  plt.show()

def optimize (n_allpasses):
  model = ga(
    function=optim_fun,
    dimension=n_allpasses,
    variable_type_mixed=np.array ([['real']] * n_allpasses),
    variable_boundaries=np.array ([[140, 2600]] * n_allpasses)
    )
  return model

model = optimize (6)
model.run()
