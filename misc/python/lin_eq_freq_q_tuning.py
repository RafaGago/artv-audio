#!/usr/bin/env python3

import math


def get_pole(freq, q, sr):
    # bandpass
    w0 = 2.0 * math.pi * freq / sr
    cosw0 = math.cos(w0)
    alpha = math.sin(w0) / (2.0 * q)

    a0 = 1 + alpha
    norm = 1 / a0

    a1 = (-2 * cosw0) * norm
    a2 = (1 - alpha) * norm

    a = 1
    b = a1
    c = a2

    inv_2a = 1.0 / (2.0 * a)
    re = -b * inv_2a
    im = (b * b) - 4.0 * a * c

    if im < 0:
        # complex conjugate
        im = math.sqrt(-im) * inv_2a
    else:
        # real poles
        # returning the one with the biggest magnitude
        re = math.abs(re) + math.sqrt(im) * inv_2a
        im = 0
    return complex(re, im)


def get_stages(pole, snr_db):
    pole_mod = abs(pole)
    return math.log2(-snr_db / (20.0 * math.log10(pole_mod)))


sr = 48000
f = 800
q = 5
db = 90

print(get_stages(get_pole(f, q, sr), db))
