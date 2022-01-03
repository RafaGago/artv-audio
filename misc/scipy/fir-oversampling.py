#!/usr/bin/env python

# https://scipy-cookbook.readthedocs.io/items/FIRFilter.html
import numpy
import scipy.signal
import pylab
import argparse
from sys import exit

parser = argparse.ArgumentParser(
    description = 'Calculate FIR oversampler coefficients.')

parser.add_argument(
    '-s', '--samplerate',
    type=float,
    help='sample rate',
    default=88200.)
parser.add_argument(
    '-c', '--cutoff-hz',
    type=float,
    help='Cutoff frequency (Hz)',
    default=22050.)
parser.add_argument(
    '-t', '--transition-width-hz',
    type=float,
    help='Transition width from passband to stopband (Hz)',
    default=2000.)
parser.add_argument(
    '-a', '--att-db',
    type=float,
    help='Stopband attenuation dB',
    default=110.)
parser.add_argument(
    '--forced-order',
    type=int,
    help='Force order on the filter. When set att_db will be ignored',
    default=None)
parser.add_argument(
    '--plot',
    help='Plot graphs',
    action='store_true',
    default=False)
parser.add_argument(
    '--minphase',
    help='Calculate low-effort min-phase converted version. The parameters won\'t match exactly.',
    action='store_true',
    default=False)
parser.add_argument(
    '--float',
    help='Print the coefficients as float.',
    action='store_true',
    default=False)

args = parser.parse_args()

samplerate = args.samplerate
cutoff_hz = args.cutoff_hz
att_db = args.att_db
transwidth_hz = args.transition_width_hz

norm_transwidth = transwidth_hz  / (samplerate * 0.5)
# TODO: For some unknown reason (or bug) this needs to be multiplied by two to
# match what this generator does:
# http://arc.id.au/FilterDesign.html
norm_transwidth *= 2

if not args.minphase:
    if (args.forced_order):
        att_db = scipy.signal.kaiser_atten (args.forced_order, norm_transwidth)

    order, kaiser_beta = scipy.signal.kaiserord (att_db, norm_transwidth)

    coeffs = scipy.signal.firwin(
        order,
        cutoff_hz,
        fs=samplerate,
        pass_zero='lowpass',
        scale=True,
        window=('kaiser', kaiser_beta))

    varname = 'linphase_fir_'
else:
    if (args.forced_order):
        att_db = scipy.signal.kaiser_atten(
            args.forced_order * 2, norm_transwidth * 2)
        att_db_factor = 1.
    else:
        att_db_factor = 2.

    order, kaiser_beta = scipy.signal.kaiserord(
        att_db * att_db_factor, norm_transwidth * 2)

    coeffs = scipy.signal.firwin(
        order * 2,
        cutoff_hz,
        fs=samplerate,
        pass_zero='lowpass',
        scale=True,
        window=('kaiser', kaiser_beta))

    print (len (coeffs))
    coeffs = scipy.signal.minimum_phase (list (coeffs))
    print (len (coeffs))
    varname = 'minphase_fir_'

coeffs[abs (coeffs) <= 1e-10] = 0.

varname += f'sr{int (samplerate)}hz_'
varname += f'fc{int (cutoff_hz)}hz_'
varname += f'tw{int (transwidth_hz)}hz_'
varname += f'att{int (att_db)}db'

vartype = 'float' if args.float else 'double'
coeff_substr = 'f' if args.float else ''

print (f'static std::array<{vartype}, {len (coeffs)}> {varname} = {{')
for coeff in coeffs:
    if coeff == 0.:
        print(f'  0.{coeff_substr},')
    else:
        print(f'  {coeff:.17g}{coeff_substr},')
print (f'}};')

print (f'cutoff (Hz): {cutoff_hz}')
print (f'transition width (Hz): {transwidth_hz}')
print (f'gain at {cutoff_hz + transwidth_hz}Hz (dB) {att_db}')
print (f'delay at destination sample rate: { (order - 1) / 2}')
print (f'order/num coeffs: {order}')
#assert (len (coeffs) == order)

if not args.plot:
    exit (0);
#------------------------------------------------
# Plot the magnitude response of the filter.
#------------------------------------------------
pylab.figure(1)
pylab.clf()
pylab.title ('Frequency Response')
pylab.xlabel ('Frequency (Hz)')
pylab.ylabel ('Gain (dB)')
pylab.ylim (-220., 3.)
pylab.grid (True)

w, h = scipy.signal.freqz (coeffs, worN=8000)
pylab.plot(
    (w / numpy.pi) * (samplerate * 0.5),
    20 * numpy.log10 (numpy.absolute (h)),
    linewidth=2)

#------------------------------------------------
# Plot the group delay of the filter.
#------------------------------------------------
pylab.figure(2)
pylab.title ('Group Delay At Oversampled Rate')
pylab.xlabel ('Frequency (Hz)')
pylab.ylabel ('Samples')
pylab.grid (True)

w, gd = scipy.signal.group_delay ((coeffs, 1))
pylab.plot(
    (w / numpy.pi) * (samplerate * 0.5), numpy.absolute (gd), linewidth=2)

pylab.show()
