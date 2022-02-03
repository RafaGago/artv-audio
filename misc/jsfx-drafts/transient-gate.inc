desc:Transient gate

slider1:recovery=        0.1<0, 1, 0.01>              Recovery
slider2:curve_decay=     0.1<0, 1, 0.01>              Curve Decay
slider3:hf_freq=         60.<10, 200, 0.1>            Detector HP(hz)
slider4:output=          0.<0, 1, 1{transient,tonal}> Output

in_pin:left input
in_pin:right input
out_pin:left output
out_pin:right output

import _core.jsfx-inc
import _dynamics.jsfx-inc
import _filters.jsfx-inc

@init

@slider

tgate.transient_gate_reset(
  DEFAULT_TGATE_K, recovery, curve_decay, hf_freq, 200, srate
  );

@sample

gate = tgate.transient_gate_tick (spl0, spl1);
gate = output == 1 ? (1 - gate) : gate;

spl0 *= gate;
spl1 *= gate;
