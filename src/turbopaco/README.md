Turbo Paco
==========

Old-style Algorithmic Reverb, based on some own and other classic algorithms
reverse engineered and posted by the user "Acreil" on this forum thread:

https://gearspace.com/board/geekzone/380233-reverb-subculture.html

This Effect is mostly about Allpass loop and Feedforward comb reverbs with few
exceptions. There are even some Choruses.

The controls are pretty much self-explanatory except Dyn Time. When Dyn Time is
positive the dynamics processor does ducking. When it's negative it does inverse
ducking, which results in a kind of gate.

The design philosophy is to have a limited set of controls and many algorithms.
Now there aren't a lot, but more can be added.

The algorithms are implemented on:

* 16-bit fixed point delay lines with noise shaping with a 32-bit accumulator for
  computation.
* 16-bit custom-format floating point encoding delay lines with single precision
  floating point computation.
* Single precision floating point.

The rate of the the sample rate converter can be changed between some preset
frequencies, some "overclocked" and others "underclocked".

The reverb is block processed and uses low samplerates, so some of the
algorithms are very CPU friendly.

My own algorithms are rookie stuff. Doing reverb is complex. I tried my best for
them to sound acceptable. On the algorithm dropdown the algorighms prefixed with
"AC" come from the user "Acreil" on gearspace. Those with no prefix are mine.

Thanks to all fellas on KVR's DSP forum for the guidance :)

https://www.kvraudio.com/forum/viewforum.php?f=33

Enjoy!!!
