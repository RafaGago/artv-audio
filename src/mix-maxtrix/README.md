Mix maxtrix
===========

This plugin has 8 stereo buses in and out.

It has 8 mixer channels, any of them can take any combination (sum) of the 8
plugin input buses and output it to any combination its 8 plugin output buses.

Additionally each mixer channel can send its output to the neighboring mixer
channel inputs.

On top of that the 8 mixers can be divided on groups of 4 and 2 and placed in
series, as if two or four instance of the plugin were loaded on the host in
series.

Each mixer channel has:

- A FX slot, with multiple FX to select.
- A Knob section with the wet signal's mix amount, pan and MS conversion/
  Width.
- A Knob section controlling FX parameters.
- A Volume fader and level meter.
- A Knob section with the dry signal's pan and MS conversion/Width and an
  overall channel pan control.
- Mute/Solo.
- A utility section with phase inversion and some channel operations.

As mixer has an FX and it is possible to route them very flexibly, this plugin
can cover multiple roles. Depending on how it's used this plugin can be:

- A bus router, splitter, joiner.
- A MS/Converter.
- An in-track Mixer.
- A FX processor to place after multichannel Drum Machines/Synths.
- A crossover of maximum 8 bands.
- A multiband FX processor.
- A channel strip with easy to reorder modules.
- A creative multi-FX.
- An Open Source plugin bundles. So you don't need to install the plugins that
  are bundled separately.
- A C++ port of some popular Reaper Jesusonic FX.

History
========

This was started because of the limitations on Reaper when managing multiple
buses and FX on the same single track. I needed bus mixers, splitters and
joiners.

Reaper provided some on JSFX format, but while functional, I think that the
slider-based UI and sometimes the slider/knob ranges have room for improvement.

The initial idea of this FX was just to be an hybrid of mixer and routing matrix
or Reaper tracks, where any input stereo pair can be sent to any output pair and
vice-versa. This while providing a mixer with pan, volume and stereo width
knobs/sliders with ranges that make sense.

Then at some point I thought that maybe some simple bread-and-butter FX could be
added to declutter the amount of FX used on each Reaper track on different
buses, as I consider the mixer view unergonomic for this.

After adding some basic FX one thing led to another and I kept adding features
to the point I wasn't doing music anymore :)

Deliberate omissions
====================

- Preset browser. I don't use presets, this is a project done for fun, so I see
  doing this as a time consuming chore.

Included FX List
================

Compand/Expand
--------------

- Aw Consoles: AirWindows Console compand/expand mixers. To use them feed many
  inputs on a mixer with this FX selected. Al the inputs will be companded,
  summed and expanded.

Delay
-----

- BBD: "BBD Delay" by "Witti". Ported from JSFX.

- Delay Utility: Simple module to delay a signal a fixed amount of samples.

- Echo Cycles: By "Geraint Luff". Ported from JSFX.

- Spring Box: By "Geraint Luff". Ported from JSFX.

Distortion
----------

- Busscolors4: By "AirwWndows".

- Nonlinear: By "Lubomir I. Ivanov". Ported from JSFX.

- Sandwitch Amp: By "Geraint Luff". Ported from JSFX.

- Signal Crusher: By "Chokehold". Ported from JSFX.

Dynamics
--------

- 1175: By "SStillwell". Ported from JSFX.

- Event Horizon 2: By "SStillwell". Ported from JSFX.

- Consolidator: By "Chokehold". Ported from JSFX.

- GateExpander: By "Chokehold". Ported from JSFX.

- FairlyChildish: By "SStillwell". Ported from JSFX.

- MajorTom: By "SStillwell". Ported from JSFX.

- Slax: By "Sonic Anomaly". Ported from JSFX.

- Track Comp: By "Chokehold". Ported from JSFX.

- Transience: By "Saike". Ported from JSFX.

Exciter
--------

- 4x4: By "SStillwell". Ported from JSFX.

- BassProfessor: By "Sonic Anomaly". Ported from JSFX.

- HugeBooty: By "SStillwell". Ported from JSFX.

- SonicEnhancer: By "Lubomir I. Ivanov". Ported from JSFX.

EQ
---

- RBJ 1073: By "SStillwell". Ported from JSFX.

- Luftikus: By "LJKB".

- 4-band EQ: Simple uncramped parameteric EQ based on Adrew Simper's SVF.

- Filters x2: A module containing synth filters for using as a very colored

  EQ. The fiters are from Saike's Yutani synth, ported from JSFX.

Modulation
-----------

- Chow Phaser: By "Chow".

- Phaser: A phaser that can be used for mixing duties. For using as FX add some
  Dry signal.

- Ripple Phaser: By "Geraint Luff". Ported from JSFX.

- Ze Big Chorus: By "Smashed Transistors". Ported from JSFX.

- Ze Little Chorus: By "Smashed Transistors". Ported from JSFX.

Reverb
-------

- Atlantis Reverb: By "Lubomir I. Ivanov". Ported from JSFX.

- Dragonfly Early: By "Michael Willis".

- Dragonfly Hall: By "Michael Willis".

- Dragonfly Plate: By "Michael Willis".

- Dragonfly Room: By "Michael Willis".

- FDN Verb Riser: By "Shabtronic". Ported from JSFX.

- TAL Reverb 2: By TAL.

Stereo
-------

- Stereo Bub 3: By "Saike". Ported from JSFX.

- Stereo Tilt: By "Lubomir I. Ivanov". Ported from JSFX.
