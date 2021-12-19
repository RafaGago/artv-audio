This folder contains building blocks for creating DSP objects. Most of all,
barring some exceptions, is organized on classes that:

- Only have static methods.
- Don't manage the memory they operate with.
- Process in a sample basis.
- Have the same function names and members: an unwritten interface.

The idea is to:

-To have full control of the memory layout when implementing some DSP using
 various blocks. E.g. when having multiple DSP parts together in one effect or
 module, it might be interesting to store all coefficients that aren't modified
 together on the same cache line and to have all the state together afterwards.
 By using classes reasoning as this level is impossible.

-To allow bulk smoothing of the coefficients of multiple different DSP processes
 externally by using one-pole smoothers.

The obvious caveat is that there is memory to manually manage.

As all static objects here follow similar layouts, there is a helper class to
assist in making the parts into classes.
