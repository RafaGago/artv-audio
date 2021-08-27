This folder contains building blocks for creating DSP objects. Most of all,
barring some exceptions, is organized on classes that:

- Only have static methods.
- Don't manage the memory they operate with.
- Process in a sample basis.
- Have the same function names and members, an unwritten interface.

The idea is to:

-To have full control of the memory layout when implementing some DSP using
 various blocks.

-Consequence of the above, to allow bulk cross-DSP/object coefficient smoothing.

The obvious caveat is that there is memory to manage.

As all static objects here follow the same layout, it could be easy to convert
them to stateful/OO with a class that takes the block to convert as a template
parameter, as every class has information of the number of coefficients and
states it has.
