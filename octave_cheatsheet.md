Octave cheatsheet
=================

As I never used Matlab before, I compile some useful commands here.

Symbolic math
-------------

Package import:

    package load symbolic

Define equations by:

    x = sym ('x')
    f = x / sqrt (x^2 1)

Find integrals:

    i_f = int (f)

To C code:

    ccode (f)

Plot:

    ezplot (f)

Limit:

    limit (f, x)

Horner, optimize number of operations:

    fn = x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x
    fd = 28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135
    f = fn / fd
    ccode (horner (fn, x))
    ccode (horner (fd, x))
