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
assist in making the parts into classes. It is more powerful because it e.g.
allows placing all the coefficients of different unrelated classes together and
to run bulk coefficient smoothing.

The basic layout of a DSP part is:

```
struct part {
    // number of smoothable coefficients
    static constexpr uint n_coeffs;

    // number of internal/integer/unsmoothable coefficients usually 0.
    static constexpr uint n_coeffs_int;

    // number of states for a given type
    static constexpr uint n_states;

    // function (or overloads) for resetting coefficients.
    //
    // V is one of the clang native vector types, e.g. for a SSE double
    // double __attribute__ ((vector_size (16), __may_alias__))
    // for single precission just use vectors of one element. There are wrappers
    // on simd.hpp
    //
    // If this function belongs to a class where "n_coeffs_int" is 0, then the
    // second parameter, co_int, is omitted.

    template <class V>
    static void reset_coeffs(
        crange<V> co,      //coefficients, size >= "n_coeffs".
        crange<V> co_int, // coefficients, size >= "n_coeffs_int". Optional.
        ...               // other relevant parameters
        );

    // function to reset the states. Usually will just run memset to 0.

    template <class V>
    static void reset_states(
        crange<V> st, //states, size >= "n_states".
        ...           // other relevant parameters
        );

    // function to run the DSP process. Most of the time will process one sample
    // and rely on compiler inlining and good judgement when combining blocks at
    // a higher level. When it makes sense the function might process blocks.
    //
    // If this function belongs to a class where "n_coeffs_int" is 0, then the
    // second parameter, co_int, is omitted.

    template <class V>
    static X<V> tick(
        const crange<V> co,     // coeffs, size >= "n_coeffs".
        const crange<V> co_int, // coeffs, size >= "n_coeffs_int". Optional.
        crange<V>       st,     //states, size >= "n_states".
        X<V>            in,     // input data, e.g V, vec_complex<V>, etc.
        ...                     // other relevant parameters
        );

};
```
