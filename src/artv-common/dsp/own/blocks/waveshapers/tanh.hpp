#pragma once

#include <algorithm>
#include <cmath>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
// clang-format off
/*
https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/
(x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x) / (28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135)

Octave reminder:

x = sym('x') # symbolic
fn = x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x
fd = 28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135
f = fn / fd
ccode (horner (fn, x))
ccode (horner (fd, x))

At the end this was unused.
*/

// clang-format on
struct tanh_waveshaper_functions {
  //----------------------------------------------------------------------------
  template <class T>
  static T fn (T x)
  {
    tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int_fn (T x)
  {
    return x + log ((T) 1. + exp ((T) -2. * x));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int2_fn (T x)
  {
    // not solvable analitically in a practical/easy way.
    static_assert (false, "TBI?");
    return (T) 0.;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using tanh_waveshaper_adaa = adaa::waveshaper<tanh_waveshaper_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
