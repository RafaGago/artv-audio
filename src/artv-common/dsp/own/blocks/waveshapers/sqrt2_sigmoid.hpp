#pragma once

#include <cmath>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// cascading the sqrt sigmoid
//
// code as with the "ccode" function
struct sqrt2_sigmoid_functions {

  template <class T>
  static T fn (T x)
  {
    // "x / sqrt (2 * x^2 + 1)" is the result of cascading
    // "x / sqrt (x^2 + 1)"
    // its positive limit is "sqrt(2)/2", so it is adjusted to -1 / 1.
    return ((T) M_SQRT2 * x) / sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int_fn (T x)
  {
    return (T) (1.0 / 2.0) * (T) M_SQRT2 * sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int2_fn (T x)
  {
    return (1.0 / 2.0) * (T) M_SQRT2
      * ((T) (1.0 / 2.0) * x * sqrt ((T) 2. * x * x + 1.)
         + (T) (1.0 / 4.0) * M_SQRT2 * asinh (M_SQRT2 * x));
  }
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt2_waveshaper_adaa = adaa::waveshaper<sqrt2_sigmoid_functions, order>;

//------------------------------------------------------------------------------
} // namespace artv
