#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// cascading the sqrt sigmoid
//
// code as with the "ccode" function
struct sqrt2_sigmoid_functions {

  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    // "x / sqrt (2 * x^2 + 1)" is the result of cascading
    // "x / sqrt (x^2 + 1)"
    // its positive limit is "sqrt(2)/2", so it is adjusted to -1 / 1.
    return ((T) M_SQRT2 * x) / sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    V num = (T) M_SQRT2 * x;
    V den = vec_sqrt ((T) 2. * x * x + (T) 1.);
    return num / den;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return (T) (1.0 / 2.0) * (T) M_SQRT2 * sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    V k = vec_set<V> ((T) (1.0 / 2.0) * (T) M_SQRT2);
    return k * vec_sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return (1.0 / 2.0) * (T) M_SQRT2
      * ((T) (1.0 / 2.0) * x * sqrt ((T) 2. * x * x + 1.)
         + (T) (1.0 / 4.0) * M_SQRT2 * asinh (M_SQRT2 * x));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    V sum1 = vec_set ((T) (1.0 / 2.0));
    sum1 *= x * vec_sqrt ((T) 2. * x * x + (T) 1.);

    V sum2 = vec_set ((T) (1.0 / 4.0) * (T) M_SQRT2);
    sum2 *= vec_asinh ((T) M_SQRT2 * x);

    return ((T) (1.0 / 2.0) * (T) M_SQRT2) * (sum1 + sum2);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt2_waveshaper_adaa = adaa::waveshaper<sqrt2_sigmoid_functions, order>;

//------------------------------------------------------------------------------
} // namespace artv
