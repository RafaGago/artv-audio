#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
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
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    batch num   = batch {(T) M_SQRT2} * x;
    batch den   = xsimd::sqrt (batch {(T) 2.} * x * x + batch {(T) 1.});
    return num / den;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return (T) (1.0 / 2.0) * (T) M_SQRT2 * sqrt ((T) 2. * x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    batch k {(T) (1.0 / 2.0) * (T) M_SQRT2};
    return k * xsimd::sqrt (batch {(T) 2.} * x * x + batch {(T) 1.});
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
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;

    batch sum1 {(T) (1.0 / 2.0)};
    sum1 *= x * xsimd::sqrt (batch {(T) 2.} * x * x + batch {(T) 1.});

    batch sum2 = {(T) (1.0 / 4.0) * (T) M_SQRT2};
    sum2 *= xsimd::asinh (batch {(T) M_SQRT2} * x);

    return batch {(T) (1.0 / 2.0) * (T) M_SQRT2} * (sum1 + sum2);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt2_waveshaper_adaa = adaa::waveshaper<sqrt2_sigmoid_functions, order>;

//------------------------------------------------------------------------------
} // namespace artv
