#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct sqrt_sin_waveshaper_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    // return tanh (x);
    return (x / sqrt (x * x + T (1.))) + (sin (16 * x)) * (1. / 16.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return (x / xsimd::sqrt (x * x + batch {(T) 1.}))
      + xsimd::sin (x * batch {(T) 16.}) * batch {(T) 1. / 16.};
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return sqrt (x * x + T (1.)) - (cos (16 * x)) * (1. / 256.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return (xsimd::sqrt (x * x + batch {(T) 1.}))
      - xsimd::cos (x * batch {(T) 16.}) * batch {(T) 1. / 256.};
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return ((x * sqrt (x * x + T (1.))) + asinh (x)) * T (0.5)
      - (sin (16 * x)) * (1. / 4096.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return ((x * xsimd::sqrt (x * x + batch {(T) 1.}))
            - xsimd::asinh (x) * xsimd::batch {(T) 0.5})
      - xsimd::sin (x * batch {(T) 16.}) * batch {(T) 1. / 4096.};
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_sin_waveshaper_adaa
  = adaa::waveshaper<sqrt_sin_waveshaper_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
