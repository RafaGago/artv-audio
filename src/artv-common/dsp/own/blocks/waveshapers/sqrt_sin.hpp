#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct sqrt_sin_sigmoid_functions {
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
    batch sq, sn;
    sq = x / xsimd::sqrt (x * x + batch {(T) 1.});
#if !XSIMD_BROKEN_W_FAST_MATH
    sn = xsimd::sin (x * (T) 16.) * (T) (1. / 16.);
#else
    for (uint i = 0; i < batch::size; ++i) {
      sn[i] = sin (x[i] * (T) 16.) * (T) (1. / 16.);
    }
#endif
    return sq + sn;
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
    using batch = simd_batch<T, N>;
    batch sq, cs;
    sq = xsimd::sqrt (x * x + (T) 1.);
#if !XSIMD_BROKEN_W_FAST_MATH
    cs = xsimd::cos (x * (T) 16.) * ((T) 1. / 256.);
#else
    for (uint i = 0; i < batch::size; ++i) {
      cs[i] = cos (x[i] * (T) 16.) * ((T) 1. / 256.);
    }
#endif
    return sq + cs;
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
    using batch = simd_batch<T, N>;
    batch sq, as, sn;
    sq = xsimd::sqrt (x * x + (T) 1.) * x;
#if !XSIMD_BROKEN_W_FAST_MATH
    as = xsimd::asinh (x) * (T) 0.5;
    sn = xsimd::sin (x * (T) 16.) * (T) (1. / 4096.);
#else
    for (uint i = 0; i < batch::size; ++i) {
      as = asinh (x[i]) * (T) 0.5;
      sn = sin (x[i] * (T) 16.) * (T) (1. / 4096.);
    }
#endif
    return sq + as - sn;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_sin_sigmoid_adaa
  = adaa::waveshaper<sqrt_sin_sigmoid_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
