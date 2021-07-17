#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// sqrt waveshaper with first derivative AA. Based on this thread:
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=521377&sid=18fa45082ea11269f4c29c3ccf36becd
struct sqrt_sigmoid_functions {
  template <class T>
  static T fn (T x)
  {
    return x / sqrt (x * x + T (1.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int_fn (T x)
  {
    return sqrt (x * x + T (1.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int2_fn (T x)
  {
    return ((x * sqrt (x * x + T (1.))) + asinh (x)) * T (0.5);
  }
};
//------------------------------------------------------------------------------
// This is special cased, so it doesn't have to branch on the problematic
// division:
// https://www.kvraudio.com/forum/viewtopic.php?t=521377&start=30
// message by "martinvicanek"
// (sqrt(1 + x^2) - sqrt(1+ x1^2))/(x - x1) =
//    (x + x1)/(sqrt(1 + x^2) + sqrt(1 + x1^2))

class sqrt_waveshaper_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_sqrt, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    T x1v     = st[x1];
    T x1vsqrt = st[x1_sqrt];
    T xsqrt   = sqrt ((T) 1. + x * x);

    st[x1]      = x;
    st[x1_sqrt] = xsqrt;

    return (x + x1v) / (xsqrt + x1vsqrt);
  }
#if 0
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_aligned (
    crange<const T>,
    crange<T>       z,
    crange<const T> x)
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg             = simd_reg<T, simd_bytes>;
    constexpr auto n_builtins = simdreg::size;

    assert (z.size() >= n_builtins * n_states);
    assert (x.size() >= n_builtins);

    simdreg x_v {x.data(), xsimd::aligned_mode {}};
    simdreg x1_v {z.data(), xsimd::aligned_mode {}};
    // this SIMD version may not translate to faster code depending on how the
    // simple version is translated (if the machine has harware inverse
    // reciprocal).
    auto ret = (x_v + x1_v)
      / (xsimd::sqrt (x_v * x_v + (T) 1.) + xsimd::sqrt (x1_v * x1_v + (T) 1.));

    x_v.store_aligned (z.data());
    return ret;
  }
  //----------------------------------------------------------------------------
#endif
};

template <uint order>
using sqrt_waveshaper_adaa = std::conditional_t<
  order == 1,
  sqrt_waveshaper_adaa_1,
  adaa::waveshaper<sqrt_sigmoid_functions, order>>;
//------------------------------------------------------------------------------
} // namespace artv
