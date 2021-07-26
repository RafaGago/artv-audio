#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct pow2_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    return xsimd::abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return xsimd::abs (x) * x * x * (T) (1. / 3.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return xsimd::abs (x) * x * x * (T) (1. / 3.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return x * x * x * x * (T) (1. / 12.);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    return x * x * x * x * (T) (1. / 12.);
  }
  //----------------------------------------------------------------------------
};
#if 1
// TODO: is the sign preservation breaking this?
//------------------------------------------------------------------------------
// This one has a trivial simplification, it doesn't require branches.
//
// The formula is:
// ((1/6) * x**3) - ((1/6) * x1**3) / x - x1
//
// It can easily be factored to:
//
// (1/6) * (x**2 + x*x1 + x1**2)
//------------------------------------------------------------------------------
class pow2_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_pow2, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    T x1v      = st[x1];
    T x1_pow2v = st[x1_pow2];

    T xpow2 = x * x;

    st[x1]      = x;
    st[x1_pow2] = xpow2;

    return abs (xpow2 + (x1v * x) + x1_pow2v)
      * sgn_no_zero (x, (T) -1. / 6., (T) 1. / 6.);
  }
  //----------------------------------------------------------------------------
  template <size_t simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>               st,
    simd_reg<T, simd_bytes> x)
  {
    static_assert (std::is_floating_point_v<T>, "");
    using batch                      = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = batch::size;

    assert (st.size() >= n_builtins * n_states);

    T* x1v_ptr     = &st[x1 * n_builtins];
    T* x1_pow2_ptr = &st[x1_pow2 * n_builtins];

    batch x1v {x1v_ptr, xsimd::aligned_mode {}};
    batch x1_pow2v {x1_pow2_ptr, xsimd::aligned_mode {}};

    batch xpow2 = x * x;

    x.store_aligned (x1v_ptr);
    xpow2.store_aligned (x1_pow2_ptr);

    return xsimd::abs (xpow2 + (x1v * x) + x1_pow2v)
      * sgn_no_zero (x, batch {(T) -1. / 6.}, batch {(T) 1. / 6.});
  }
};
//------------------------------------------------------------------------------
template <uint order>
using pow2_adaa = std::conditional_t<
  order == 1,
  pow2_adaa_1,
  adaa::waveshaper<pow2_functions, order>>;
//------------------------------------------------------------------------------
#else
template <uint order>
using pow2_adaa = adaa::waveshaper<pow2_functions, order>;
#endif
} // namespace artv
