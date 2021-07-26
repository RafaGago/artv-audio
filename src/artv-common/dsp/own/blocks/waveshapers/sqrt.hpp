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
struct sqrt_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return sqrt (abs (x)) * sgn_no_zero (x);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    return xsimd::sqrt (xsimd::abs (x)) * sgn_no_zero (x);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return sqrt (abs (x)) * x * sgn_no_zero (x, (T) (-2. / 3.), (T) (2. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    return xsimd::sqrt (x) * x * sgn_no_zero (x, (T) (-2. / 3.), (T) (2. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return sqrt (x) * x * x * sgn_no_zero (x, (T) (-4. / 15.), (T) (4. / 15.));
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    return xsimd::sqrt (x) * x * x
      * sgn_no_zero (x, (T) (-4. / 15.), (T) (4. / 15.));
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// No unstable division for ADAA on sqrt(x). Simple math, but I prefer to write
// a reminder
//
// 1st order ADDA formula:
//
//  ((3/2)(x**(3/2) - x1**(3/2))) /  (x - x1)
//
// Multiplying numerator and denominator by: "((3/2)(x**(3/2) - x1**(3/2)))"
//
// On Octave, with the symbolic package loaded and the x and x1 vars created we:
// attack the numerator:
//
// expand (((3/2)*(x**(3/2)) - (3/2)*(x1**(3/2))) * ((3/2)*(x**(3/2)) +
// (3/2)*(x1**(3/2))))
//
// Which returns:
//
// (9/4)(x**3 - x1**3)
//
// Expanding it:
//
// factor ((9/4)*(x**3 - x1**3))
//
// Which returns:
//
// (9/4) * (x - x1) * (x**2 + x*x1 + x1**2)
//
// So (x - x1) goes away from the denominator.
//
// simplify (9*(x**2 + x*x1 + x1**2)) / (4 * ((3/2)*(x**(3/2) + x1**(3/2))))
//
// (9 * (x**2 + x*x1 + x1**2)) / (6 * (x**3/2 + x1**3/2))
//
// x**3/2 is "x(2/2)*x(1/2)" so it becomes x * sqrt (x)
//
// Final. This is heavier but it doesn't have branching to care about;
// 1 sqrt + 1 div (as the values for the prev sample are stored). It just adds
// low costs operations: sums and multiplies. The division and sqrt was going to
// be there anyways.
//
// (9 * (x**2 + x*x1 + x1**2)) / (6 * (x * sqrt (x) + x1 * sqrt (x1))
//------------------------------------------------------------------------------
// TODO: is the sign preservation breaking this?
class sqrt_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_sqrt, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    T x1v      = st[x1];
    T x1_sqrtv = st[x1_sqrt];
    T x_sqrtv  = sqrt (x) * x * sgn_no_zero (x);

    st[x1]      = x;
    st[x1_sqrt] = x_sqrtv;

    T num = (x * x + x * x1v + x1v * x1v) * ((T) 9.);
    T den = (x_sqrtv + x1_sqrtv) * ((T) 6.);

    return num / den;
  }
  //----------------------------------------------------------------------------
  template <size_t sse_bytes, class T>
  static simd_reg<T, sse_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>              st,
    simd_reg<T, sse_bytes> x)
  {
    static_assert (std::is_floating_point_v<T>, "");
    using batch                      = simd_reg<T, sse_bytes>;
    static constexpr auto n_builtins = batch::size;

    assert (st.size() >= n_builtins * n_states);

    T* x1v_ptr      = &st[x1 * n_builtins];
    T* x1v_sqrt_ptr = &st[x1_sqrt * n_builtins];

    batch x1v {x1v_ptr, xsimd::aligned_mode {}};
    batch x1vsqrt {x1v_sqrt_ptr, xsimd::aligned_mode {}};

    batch xsqrt = xsimd::sqrt (x) * x * sgn_no_zero (x);

    x.store_aligned (x1v_ptr);
    xsqrt.store_aligned (x1v_sqrt_ptr);

    batch num = (x * x + x * x1v + x1v * x1v) * ((T) 9.);
    batch den = (xsqrt + x1vsqrt) * ((T) 6.);

    return num / den;
  }
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_adaa = std::conditional_t<
  order == 1,
  sqrt_adaa_1,
  adaa::waveshaper<sqrt_functions, order>>;
//------------------------------------------------------------------------------
} // namespace artv
