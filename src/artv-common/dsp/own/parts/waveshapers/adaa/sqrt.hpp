#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct sqrt_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    auto ax = vec_abs (x);
    return vec_sqrt (ax) * vec_sgn_no_zero (x);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    auto ax = vec_abs (x);
    return vec_sqrt (ax) * ax * (T) (2. / 3.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    auto ax = vec_abs (x);
    return vec_sqrt (ax) * ax * x * (T) (4. / 15.);
  }
  //----------------------------------------------------------------------------
};
#if 0
// This simplification loses the sign info. Can be fixed? The integral is always
// positive.
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
// So (x - x1) goes away from the numerator and denominator.
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
class sqrt_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_sqrt, n_states };
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>>)
  {
  }
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V x)
  {
    static_assert (std::is_floating_point_v<T>, "");
    using batch                      = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = batch::size;

    assert (st.size() >= n_builtins * n_states);

    T* x1v_ptr      = &st[x1 * n_builtins];
    T* x1v_sqrt_ptr = &st[x1_sqrt * n_builtins];

    batch x1v {x1v_ptr, xsimd::aligned_mode {}};
    batch x1vsqrt {x1v_sqrt_ptr, xsimd::aligned_mode {}};

    auto ax = xsimd::abs (x);
    // likely avoid NaN, this value is almost out of resolution.
    ax += batch {std::numeric_limits<T>::min()};
    batch xsqrt = xsimd::sqrt (ax) * x;
    // likely avoid division by zero, this value will get out of resolution.
    xsqrt += batch {std::numeric_limits<T>::min()};

    x.store_aligned (x1v_ptr);
    xsqrt.store_aligned (x1v_sqrt_ptr);

    batch num = (x * x + x * x1v + x1v * x1v);
    num *= sgn_no_zero (x, batch {(T) -9.}, batch {(T) 9.});
    batch den = (xsqrt + x1vsqrt) * ((T) 6.);
    return (num / den);
  }
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_adaa = std::conditional_t<
  order == 1,
  sqrt_adaa_1,
  adaa::waveshaper<sqrt_functions, order>>;
//------------------------------------------------------------------------------
#else
//------------------------------------------------------------------------------
template <uint order>
using sqrt_adaa = adaa::waveshaper<sqrt_functions, order>;
//------------------------------------------------------------------------------
#endif

} // namespace artv
