#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
// clang-format off
/*
There are many approximations out there suitable for a waveshaper, but is there
one that is integrable analitically? Can the integral be approximated instead?

https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/
(x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x) / (28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135)

https://www.kvraudio.com/forum/viewtopic.php?f=33&t=388650&start=45

inline double vox_fasttanh2( const double x )
{
	const double ax = fabs( x );
	const double x2 = x * x;

	return( x * ( 2.45550750702956 + 2.45550750702956 * ax + ( 0.893229853513558 + 0.821226666969744 * ax ) * x2 ) / ( 2.44506634652299 + ( 2.44506634652299 + x2 ) * fabs( x + 0.814642734961073 * x * ax )));
}

Octave reminder:

x = sym('x') # symbolic
fn = x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x
fd = 28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135
f = fn / fd
ccode (horner (fn, x))
ccode (horner (fd, x))

At the end this was unused.


Tan based on exp?

tanh = exp(x) - exp(-x) / exp(x) + exp (-x)
cosh = 0.5 * (exp(x) + exp(-x))

Notice that:
exp(-x) = 1/exp(x)

So for ADAA it is needed:

-exp (x)
-exp (-x)
-exp (x1)
-exp (-x1)
-exp (x1*x*0.5).
-exp (-(x1*x*0.5)).

exp(x1) and exp(-x1) are stored from the previous round.

So 4 "exp" calls are needed, alternatively 2 exp calls and 2 divisions.
*/

// clang-format on
struct tanh_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
#if !XSIMD_BROKEN_W_FAST_MATH
    return xsimd::tanh (x);
#else
    simd_batch<T, N> r;
    for (uint i = 0; i < simd_batch<T, N>::size; ++i) {
      r[i] = tanh (x[i]);
    }
#endif
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return log (cosh (x));
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
#if !XSIMD_BROKEN_W_FAST_MATH
    return xsimd::log (xsimd::cosh (x));
#else
    simd_batch<T, N> r;
    for (uint i = 0; i < simd_batch<T, N>::size; ++i) {
      r[i] = log (cosh (x[i]));
    }
#endif
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    // not solveable analitically in a practical/easy way.
    static_assert (!std::is_same_v<T, T>, "TBI?");
    return (T) 0.;
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    // not solveable analitically in a practical/easy way.
    static_assert (!std::is_same_v<T, T>, "TBI?");
    return {(T) 0.};
  }
  //----------------------------------------------------------------------------
};
#if 1
//------------------------------------------------------------------------------
// Try to use the exp() based hyperbolic/trigonometric functions
//
// cosh = exp(x) + exp(-x) / 2;
// tanh = (exp(2x) + 1) / (exp(2x) - 1);
//
// Regular version is rewritten as (notice that two "-log(2)" cancel):
//
// (log (exp(x) + 1/exp(x)) - log (exp(x1) + 1/exp(x1))) / (x - x1)
//
// Fallback (notice that the 2x on exp cancels by the "1/2" on  (x + x1) / 2):
//
// (exp(x)*exp(x1) + 1) / (exp(x)*exp(x1) - 1);
//------------------------------------------------------------------------------
class tanh_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_exp, x1_int, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void init_states (crange<T> st)
  {
    static_assert (std::is_floating_point_v<T>, "");

    assert (st.size() >= n_states);

    st[x1]     = (T) 0.;
    st[x1_exp] = (T) 1.; // exp (0)
    st[x1_int] = (T) 0.6931471805599453; // log (2)
  }
  //----------------------------------------------------------------------------
  template <size_t simd_bytes, class T>
  static void init_states_multi_aligned (crange<T> st)
  {
    static_assert (std::is_floating_point_v<T>, "");
    using batch                      = simd_reg<T, simd_bytes>;
    static constexpr auto n_builtins = batch::size;

    assert (st.size() >= n_builtins * n_states);

    T* x1v_ptr     = &st[x1 * n_builtins];
    T* x1_expv_ptr = &st[x1_exp * n_builtins];
    T* x1_intv_ptr = &st[x1_int * n_builtins];

    auto x1v     = batch {(T) 0.};
    auto x1_expv = batch {(T) 1.}; // exp (0)
    auto x1_intv = batch {(T) 0.6931471805599453}; // log (2)

    x1v.store_aligned (x1v_ptr);
    x1_expv.store_aligned (x1_expv_ptr);
    x1_intv.store_aligned (x1_intv_ptr);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    static_assert (std::is_floating_point_v<T>, "");

    assert (st.size() >= n_states);

    T x1v     = st[x1];
    T x1_expv = st[x1_exp];
    T x1_intv = st[x1_int];

    // avoid exp going out of range by clipping the input. clipping at a point
    // that allows "x_exp * x1_expv" to not result on infinity.
    if constexpr (std::is_same_v<T, double>) {
      x = std::clamp (x, -300., 300.);
    }
    else {
      x = std::clamp (x, -150., 150.); // untested...
    }

    T x_exp = exp (x);
    T x_int = log (x_exp + ((T) 1. / x_exp));

    st[x1]     = x;
    st[x1_exp] = x_exp;
    st[x1_int] = x_int;

    T diff = x - x1v;
    T ret;
    if (abs (diff) >= adaa::epsilon (T {})) {
      ret = (x_int - x1_intv) / diff;
    }
    else {
      ret = (x_exp * x1_expv - (T) 1.) / (x_exp * x1_expv + (T) 1);
    }
    return ret;
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

    if constexpr (std::is_same_v<T, double>) {
      x = xsimd::clip (x, batch {(T) -300.}, batch {(T) 300.});
    }
    else {
      x = xsimd::clip (x, batch {(T) -150.}, batch {(T) 150.}); // untested
    }

    T* x1v_ptr     = &st[x1 * n_builtins];
    T* x1_expv_ptr = &st[x1_exp * n_builtins];
    T* x1_intv_ptr = &st[x1_int * n_builtins];

    batch x1v {x1v_ptr, xsimd::aligned_mode {}};
    batch x1_expv {x1_expv_ptr, xsimd::aligned_mode {}};
    batch x1_intv {x1_intv_ptr, xsimd::aligned_mode {}};

    static_assert (!std::is_same_v<T, T>, "XSIMD broken on ffast-math. TODO");

#if !XSIMD_BROKEN_W_FAST_MATH
    batch x_exp = xsimd::exp (x);
    batch x_int = xsimd::log (x_exp + ((T) 1. / x_exp));
#else
    batch x_exp, x_int;
    for (uint i = 0; i < n_builtins; ++i) {
      x_exp[i] = exp (x);
      x_int[i] = log (x_exp + ((T) 1. / x_exp));
    }
#endif

    x.store_aligned (x1v_ptr);
    x_exp.store_aligned (x1_expv_ptr);
    x_int.store_aligned (x1_intv_ptr);

    batch num      = x_int - x1_intv;
    batch den      = x - x1v;
    batch fallback = x_exp * x1_expv - (T) 1.;

    auto no_fallback = xsimd::abs (den) > batch {adaa::epsilon (T {})};

    num = xsimd::select (no_fallback, num, fallback);
    den = xsimd::select (no_fallback, den, fallback + (T) 2.);
    return num / den;
  }
};
//------------------------------------------------------------------------------
template <uint order>
using tanh_adaa = std::conditional_t<
  order == 1,
  tanh_adaa_1,
  adaa::waveshaper<tanh_functions, order>>;
//------------------------------------------------------------------------------
#else
//------------------------------------------------------------------------------
template <uint order>
using tanh_adaa = adaa::waveshaper<tanh_functions, order>;
//------------------------------------------------------------------------------
#endif

//------------------------------------------------------------------------------
} // namespace artv
