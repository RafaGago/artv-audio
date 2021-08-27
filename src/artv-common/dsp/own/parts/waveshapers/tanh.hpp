#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    return vec_tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    return vec_log (vec_cosh (x));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int2_fn (V x)
  {
    // not solveable analitically in a practical/easy way.
    static_assert (!std::is_same_v<V, V>, "TBI?");
    return x;
  }
  //----------------------------------------------------------------------------
};
#if 0
// This seems right analitically, but it seems to have either a bug or precision
// issues
//------------------------------------------------------------------------------
// This tries to use the exp() based hyperbolic/trigonometric functions and
// buffering results.
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
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>>)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (traits.size * n_states));

    T* x1v_ptr     = &st[x1 * traits.size];
    T* x1_expv_ptr = &st[x1_exp * traits.size];
    T* x1_intv_ptr = &st[x1_int * traits.size];

    auto x1v     = vec_set ((T) 0.);
    auto x1_expv = vec_set ((T) 1.); // exp (0)
    auto x1_intv = vec_set ((T) 0.6931471805599453); // log (2)

    vec_store (x1v_ptr, x1v);
    vec_store (x1_expv_ptr, x1_expv);
    vec_store (x1_intv_ptr, x1_intv);
  }
  //----------------------------------------------------------------------------
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    using T = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= traits.size * n_states);

    // avoid exp going out of range by clipping the input. clipping at a point
    // that allows "x_exp * x1_expv" to not result on infinity.
    if constexpr (std::is_same_v<T, double>) {
      x = vec_clamp (x, (T) -300., (T) 300.);
    }
    else {
      x = vec_clamp (x, (T) -150.,
                     (T) 150.); // untested
    }

    T* x1v_ptr     = &st[x1 * traits.size];
    T* x1_expv_ptr = &st[x1_exp * traits.size];
    T* x1_intv_ptr = &st[x1_int * traits.size];

    V x1v     = vec_load<V> (x1v_ptr);
    V x1_expv = vec_load<V> (x1_expv_ptr);
    V x1_intv = vec_load<V> (x1_intv_ptr);

    V x_exp = vec_exp (x);
    V x_int = vec_log (x_exp + ((T) 1. / x_exp));

    vec_store (x1v_ptr, x);
    vec_store (x1_expv_ptr, x_exp);
    vec_store (x1_intv_ptr, x_int);

    V num      = x_int - x1_intv;
    V den      = x - x1v;
    V fallback = x_exp * x1_expv - (T) 1.;

    auto no_fallback = vec_abs (den) > adaa::epsilon (T {});

    num = no_fallback ? num : fallback;
    den = no_fallback ? den : fallback + (T) 2.;
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
