#pragma once

#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_bits.hpp"
#include "artv-common/misc/vec_util.hpp"
#include <limits>

namespace artv {
// approximations (might belong sowhere else)
//------------------------------------------------------------------------------
// https://www.kvraudio.com/forum/viewtopic.php?t=530167
// For use with single-precision floating point vectors mostly.
//
// The input is taken as if pi where multiplied to it, so:
// vec_tan_pix_prewarp(x) ~= tan(x * M_PI)
// range [1e-35, 0.4999999]
// Error is around 3.6 ULPs, 2.3e-07 relative error.
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline V vec_tan_pix_prewarp_2dat (V pi_x)
{
  // R(x^2)*x is relative error minmax rational fit to tan(x*pi):
  // R(x) =
  // (-3.13832354545593x+3.31913113594055)/(x^2+-4.47475242614746x+1.05651223659515)

  // We use this useful identity to compute on finite range of values:
  // tan(x - pi/2)  == 1 / -tan(x)
  // Since we compute tan(x*pi), identity is tanpi(x-0.5) == 1 / -tanpi(x)
  using T                = vec_value_type_t<V>;
  constexpr bool k_round = sizeof (T) == 4;
  constexpr T    k1      = k_round ? 3.319131136 : 3.31913113594055;
  constexpr T    k2      = k_round ? -3.138323545 : -3.13832354545593;
  constexpr T    k3      = k_round ? -4.474752426 : -4.47475242614746;
  constexpr T    k4      = k_round ? 1.056512237 : 1.05651223659515;

  for (uint i = 0; i < vec_traits_t<V> {}.size; ++i) {
    assert (pi_x[i] >= std::numeric_limits<T>::epsilon());
    assert (pi_x[i] <= 0.4999999);
  }

  // Compute the mask, whether we compute using the identity (x>0.25) or
  // directly (x<0.25)
  // Force the value into the finite range (-0.25,0.25)
  // Numerator and Denominator use Horner's scheme.
  auto mask = pi_x > T {0.25};
  pi_x -= mask ? vec_set<V> (T {0.5}) : vec_set<V> (T {0.});

  V f2  = pi_x * pi_x;
  V num = vec_abs (pi_x) * (k1 + f2 * k2);
  V den = f2 * (f2 + k3) + k4;
  vec_xorswap (mask, num, den);
  return num / den;
}
//------------------------------------------------------------------------------
// https://www.kvraudio.com/forum/viewtopic.php?p=7524831#p7524831
// Error is around 15-17 ULPs for these functions on
// (std::numeric_limits<float>::min(), std::numeric_limits<float>::max()).
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_log2_2dat (V x)
{

  using T  = vec_value_type_t<V>;
  using U  = typename vec_traits_t<V>::same_size_uint_type;
  using S  = typename vec_traits_t<V>::same_size_int_type;
  using TU = vec_value_type_t<U>;

  static_assert (std::is_same_v<T, float>, "Only implemented for float");
  constexpr TU mantissa_mask = (TU) ~((1u << 23u) - 1u);

  for (uint i = 0; i < vec_traits_t<V> {}.size; ++i) {
    // No negative and NaN checks! No zero.
    assert (x[i] >= std::numeric_limits<T>::epsilon());
  }

  S spl_exp = vec_bit_cast<S> (x * T {1.41421356237});
  spl_exp -= vec_bit_cast<S> (T {1.0});
  spl_exp &= mantissa_mask;
  V mantissa      = vec_bit_cast<V> (vec_bit_cast<S> (x) - spl_exp);
  V log2_exponent = vec_cast<V> (ashr (spl_exp, 23));

  V num = ((mantissa + T {1.929443550e+01}) * mantissa) + T {1.011593342e+01};
  num *= mantissa - T {1.};

  V den = ((T {6.316540241e+00} * mantissa) + T {1.266638851e+01}) * mantissa;
  den += T {2.095932245e+00};

  return (num / den) + log2_exponent;
}

template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_log_2dat (V x)
{
  return vec_log2_2dat (x) * vec_set<V> (0.69314718056); // round(1 / log2(e))
}

template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_log10_2dat (V x)
{
  return vec_log2_2dat (x) * vec_set<V> (0.301029995664); // round(1 / log2(10))
}
//------------------------------------------------------------------------------
// exp_mp https://www.kvraudio.com/forum/viewtopic.php?t=510794&hilit=exp
// exp_mp: 3.27 ulp, simplified range reduction.
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_exp_2dat (V x)
{
  using T  = vec_value_type_t<V>;
  using U  = typename vec_traits_t<V>::same_size_uint_type;
  using S  = typename vec_traits_t<V>::same_size_int_type;
  using TU = vec_value_type_t<U>;

  static_assert (std::is_same_v<T, float>, "Only implemented for float");

  //[-0.5,0.5] 2^x approx polynomial ~ 2.4 ulp
  V f = x * T {1.442695041}; // f = log2(e) *x
  S i = vec_sint_round (f);
  f -= vec_cast<V> (i);

  // estrin-horner evaluation scheme
  V f2 = f * f;
  V p  = f2 * (T {9.675540961e-03} + f * T {1.327647245e-03});
  p += T {2.402212024e-01} + f * T {5.550713092e-02};
  p *= f2;
  p += T {1.000000119e+00} + f * T {6.931469440e-01};

  i = vec_bit_cast<S> (vec_bit_cast<U> (i) << 23); // logical shift
  return vec_bit_cast<V> (i + vec_bit_cast<S> (p));
}

template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_exp2_2dat (V x)
{
  return vec_exp_2dat (x * vec_set<V> (M_LN2));
}

template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_exp10_2dat (V x)
{
  return vec_exp_2dat (x * vec_set<V> (M_LN10));
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_pow_2dat (V x, V y)
{
  return vec_exp2_2dat (vec_log2_2dat (x) * y);
}
//------------------------------------------------------------------------------
// https://www.kvraudio.com/forum/viewtopic.php?t=519728&hilit=exp
// 3.0507440678775311 ulps on [3 * std::numeric_limits<float>::min(), 300]
// ~3e-7 relative
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_tanh_ln2_div2_2dat (V x)
{
  using T  = vec_value_type_t<V>;
  using U  = typename vec_traits_t<V>::same_size_uint_type;
  using S  = typename vec_traits_t<V>::same_size_int_type;
  using TU = vec_value_type_t<U>;

  static_assert (std::is_same_v<T, float>, "Only implemented for float");

  V f = vec_clamp (x, T {30.});
  S i = vec_sint_round (f);
  f -= vec_cast<V> (i);
  V f2 = f * f;
  V p  = f2 * (T {1.339077600e-03} + f * T {1.540359954e-04});
  p += T {5.550357327e-02} + f * T {9.618237615e-03};
  p *= f2;
  p += T {6.931471825e-01} + f * T {2.402264923e-01};
  p *= f;

  i <<= 23;
  V biased_expm = vec_bit_cast<V> (vec_bit_cast<S> (T {1.}) - i);
  V exp_cor     = T {1.} - biased_expm;
  V exp_cor_p   = T {1.} + biased_expm;
  V exp2xm1     = p + exp_cor;
  V exp2xp1     = p + exp_cor_p;
  return exp2xm1 / exp2xp1;
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
V vec_tanh_2dat (V x)
{
  return vec_tanh_ln2_div2_2dat (vec_value_type_t<V> {M_LN2 * 0.5});
}
//------------------------------------------------------------------------------

} // namespace artv
