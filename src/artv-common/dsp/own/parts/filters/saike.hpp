#pragma once

// Saike's excellent filters.

#include <cmath>

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"
//------------------------------------------------------------------------------
namespace artv { namespace saike {

namespace detail {

template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
inline V f_g (V v)
{
  using T = vec_value_type_t<V>;
  return vec_max ((T) -1., vec_min ((T) 1., v));
}
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
inline V f_dg (V v)
{
  // return 1. - 1. * (abs (v) > 1.);
  using T = vec_value_type_t<V>;
  return (vec_abs (v) > (T) 1.) ? vec_set<V> (0.) : vec_set<V> (1.);
}

template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
inline V f_g_asym (V v)
{
  //  v > 0 ? std::min (1., v) : std::max (-1., v * .25);
  using T = vec_value_type_t<V>;
  return (v > (T) 0) ? vec_min ((T) 1., v) : vec_max ((T) -1., v * (T) .25);
}

template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
inline V f_dg_asym (V v)
{
  // v > 0 ? 1. - 1. * (std::abs (v) > 1) : .25 - .25 * (std::abs (v) > 4.);
  using T = vec_value_type_t<V>;
  return v > (T) 0 ? ((v > (T) 1.) ? vec_set<V> (0.) : vec_set<V> (1.))
                   : ((v < (T) -4.) ? vec_set<V> (0.) : vec_set<V> (.25));
}
//##############################################################################
// MS20 /K35: base
//##############################################################################
struct ms20_base {
  // From Saike's FM Filter2 (Yutani).
  // SHA b1c1952c5f798d968ff15714bca6a5c694e06d82
  // https://github.com/JoepVanlier/JSFX
  //----------------------------------------------------------------------------
  enum coeffs { hh, k, n_coeffs };
  enum state { y1, y2, d1, d2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    freq           = vec_min (sr < (T) 80000. ? (T) 14000. : (T) 20000., freq);
    T oversampling = (sr * (T) (1. / 44100.));
    V norm_freq    = (freq * (T) (1. / 20000.));

    V f0 = norm_freq * M_PI / oversampling;
    V hv = vec_tan (f0 / ((T) 2.1 * oversampling)) * (T) 2.1 * oversampling;
    vec_store (&c[hh * traits.size], hv * (T) 0.5);
    vec_store (&c[k * traits.size], reso * (T) 2.);
  }

  static constexpr double epsilon  = 0.00000001;
  static constexpr uint   max_iter = 6;
};
} // namespace detail
//##############################################################################
// MS20 /K35: "MS20_nonlin_tanh" in the original sources.
//##############################################################################
struct ms20_lowpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V gd2k      = detail::f_g (d2_v * k_v);
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v + in - gd2k);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2          = k_v * y2_v;
      V gky2         = detail::f_g (ky2);
      V dgky2        = detail::f_dg (ky2);
      V sig1         = in - y1_v - gky2;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = -hhthsig1sqm1 + (T) 1.;
      V b            = -k_v * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (k_v * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return d2_v;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_highpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V kc        = .85 * k_v;
    V gkd2px    = detail::f_g (kc * (d2_v + in));
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v - gkd2px);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v - in + gkd2px);

    for (uint i = 0; i < max_iter; ++i) {
      V kxpy2        = kc * (in + y2_v);
      V gkxpy2       = detail::f_g (kxpy2);
      V dgky2px      = detail::f_dg (kxpy2);
      V sig1         = -y1_v - gkxpy2;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = -in + y1_v - y2_v + gkxpy2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = -hhthsig1sqm1 + (T) 1.;
      V b            = -kc * hhthsig1sqm1 * dgky2px;
      V c            = hhthsig2sqm1;
      V d            = (kc * dgky2px - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return y2_v + in;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_bandpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V kc        = (T) .95 * k_v;
    V gd2k      = detail::f_g (d2_v * kc);
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v - in - gd2k);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + in + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2          = kc * y2_v;
      V gky2         = detail::f_g (ky2);
      V dgky2        = detail::f_dg (ky2);
      V sig1         = -in - y1_v - gky2;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = in + y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = (T) 1. - hhthsig1sqm1;
      V b            = -kc * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (kc * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return d2_v;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_notch : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V gd2k      = detail::f_g (d2_v * k_v);
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v - in - gd2k);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + in + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2          = k_v * y2_v;
      V gky2         = detail::f_g (ky2);
      V dgky2        = detail::f_dg (ky2);
      V sig1         = -in - y1_v - gky2;
      V thsig1       = vec_tanh (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = in + y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = (T) 1. - hhthsig1sqm1;
      V b            = -k_v * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (k_v * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return in - y2_v;
  }
  //----------------------------------------------------------------------------
};
//##############################################################################
// MS20 /K35: "MS20_nonlin_asym_tanh" in the original sources.
//##############################################################################
//----------------------------------------------------------------------------
struct ms20_asym_lowpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V gd2k = detail::f_g_asym (d2_v * k_v);
    V sig1 = -d1_v + in - gd2k;
    V qq   = (sig1 < 0.) ? (T) .6 : (T) 1.0;

    sig1 *= qq;
    V tanhterm1 = vec_tanh_approx_vaneev (sig1);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2   = k_v * y2_v;
      V gky2  = detail::f_g_asym (ky2);
      V dgky2 = detail::f_dg_asym (ky2);
      sig1    = in - y1_v - gky2;
      qq      = (sig1 < 0.) ? (T) .6 : (T) 1.0;
      sig1 *= qq;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = -qq * hhthsig1sqm1 + (T) 1.;
      V b            = -qq * k_v * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (k_v * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return d2_v;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_highpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V kc        = k_v;
    V gkd2px    = detail::f_g_asym (kc * (d2_v + in));
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v - gkd2px);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v - in + gkd2px);

    for (uint i = 0; i < max_iter; ++i) {
      V kxpy2        = kc * (in + y2_v);
      V gkxpy2       = detail::f_g_asym (kxpy2);
      V dgky2px      = detail::f_dg_asym (kxpy2);
      V sig1         = -y1_v - gkxpy2;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = -in + y1_v - y2_v + gkxpy2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = -hhthsig1sqm1 + (T) 1.;
      V b            = -kc * hhthsig1sqm1 * dgky2px;
      V c            = hhthsig2sqm1;
      V d            = (kc * dgky2px - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1.0 / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return y2_v + in;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_bandpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V kc   = k_v;
    V gd2k = detail::f_g_asym (d2_v * kc);
    V sig1 = -d1_v - in - gd2k;
    V qq   = (sig1 < 0.) ? (T) .6 : (T) 1.0;
    sig1 *= qq;
    V tanhterm1 = vec_tanh_approx_vaneev (sig1);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + in + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2   = kc * y2_v;
      V gky2  = detail::f_g_asym (ky2);
      V dgky2 = detail::f_dg_asym (ky2);
      sig1    = -in - y1_v - gky2;
      qq      = (sig1 < 0.) ? (T) .6 : (T) 1.0;
      sig1 *= qq;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = in + y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = (T) 1. - qq * hhthsig1sqm1;
      V b            = -qq * kc * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (kc * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return d2_v;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_notch : public detail::ms20_base {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr)
  {
    detail::ms20_base::reset_coeffs (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V hh_v = vec_load<V> (&co[hh * traits.size]);
    V k_v  = vec_load<V> (&co[k * traits.size]);
    V y1_v = vec_load<V> (&st[y1 * traits.size]);
    V y2_v = vec_load<V> (&st[y2 * traits.size]);
    V d1_v = vec_load<V> (&st[d1 * traits.size]);
    V d2_v = vec_load<V> (&st[d2 * traits.size]);

    V gd2k      = detail::f_g_asym (d2_v * k_v);
    V tanhterm1 = vec_tanh_approx_vaneev (-d1_v - in - gd2k);
    V tanhterm2 = vec_tanh_approx_vaneev (d1_v - d2_v + in + gd2k);

    for (uint i = 0; i < max_iter; ++i) {
      V ky2          = k_v * y2_v;
      V gky2         = detail::f_g_asym (ky2);
      V dgky2        = detail::f_dg_asym (ky2);
      V sig1         = -in - y1_v - gky2;
      V thsig1       = vec_tanh_approx_vaneev (sig1);
      V thsig1sq     = thsig1 * thsig1;
      V sig2         = in + y1_v - y2_v + gky2;
      V thsig2       = vec_tanh_approx_vaneev (sig2);
      V thsig2sq     = thsig2 * thsig2;
      V hhthsig1sqm1 = hh_v * (thsig1sq - (T) 1.);
      V hhthsig2sqm1 = hh_v * (thsig2sq - (T) 1.);
      V f1           = y1_v - d1_v - hh_v * (tanhterm1 + thsig1);
      V f2           = y2_v - d2_v - hh_v * (tanhterm2 + thsig2);
      V res          = vec_abs (f1) + vec_abs (f2);
      V a            = 1. - hhthsig1sqm1;
      V b            = -k_v * hhthsig1sqm1 * dgky2;
      V c            = hhthsig2sqm1;
      V d            = (k_v * dgky2 - (T) 1.) * hhthsig2sqm1 + (T) 1.;
      V norm         = (T) 1. / (a * d - b * c);
      y1_v -= (d * f1 - b * f2) * norm;
      y2_v -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    d1_v = y1_v;
    d2_v = y2_v;

    vec_store (&st[y1 * traits.size], y1_v);
    vec_store (&st[y2 * traits.size], y2_v);
    vec_store (&st[d1 * traits.size], d1_v);
    vec_store (&st[d2 * traits.size], d2_v);
    return in - y2_v;
  }
  //----------------------------------------------------------------------------
};
//##############################################################################
// Steiner: base
//##############################################################################
namespace detail {

struct steiner_base {
  // From Saike's FM Filter2 (Yutani).
  // SHA b1c1952c5f798d968ff15714bca6a5c694e06d82
  // https://github.com/JoepVanlier/JSFX
  //----------------------------------------------------------------------------
  enum coeffs { h, k, kh, hsq, lp, bp, hp, normalizing_const, n_coeffs };
  enum state { x, v1, v2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    reset_coeffs (co, freq, reso, 0., sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    bandpass_tag)
  {
    reset_coeffs (co, freq, reso, 0.25, sr);
  }

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    reset_coeffs (co, freq, reso, 0.5, sr);
  }

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    notch_tag)
  {
    reset_coeffs (co, freq, reso, 0.75, sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         morph,
    vec_value_type_t<V>         sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    freq           = vec_min (sr < (T) 80000. ? (T) 14000. : (T) 20000., freq);
    T oversampling = (sr * (T) (1. / 44100.));
    V norm_freq    = (freq * (T) (1. / 20000.));

    V f0 = norm_freq / oversampling;

    V h_v   = vec_tan ((T) 0.5 * M_PI * f0);
    V k_v   = (T) 3.98 * reso;
    V hsq_v = h_v * h_v; // calculated on "tick "instead to save memory?
    V kh_v  = k_v * h_v; // calculated on "tick "instead to save memory?

    vec_store (&co[h * traits.size], h_v);
    vec_store (&co[k * traits.size], k_v);
    vec_store (&co[hsq * traits.size], hsq_v);
    vec_store (&co[kh * traits.size], kh_v);
    vec_store (
      &co[normalizing_const * traits.size],
      (T) 1.0 / (-kh_v + hsq_v + (T) 2. * h_v + (T) 1.));

    /*alpha = 20.94153124476462;
    beta = 0.057872340425531923;
    vref = log(h / (alpha * beta)) / alpha;*/

    if (morph < 0.25) {
      vec_store (&co[lp * traits.size], vec_set<V> ((T) 1.0 - morph * (T) 4));
      vec_store (&co[bp * traits.size], vec_set<V> (morph * (T) 4));
      vec_store (&co[hp * traits.size], vec_set<V> ((T) 0));
    }
    else if (morph < 0.5) {
      vec_store (
        &co[hp * traits.size], vec_set<V> ((morph - (T) 0.25) * (T) 4.));
      vec_store (
        &co[bp * traits.size],
        vec_set<V> ((T) 1.0 - (morph - (T) 0.25) * (T) 4.));
      vec_store (&co[lp * traits.size], vec_set<V> ((T) 0.));
    }
    else if (morph < (T) 0.75) {
      vec_store (&co[hp * traits.size], vec_set<V> ((T) 1.0));
      vec_store (
        &co[bp * traits.size],
        vec_set<V> ((T) -2.0 * (morph - (T) 0.5) * (T) 4.));
      vec_store (
        &co[lp * traits.size], vec_set<V> ((morph - (T) 0.5) * (T) 4.));
    }
    else {
      vec_store (
        &co[hp * traits.size],
        vec_set<V> ((T) 1.0 - (morph - (T) 0.75) * (T) 4.));
      vec_store (
        &co[bp * traits.size],
        vec_set<V> ((T) -2.0 * ((T) 1 - (morph - (T) 0.75) * (T) 4.)));
      vec_store (&co[lp * traits.size], vec_set<V> ((T) 1));
    }
  }

  static constexpr double epsilon  = 0.00000001;
  static constexpr uint   max_iter = 26;
};
} // namespace detail

//##############################################################################
// Steiner1: Crazy filter!!! (eval_steiner).
//##############################################################################
struct steiner_1 : public detail::steiner_base {
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V h_v                 = vec_load<V> (&co[h * traits.size]);
    V k_v                 = vec_load<V> (&co[k * traits.size]);
    V hsq_v               = vec_load<V> (&co[hsq * traits.size]);
    V kh_v                = vec_load<V> (&co[kh * traits.size]);
    V lp_v                = vec_load<V> (&co[lp * traits.size]);
    V bp_v                = vec_load<V> (&co[bp * traits.size]);
    V hp_v                = vec_load<V> (&co[hp * traits.size]);
    V normalizing_const_v = vec_load<V> (&co[normalizing_const * traits.size]);

    V x_v  = vec_load<V> (&st[x * traits.size]);
    V v1_v = vec_load<V> (&st[v1 * traits.size]);
    V v2_v = vec_load<V> (&st[v2 * traits.size]);

    V x_xn = x_v + in;
    V v1n
      = (-kh_v * v1_v + (bp_v - hp_v) * h_v * x_xn
         + hsq_v * ((lp_v - hp_v) * x_xn - v1_v) + (T) 2. * h_v * v2_v + v1_v)
      * normalizing_const_v;
    V v2n = (kh_v * ((bp_v - hp_v) * x_xn + v2_v - (T) 2.0 * v1_v)
             + (lp_v - bp_v) * hsq_v * x_xn + (lp_v - bp_v) * h_v * x_xn
             - hsq_v * v2_v + v2_v)
      * normalizing_const_v;
    V kp1      = k_v + (T) 1.0;
    V fb       = detail::f_g ((k_v + (T) 1.) * (hp_v * x_v + v1_v));
    V s1_fixed = (T) 2.0 * bp_v * x_v - hp_v * x_v - v1_v + v2_v + fb;
    V s2_fixed = (T) -4.0 * bp_v * x_v + hp_v * x_v - (T) 2.0 * (v2_v + fb)
      + v1_v + lp_v * x_v;

    for (uint i = 0; i < max_iter; ++i) {
      V fb_s1      = kp1 * (hp_v * in + v1n);
      V fb_clipped = detail::f_g (fb_s1);
      V s1         = s1_fixed - hp_v * in - v1n + v2n + fb_clipped;
      V s2
        = s2_fixed + hp_v * in + lp_v * in + v1n - (T) 2.0 * (v2n + fb_clipped);
      V f1  = -v1n + v1_v + h_v * s1;
      V f2  = -v2n + v2_v + h_v * s2;
      V res = vec_abs (f1) + vec_abs (f2);
      V a   = h_v * (kp1 * detail::f_dg (fb_s1) - 1.0) - 1.0;
      V b   = h_v;
      V c = h_v * ((T) 1.0 - (T) 2.0 * (k_v + (T) 1.0) * detail::f_dg (fb_s1));
      V d = (T) -2.0 * h_v - (T) 1.0;
      V norm = (T) 1.0 / (a * d - b * c);
      v1n -= (d * f1 - b * f2) * norm;
      v2n -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    x_v  = in;
    v1_v = v1n;
    v2_v = v2n;

    vec_store (&st[x * traits.size], x_v);
    vec_store (&st[v1 * traits.size], v1_v);
    vec_store (&st[v2 * traits.size], v2_v);
    return (v1n + hp_v * in);
  }
  //----------------------------------------------------------------------------
};
//##############################################################################
// Steiner1: Crazy still!!! (eval_steiner_asym).
//##############################################################################
struct steiner_2 : public detail::steiner_base {
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V h_v                 = vec_load<V> (&co[h * traits.size]);
    V k_v                 = vec_load<V> (&co[k * traits.size]);
    V hsq_v               = vec_load<V> (&co[hsq * traits.size]);
    V kh_v                = vec_load<V> (&co[kh * traits.size]);
    V lp_v                = vec_load<V> (&co[lp * traits.size]);
    V bp_v                = vec_load<V> (&co[bp * traits.size]);
    V hp_v                = vec_load<V> (&co[hp * traits.size]);
    V normalizing_const_v = vec_load<V> (&co[normalizing_const * traits.size]);

    V x_v  = vec_load<V> (&st[x * traits.size]);
    V v1_v = vec_load<V> (&st[v1 * traits.size]);
    V v2_v = vec_load<V> (&st[v2 * traits.size]);

    V x_xn = x_v + in;
    V v1n
      = (-kh_v * v1_v + (bp_v - hp_v) * h_v * x_xn
         + hsq_v * ((lp_v - hp_v) * x_xn - v1_v) + (T) 2. * h_v * v2_v + v1_v)
      * normalizing_const_v;
    V v2n = (kh_v * ((bp_v - hp_v) * x_xn + v2_v - (T) 2.0 * v1_v)
             + (lp_v - bp_v) * hsq_v * x_xn + (lp_v - bp_v) * h_v * x_xn
             - hsq_v * v2_v + v2_v)
      * normalizing_const_v;
    V kp1      = k_v + (T) 1.0;
    V fb       = detail::f_g_asym ((k_v + (T) 1.) * (hp_v * x_v + v1_v));
    V s1_fixed = (T) 2.0 * bp_v * x_v - hp_v * x_v - v1_v + v2_v + fb;
    V s2_fixed = (T) -4.0 * bp_v * x_v + hp_v * x_v - (T) 2.0 * (v2_v + fb)
      + v1_v + lp_v * x_v;

    for (uint i = 0; i < max_iter; ++i) {
      V fb_s1      = kp1 * (hp_v * in + v1n);
      V fb_clipped = detail::f_g_asym (fb_s1);
      V s1         = s1_fixed - hp_v * in - v1n + v2n + fb_clipped;
      V s2
        = s2_fixed + hp_v * in + lp_v * in + v1n - (T) 2.0 * (v2n + fb_clipped);
      V f1  = -v1n + v1_v + h_v * s1;
      V f2  = -v2n + v2_v + h_v * s2;
      V res = vec_abs (f1) + vec_abs (f2);
      V a   = h_v * (kp1 * detail::f_dg_asym (fb_s1) - (T) 1.0) - (T) 1.0;
      V b   = h_v;
      V c   = h_v
        * ((T) 1.0 - (T) 2.0 * (k_v + (T) 1.0) * detail::f_dg_asym (fb_s1));
      V d    = (T) -2.0 * h_v - (T) 1.0;
      V norm = (T) 1.0 / (a * d - b * c);
      v1n -= (d * f1 - b * f2) * norm;
      v2n -= (a * f2 - c * f1) * norm;

      auto stillerror = res > (T) epsilon;
      if (vec_is_all_zeros (stillerror)) {
        break;
      }
    }
    x_v  = in;
    v1_v = v1n;
    v2_v = v2n;

    vec_store (&st[x * traits.size], x_v);
    vec_store (&st[v1 * traits.size], v1_v);
    vec_store (&st[v2 * traits.size], v2_v);
    return (v1n + hp_v * in);
  }
  //----------------------------------------------------------------------------
};
//##############################################################################
// Moog.
//##############################################################################
struct moog_1 {
  // From Saike's FM Filter2 (Yutani).
  // SHA b1c1952c5f798d968ff15714bca6a5c694e06d82
  // https://github.com/JoepVanlier/JSFX
  //----------------------------------------------------------------------------
  enum coeffs {
    k,
    q0s,
    r1s,
    k0s,
    k0g,
    rg1,
    rg2,
    rg3,
    rg4,
    qg1,
    qg2,
    qg3,
    qg4,
    choice,
    frac,
    n_coeffs
  };
  enum state {
    sf1,
    sf2,
    sf3,
    sf4,
    sg1,
    sg2,
    sg3,
    sg4,
    si1,
    si2,
    si3,
    si4,
    n_states
  };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    reset_coeffs (co, freq, reso, 0, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    bandpass_tag)
  {
    reset_coeffs (co, freq, reso, 1, vec_set<V> (0.), sr);
  }

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    reset_coeffs (co, freq, reso, 2, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    notch_tag)
  {
    reset_coeffs (co, freq, reso, 3, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V k_v      = vec_load<V> (&co[k * traits.size]);
    V q0s_v    = vec_load<V> (&co[q0s * traits.size]);
    V r1s_v    = vec_load<V> (&co[r1s * traits.size]);
    V k0s_v    = vec_load<V> (&co[k0s * traits.size]);
    V k0g_v    = vec_load<V> (&co[k0g * traits.size]);
    V rg1_v    = vec_load<V> (&co[rg1 * traits.size]);
    V rg2_v    = vec_load<V> (&co[rg2 * traits.size]);
    V rg3_v    = vec_load<V> (&co[rg3 * traits.size]);
    V rg4_v    = vec_load<V> (&co[rg4 * traits.size]);
    V qg1_v    = vec_load<V> (&co[qg1 * traits.size]);
    V qg2_v    = vec_load<V> (&co[qg2 * traits.size]);
    V qg3_v    = vec_load<V> (&co[qg3 * traits.size]);
    V qg4_v    = vec_load<V> (&co[qg4 * traits.size]);
    V frac_v   = vec_load<V> (&co[frac * traits.size]);
    V choice_v = vec_load<V> (&co[choice * traits.size]);
    V sf1_v    = vec_load<V> (&st[sf1 * traits.size]);
    V sf2_v    = vec_load<V> (&st[sf2 * traits.size]);
    V sf3_v    = vec_load<V> (&st[sf3 * traits.size]);
    V sf4_v    = vec_load<V> (&st[sf4 * traits.size]);
    V sg1_v    = vec_load<V> (&st[sg1 * traits.size]);
    V sg2_v    = vec_load<V> (&st[sg2 * traits.size]);
    V sg3_v    = vec_load<V> (&st[sg3 * traits.size]);
    V sg4_v    = vec_load<V> (&st[sg4 * traits.size]);
    V si1_v    = vec_load<V> (&st[si1 * traits.size]);
    V si2_v    = vec_load<V> (&st[si2 * traits.size]);
    V si3_v    = vec_load<V> (&st[si3 * traits.size]);
    V si4_v    = vec_load<V> (&st[si4 * traits.size]);

    V yo  = vec_tanh_approx_vaneev (k0g_v * (in + sg1_v));
    V a   = yo;
    V yi  = yo;
    V yd  = k0s_v * (yi + sf1_v);
    V y   = yd + si1_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    V b   = yo;
    si1_v = yd + y;
    sf1_v = r1s_v * yi - q0s_v * yo;
    yi    = yo;
    yd    = k0s_v * (yi + sf2_v);
    y     = yd + si2_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    V c   = yo;
    si2_v = yd + y;
    sf2_v = r1s_v * yi - q0s_v * yo;
    yi    = yo;
    yd    = k0s_v * (yi + sf3_v);
    y     = yd + si3_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    V d   = yo;
    si3_v = yd + y;
    sf3_v = r1s_v * yi - q0s_v * yo;
    yi    = yo;
    yd    = k0s_v * (yi + sf4_v);
    y     = yd + si4_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    si4_v = yd + y;
    sf4_v = r1s_v * yi - q0s_v * yo;
    V yf  = k_v * y;
    sg1_v = rg1_v * in + qg1_v * yf + sg2_v;
    sg2_v = rg2_v * in + qg2_v * yf + sg3_v;
    sg3_v = rg3_v * in + qg3_v * yf + sg4_v;
    sg4_v = rg4_v * in + qg4_v * yf;

    V f1, f2;
    switch ((uint) choice_v[0]) {
    case 0:
      f1 = y * ((T) 1. + k_v);
      f2 = (T) vt2 * ((T) 2. * c - (T) 2. * b);
      break;
    case 1:
      f1 = (T) vt2 * ((T) 2. * c - (T) 2. * b);
      f2 = (T) vt2 * (a - (T) 4. * b + (T) 6. * c - (T) 4. * d + yo);
      break;
    case 2:
      frac_v *= frac_v;
      frac_v *= frac_v;
      f1 = (T) vt2 * (a - (T) 4. * b + (T) 6. * c - (T) 4. * d + yo);
      f2 = (T) vt2 * (a - (T) 4. * b + (T) 6. * c - (T) 4. * d);
      break;
    case 3:
      f1 = (T) vt2 * (a - (T) 4. * b + (T) 6. * c - (T) 4. * d);
      f2 = y * ((T) 1. + k_v);
      break;
    default:
      f1 = vec_set<V> ((T) 0);
      f2 = vec_set<V> ((T) 0);
      break;
    }

    vec_store (&st[sf1 * traits.size], sf1_v);
    vec_store (&st[sf2 * traits.size], sf2_v);
    vec_store (&st[sf3 * traits.size], sf3_v);
    vec_store (&st[sf4 * traits.size], sf4_v);
    vec_store (&st[sg1 * traits.size], sg1_v);
    vec_store (&st[sg2 * traits.size], sg2_v);
    vec_store (&st[sg3 * traits.size], sg3_v);
    vec_store (&st[sg4 * traits.size], sg4_v);
    vec_store (&st[si1 * traits.size], si1_v);
    vec_store (&st[si2 * traits.size], si2_v);
    vec_store (&st[si3 * traits.size], si3_v);
    vec_store (&st[si4 * traits.size], si4_v);

    return (f2 * frac_v + f1 * (1.0 - frac_v)) * (T) vt2i;
  }
  //----------------------------------------------------------------------------
private:
  static constexpr double vt2  = 0.052;
  static constexpr double vt2i = 19.23076923076923;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    uint                        filter_type, // choice on the original code
    V                           fraction, // frac on the original code
    vec_value_type_t<V>         sr)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.

    // needs generous sample/rate oversampling
    freq = vec_min (sr < (T) 176200. ? (T) 12000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 96000. ? (T) 11000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 88200. ? (T) 6000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 48000. ? (T) 5500. : (T) 20000., freq);

    V fc = freq;
    T fs = sr;

    V k_v = reso * (T) 3.9999999999999987;
    vec_store (&co[k * traits.size], k_v);

    V g = vec_tan ((T) M_PI * fc / fs)
      / vec_sqrt (
            (T) 1.0 + vec_sqrt (k_v)
            - (T) 2. * vec_pow (k_v, (T) 0.25) * (T) M_SQRT1_2);

    V p0s = 1.0 / (1.0 + g);
    vec_store (&co[q0s * traits.size], 1.0 - g); // save memory?
    vec_store (&co[r1s * traits.size], -g);
    vec_store (&co[k0s * traits.size], (T) vt2 * g * p0s);

    V nmp = (1.0 - p0s);
    V gn  = nmp * nmp * nmp * nmp;
    V kgn = k_v * gn;
    V p0g = 1.0 / (1.0 + kgn);

    vec_store (&co[k0g * traits.size], -(T) vt2i * p0g);
    vec_store (&co[rg1 * traits.size], -4.0 * kgn); // save memory?
    vec_store (&co[rg2 * traits.size], -6.0 * kgn); // save memory?
    vec_store (&co[rg3 * traits.size], -4.0 * kgn); // save memory?
    vec_store (&co[rg4 * traits.size], -1.0 * kgn); // save memory?

    V tmp = p0s * (g - 1.0);
    V acc = tmp;

    vec_store (&co[qg1 * traits.size], -4.0 * (kgn + acc));
    acc *= tmp;
    vec_store (&co[qg2 * traits.size], -6.0 * (kgn + acc));
    acc *= tmp;
    vec_store (&co[qg3 * traits.size], -4.0 * (kgn + acc));
    acc *= tmp;
    vec_store (&co[qg4 * traits.size], -1.0 * (kgn + acc));

    // this integer shouldn't be part of a filter, it will delay when smoothing
    // adding 0.5 so it can effectively reach its value when casted to uint.
    auto cho = vec_set<V> ((T) filter_type + (T) 0.5);
    vec_store (&co[choice * traits.size], cho);
    vec_store (&co[frac * traits.size], fraction);
  }
};
//##############################################################################
// Moog2. (ladder)
//##############################################################################
struct moog_2 {
  // From Saike's FM Filter2 (Yutani).
  // SHA b1c1952c5f798d968ff15714bca6a5c694e06d82
  // https://github.com/JoepVanlier/JSFX
  //----------------------------------------------------------------------------
  enum coeffs {
    k,
    q0s,
    r1s,
    k0s,
    k0g,
    rg1,
    rg2,
    qg1,
    qg2,
    choice,
    frac,
    n_coeffs
  };
  enum state { sf1, sf2, sg1, sg2, si1, si2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    reset_coeffs (co, freq, reso, 0, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    bandpass_tag)
  {
    reset_coeffs (co, freq, reso, 1, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    reset_coeffs (co, freq, reso, 2, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    vec_value_type_t<V>         sr,
    notch_tag)
  {
    reset_coeffs (co, freq, reso, 3, vec_set<V> (0.), sr);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states;
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));
    assert (co.size() >= (n_coeffs * traits.size));

    V k_v      = vec_load<V> (&co[k * traits.size]);
    V q0s_v    = vec_load<V> (&co[q0s * traits.size]);
    V r1s_v    = vec_load<V> (&co[r1s * traits.size]);
    V k0s_v    = vec_load<V> (&co[k0s * traits.size]);
    V k0g_v    = vec_load<V> (&co[k0g * traits.size]);
    V rg1_v    = vec_load<V> (&co[rg1 * traits.size]);
    V rg2_v    = vec_load<V> (&co[rg2 * traits.size]);
    V qg1_v    = vec_load<V> (&co[qg1 * traits.size]);
    V qg2_v    = vec_load<V> (&co[qg2 * traits.size]);
    V frac_v   = vec_load<V> (&co[frac * traits.size]);
    V choice_v = vec_load<V> (&co[choice * traits.size]);
    V sf1_v    = vec_load<V> (&st[sf1 * traits.size]);
    V sf2_v    = vec_load<V> (&st[sf2 * traits.size]);
    V sg1_v    = vec_load<V> (&st[sg1 * traits.size]);
    V sg2_v    = vec_load<V> (&st[sg2 * traits.size]);
    V si1_v    = vec_load<V> (&st[si1 * traits.size]);
    V si2_v    = vec_load<V> (&st[si2 * traits.size]);

    V yo  = vec_tanh_approx_vaneev (k0g_v * (in + sg1_v));
    V a   = yo;
    V yi  = yo;
    V yd  = k0s_v * (yi + sf1_v);
    V y   = yd + si1_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    V b   = yo;
    si1_v = yd + y;
    sf1_v = r1s_v * yi - q0s_v * yo;
    yi    = yo;
    yd    = k0s_v * (yi + sf2_v);
    y     = yd + si2_v;
    yo    = vec_tanh_approx_vaneev ((T) vt2i * y);
    si2_v = yd + y;
    sf2_v = r1s_v * yi - q0s_v * yo;
    V yf  = k_v * y;
    sg1_v = rg1_v * in + qg1_v * yf + sg2_v;
    sg2_v = rg2_v * in + qg2_v * yf;

    V f1, f2;
    switch ((uint) choice_v[0]) {
    case 0:
      f1 = y * ((T) 1. + k_v);
      f2 = (T) vt2 * ((T) 2. * b - (T) 2. * yo) * (T) 8.;
      break;
    case 1:
      f1 = (T) vt2 * ((T) 2. * b - (T) 2. * yo) * (T) 8.;
      f2 = (T) vt2 * (a - b);
      break;
    case 2:
      frac_v *= frac_v;
      frac_v *= frac_v;
      f1 = (T) vt2 * (a - b);
      f2 = -(T) vt2 * ((T) 2. * b - (T) 2. * yo - a);
      break;
    case 3:
      f1 = -(T) vt2 * ((T) 2. * b - (T) 2. * yo - a);
      f2 = y * ((T) 1. + k_v);
      break;
    default:
      f1 = vec_set<V> ((T) 0);
      f2 = vec_set<V> ((T) 0);
      break;
    }

    vec_store (&st[sf1 * traits.size], sf1_v);
    vec_store (&st[sf2 * traits.size], sf2_v);
    vec_store (&st[sg1 * traits.size], sg1_v);
    vec_store (&st[sg2 * traits.size], sg2_v);
    vec_store (&st[si1 * traits.size], si1_v);
    vec_store (&st[si2 * traits.size], si2_v);

    return (f2 * frac_v + f1 * (1.0 - frac_v)) * (T) vt2i;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr double vt2  = 0.052;
  static constexpr double vt2i = 19.23076923076923;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           reso,
    uint                        filter_type, // choice on the original code
    V                           fraction, // frac on the original code
    vec_value_type_t<V>         sr)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.

    // needs generous sample/rate oversampling
    freq = vec_min (sr < (T) 176200. ? (T) 12000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 96000. ? (T) 11000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 88200. ? (T) 6000. : (T) 20000., freq);
    freq = vec_min (sr < (T) 48000. ? (T) 5500. : (T) 20000., freq);

    V fc = freq;
    T fs = sr;

    V k_v   = reso * 120.;
    V g     = vec_tan (3.141592653589793 / fs * fc) / vec_sqrt (1. + k_v);
    V p0s   = 1.0 / (1.0 + g);
    V q0s_v = 1.0 - g;
    V r1s_v = -g;
    V k0s_v = (T) vt2 * g * p0s;
    V nmp   = (1.0 - p0s);
    V gn    = vec_pow (nmp, 2.);
    V kgn   = k_v * gn;
    V p0g   = 1.0 / (1.0 + kgn);
    V k0g_v = -(T) vt2i * p0g;
    V rg1_v = -2.0 * kgn;
    V rg2_v = -1.0 * kgn;
    V tmp   = p0s * (g - 1.0);
    V acc   = tmp;
    V qg1_v = -2.0 * (kgn + acc);
    acc *= tmp;
    V qg2_v = -1.0 * (kgn + acc);

    vec_store (&co[k * traits.size], k_v);
    vec_store (&co[q0s * traits.size], q0s_v);
    vec_store (&co[r1s * traits.size], r1s_v);
    vec_store (&co[k0s * traits.size], k0s_v);
    vec_store (&co[k0g * traits.size], k0g_v);
    vec_store (&co[rg1 * traits.size], rg1_v);
    vec_store (&co[rg2 * traits.size], rg2_v);
    vec_store (&co[qg1 * traits.size], qg1_v);
    vec_store (&co[qg2 * traits.size], qg2_v);
    // this integer shouldn't be part of a filter, it will delay when smoothing
    // adding 0.5 so it can effectively reach its value when casted to uint.
    auto cho = vec_set<V> ((T) filter_type + (T) 0.5);
    vec_store (&co[choice * traits.size], cho);
    vec_store (&co[frac * traits.size], fraction);
  }
};
//------------------------------------------------------------------------------
}} // namespace artv::saike
