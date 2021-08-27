#pragma once

// Port of Presence e.q by liteon

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

//------------------------------------------------------------------------------
namespace artv { namespace liteon {
// Presence (moorer) JSFX
struct presence_high_shelf {
  //----------------------------------------------------------------------------
  enum coeffs { a0, a1, a2, b1, b2, n_coeffs };
  enum state { y1, y2, x1, x2, n_states };
  //----------------------------------------------------------------------------
  // BW on the original JSFX is unitless BW from 0.007 to 0.4 and the lower end
  // makes it narrower (?). I scale it from 0 to 1.

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    V                           bogus_q,
    V                           gain_db,
    vec_value_type_t<V>         sr)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V a {};
    V a2plus1 {};
    V alphad {};
    V alphan {};
    V as2 {};
    V as4 {};
    V asnd {};
    V b0 {};
    V boost {};
    V bw {};
    V c {};
    V ca {};
    V cf {};
    V cs {};
    V d {};
    V delta {};
    V f {};
    V f2 {};
    V ma2plus1 {};
    V mag {};
    V recipb0 {};
    V sn {};
    V t {};
    V theta {};
    V tmp {};
    V xfmbw {};

    freq = vec_max (freq, 3100.);
    freq = vec_min (freq, 18500.);

    cf    = freq / sr;
    boost = gain_db;

    // scaling to the orignal range
    for (uint i = 0; i < traits.size; ++i) {
      assert (bogus_q[i] <= 1.);
    }

    bw = (((T) 1. - bogus_q) * (T) 0.33) + (T) 0.07;

    ca = vec_tan ((T) M_PI * (cf - (T) 0.25));
    a  = vec_exp ((boost / (T) 20.) * (T) M_LN10);
#if 0
    // original code, kept here...
    if ((boost < 6.0) && (boost > -6.0)) {
      (f = std::sqrt (a));
    }
    else {
      if (a > 1.) {
        (f = a / std::sqrt (2.));
      }
      else {
        (f = a * std::sqrt (2.));
      }
    }
#else
    auto run_sqrt = vec_abs (boost) < (T) 6.;
    V    fsqrt {};
    if (!vec_is_all_zeros (run_sqrt)) {
      fsqrt = vec_sqrt (a);
    }
    V fsqrt2 = (a > (T) 1.) ? vec_set<V> ((T) (1. / M_SQRT2))
                            : vec_set<V> ((T) M_SQRT2);
    fsqrt2 *= a;
    f = run_sqrt ? fsqrt : fsqrt2;
#endif
    t     = vec_tan ((T) (M_PI * 2.) * bw);
    as2   = ca * ca;
    as4   = as2 * as2;
    d     = 2. * as2 * t;
    sn    = (1. + as4) * t;
    cs    = (1. - as4);
    mag   = vec_sqrt (sn * sn + cs * cs);
    d     = mag;
    delta = vec_atan2 (sn, cs);

    if (!vec_is_all_zeros (vec_abs (d) <= (T) 1.)) {
      asnd = vec_asin (d);
    }
    theta = (T) 0.5 * ((T) M_PI - asnd - delta);
    tmp   = (T) 0.5 * (asnd - delta);
    theta = (tmp > (T) 0. && tmp < theta) ? tmp : theta;
    xfmbw = theta * (T) (1. / (2. * M_PI));
    c     = 1. / vec_tan ((T) (M_PI * 2.) * xfmbw);
    f2    = f * f;
    tmp   = a * a - f2;
#if 0
// On the original script SPN was undefined...
    if (std::fabs (tmp) <= spn) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - (T)1.) / tmp));
    }
#else
#if 0
    if (std::fabs (tmp) <= (T)0.) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - (T)1.) / tmp));
    }
#else
    alphad = (tmp == (T) 0.) ? c : vec_sqrt (c * c * (f2 - (T) 1.) / tmp);
#endif
#endif
    alphan   = a * alphad;
    a2plus1  = 1. + as2;
    ma2plus1 = 1. - as2;

    V a0_v, a1_v, a2_v, b1_v, b2_v;
    a0_v    = a2plus1 + alphan * ma2plus1;
    a1_v    = 4.0 * ca;
    a2_v    = a2plus1 - alphan * ma2plus1;
    b0      = a2plus1 + alphad * ma2plus1;
    b2_v    = a2plus1 - alphad * ma2plus1;
    recipb0 = 1. / b0;
    a0_v *= recipb0;
    a1_v *= recipb0;
    a2_v *= recipb0;
    b1_v = a1_v;
    b2_v *= recipb0;
    b1_v = -b1_v;
    b2_v = -b2_v;

    vec_store (&co[a0 * traits.size], a0_v);
    vec_store (&co[a1 * traits.size], a1_v);
    vec_store (&co[a2 * traits.size], a2_v);
    vec_store (&co[b1 * traits.size], b1_v);
    vec_store (&co[b2 * traits.size], b2_v);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs));
    assert (st.size() >= (traits.size * n_states));

    V a0v = vec_set<V> (co[a0]);
    V a1v = vec_set<V> (co[a1]);
    V a2v = vec_set<V> (co[a2]);
    V b1v = vec_set<V> (co[b1]);
    V b2v = vec_set<V> (co[b2]);

    V x1v = vec_load<V> (&st[x1 * traits.size]);
    V x2v = vec_load<V> (&st[x2 * traits.size]);
    V y1v = vec_load<V> (&st[y1 * traits.size]);
    V y2v = vec_load<V> (&st[y2 * traits.size]);

    V out = in;
    out *= a0v;

    // TODO: worth all that load and storing for this (?)
    out += a1v * x1v + a2v * x2v + b1v * y1v + b2v * y2v;

    vec_store (&st[x2 * traits.size], x1v);
    vec_store (&st[x1 * traits.size], in);
    vec_store (&st[y2 * traits.size], y1v);
    vec_store (&st[y1 * traits.size], out);

    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * n_states));

    V a0v = vec_load<V> (&co[a0 * traits.size]);
    V a1v = vec_load<V> (&co[a1 * traits.size]);
    V a2v = vec_load<V> (&co[a2 * traits.size]);
    V b1v = vec_load<V> (&co[b1 * traits.size]);
    V b2v = vec_load<V> (&co[b2 * traits.size]);

    V x1v = vec_load<V> (&st[x1 * traits.size]);
    V x2v = vec_load<V> (&st[x2 * traits.size]);
    V y1v = vec_load<V> (&st[y1 * traits.size]);
    V y2v = vec_load<V> (&st[y2 * traits.size]);

    V out = in;
    out *= a0v;

    // TODO: worth all that load and storing for this (?)
    out += a1v * x1v + a2v * x2v + b1v * y1v + b2v * y2v;

    vec_store (&st[x2 * traits.size], x1v);
    vec_store (&st[x1 * traits.size], in);
    vec_store (&st[y2 * traits.size], y1v);
    vec_store (&st[y1 * traits.size], out);

    return out;
  }
  //----------------------------------------------------------------------------
};

}} // namespace artv::liteon
