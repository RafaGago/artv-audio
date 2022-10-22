#pragma once

// Port of Presence e.q by liteon

#include <cmath>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

//------------------------------------------------------------------------------
namespace artv { namespace liteon {
// Presence (moorer) JSFX
struct presence_high_shelf {
  //----------------------------------------------------------------------------
  enum coeffs { a0, a1, a2, b1, b2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, y2, x1, x2, n_states };
  //----------------------------------------------------------------------------
  // BW on the original JSFX is unitless BW from 0.007 to 0.4 and the lower end
  // makes it narrower (?). I scale it from 0 to 1.
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   freq,
    V                   bogus_q,
    V                   gain_db,
    vec_value_type_t<V> t_spl)
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

    cf    = freq * t_spl;
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

    co[a0] = a0_v;
    co[a1] = a1_v;
    co[a2] = a2_v;
    co[b1] = b1_v;
    co[b2] = b2_v;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<vec_value_type_t<V> const> co, // coeffs (1 set)
    xspan<V>                         st, // states (interleaved, SIMD aligned)
    V                                in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V a0v = vec_set<V> (co[a0]);
    V a1v = vec_set<V> (co[a1]);
    V a2v = vec_set<V> (co[a2]);
    V b1v = vec_set<V> (co[b1]);
    V b2v = vec_set<V> (co[b2]);

    V out = in;
    out *= a0v;
    out += a1v * st[x1] + a2v * st[x2] + b1v * st[y1] + b2v * st[y2];

    st[x2] = st[x1];
    st[x1] = in;
    st[y2] = st[y1];
    st[y1] = out;

    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    xspan<V const> co, // coeffs (interleaved, SIMD aligned)
    xspan<V>       st, // states (interleaved, SIMD aligned)
    V              in)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V out = in;
    out *= co[a0];
    out
      += co[a1] * st[x1] + co[a2] * st[x2] + co[b1] * st[y1] + co[b2] * st[y2];

    st[x2] = st[x1];
    st[x1] = in;
    st[y2] = st[y1];
    st[y1] = out;

    return out;
  }
  //----------------------------------------------------------------------------
};

}} // namespace artv::liteon
