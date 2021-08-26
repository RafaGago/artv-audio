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
  static void high_shelf (
    crange<double> co,
    double         freq,
    double         bogus_q,
    double         gain_db,
    double         sr)
  {
    double a        = 0.;
    double a2plus1  = 0.;
    double alphad   = 0.;
    double alphan   = 0.;
    double as2      = 0.;
    double as4      = 0.;
    double asnd     = 0.;
    double b0       = 0.;
    double boost    = 0.;
    double bw       = 0.;
    double c        = 0.;
    double ca       = 0.;
    double cf       = 0.;
    double cs       = 0.;
    double d        = 0.;
    double delta    = 0.;
    double f        = 0.;
    double f2       = 0.;
    double ma2plus1 = 0.;
    double mag      = 0.;
    double recipb0  = 0.;
    double sn       = 0.;
    double t        = 0.;
    double theta    = 0.;
    double tmp      = 0.;
    double xfmbw    = 0.;

    freq = std::max (freq, 3100.);
    freq = std::min (freq, 18500.);

    cf    = freq / sr;
    boost = gain_db;

    // scaling to the orignal range
    assert (bogus_q <= 1.);
    bw = ((1. - bogus_q) * 0.33) + 0.07;

    ca = std::tan (3.141592653589793 * (cf - 0.25));
    a  = pow (10., (boost / 20.));
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
    t     = std::tan (2. * 3.141592653589793 * bw);
    as2   = ca * ca;
    as4   = as2 * as2;
    d     = 2. * as2 * t;
    sn    = (1. + as4) * t;
    cs    = (1. - as4);
    mag   = std::sqrt (sn * sn + cs * cs);
    d     = mag;
    delta = std::atan2 (sn, cs);
    asnd  = (d <= 1. && d >= -1.) ? asin (d) : 0;
    theta = 0.5 * (3.141592653589793 - asnd - delta);
    tmp   = 0.5 * (asnd - delta);
    if (tmp > 0. && tmp < theta) {
      (theta = tmp);
    }
    xfmbw = theta / (2. * 3.141592653589793);
    c     = 1. / std::tan (2. * 3.141592653589793 * xfmbw);
    f2    = f * f;
    tmp   = a * a - f2;
#if 0
    if (std::fabs (tmp) <= spn) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#else
    if (std::fabs (tmp) <= 0.) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#endif
    alphan   = a * alphad;
    a2plus1  = 1. + as2;
    ma2plus1 = 1. - as2;
    co[a0]   = a2plus1 + alphan * ma2plus1;
    co[a1]   = 4.0 * ca;
    co[a2]   = a2plus1 - alphan * ma2plus1;
    b0       = a2plus1 + alphad * ma2plus1;
    co[b2]   = a2plus1 - alphad * ma2plus1;
    recipb0  = 1. / b0;
    co[a0] *= recipb0;
    co[a1] *= recipb0;
    co[a2] *= recipb0;
    co[b1] = co[a1];
    co[b2] *= recipb0;
    co[b1] = -co[b1];
    co[b2] = -co[b2];
    ;
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void reset_states (crange<double> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double out = co[a0] * in + co[a1] * st[x1] + co[a2] * st[x2]
      + co[b1] * st[y1] + co[b2] * st[y2];
    st[x2] = st[x1];
    st[x1] = in;
    st[y2] = st[y1];
    st[y1] = out;
    return out;
  }
  //----------------------------------------------------------------------------
  static double_x2 tick (
    crange<const double>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    double_x2                     in)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    double_x2 a1v = vec_set<2> (co[a1]);
    double_x2 a2v = vec_set<2> (co[a2]);
    double_x2 b1v = vec_set<2> (co[b1]);
    double_x2 b2v = vec_set<2> (co[b2]);

    double_x2 x1v = {st[0][x1], st[1][x1]};
    double_x2 x2v = {st[0][x2], st[1][x2]};
    double_x2 y1v = {st[0][y1], st[1][y1]};
    double_x2 y2v = {st[0][y2], st[1][y2]};

    double_x2 out = {in[0], in[1]};
    out *= co[a0];

    // TODO: worth all that load and storing for this (?)
    out += a1v * x1v + a2v * x2v + b1v * y1v + b2v * y2v;

    for (uint i = 0; i < 2; ++i) {
      st[i][x2] = st[i][x1];
      st[i][x1] = in[i];
      st[i][y2] = st[i][y1];
      st[i][y1] = out[i];
    }
    return out;
  }
};

}} // namespace artv::liteon
