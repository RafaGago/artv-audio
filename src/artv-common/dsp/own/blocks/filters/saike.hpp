#pragma once

// Saike's excellent filters.

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"
//------------------------------------------------------------------------------
namespace artv { namespace saike {

namespace detail {

inline double f_g (double v)
{
  return std::max (-1., std::min (1., v));
}
inline double f_dg (double v)
{
  return 1. - 1. * (abs (v) > 1.);
}

inline double f_g_asym (double v)
{
  return v > 0 ? std::min (1., v) : std::max (-1., v * .25);
}
inline double f_dg_asym (double v)
{
  return v > 0 ? 1. - 1. * (std::abs (v) > 1) : .25 - .25 * (std::abs (v) > 4.);
}
//##############################################################################
// MS20 /K35: base
//##############################################################################
struct ms20_base {
  // From Saike's FM Filter2 (Yutani).
  // SHA b1c1952c5f798d968ff15714bca6a5c694e06d82
  // https://github.com/JoepVanlier/JSFX
  //----------------------------------------------------------------------------
  enum coeffs { h, hh, k, n_coeffs };
  enum state { y1, y2, d1, d2, n_states };
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
protected:
  static void init (crange<double> c, double freq, double reso, double sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.

    assert (c.size() >= n_coeffs);
    freq                = std ::min (sr < 80000. ? 14000. : 20000., freq);
    double oversampling = (sr / 44100.);
    double norm_freq    = (freq / 20000.);
    double f0           = norm_freq * M_PI / oversampling;
    c[h]                = tan (f0 / (2.1 * oversampling)) * 2.1 * oversampling;
    c[hh]               = 0.5 * c[h]; // is memory needed for this?
    c[k]                = 2 * reso;
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
  static void lowpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double gd2k      = detail::f_g (st[d2] * co[k]);
    double tanhterm1 = std::tanh (-st[d1] + in - gd2k);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2          = co[k] * st[y2];
      double gky2         = detail::f_g (ky2);
      double dgky2        = detail::f_dg (ky2);
      double sig1         = in - st[y1] - gky2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = -hhthsig1sqm1 + 1.;
      double b            = -co[k] * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (co[k] * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[d2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_highpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void highpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double kc        = .85 * co[k];
    double gkd2px    = detail::f_g (kc * (st[d2] + in));
    double tanhterm1 = std::tanh (-st[d1] - gkd2px);
    double tanhterm2 = std::tanh (st[d1] - st[d2] - in + gkd2px);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double kxpy2        = kc * (in + st[y2]);
      double gkxpy2       = detail::f_g (kxpy2);
      double dgky2px      = detail::f_dg (kxpy2);
      double sig1         = -st[y1] - gkxpy2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = -in + st[y1] - st[y2] + gkxpy2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = (thsig1sq - 1.);
      double hhthsig2sqm1 = (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = -hhthsig1sqm1 + 1.;
      double b            = -kc * hhthsig1sqm1 * dgky2px;
      double c            = hhthsig2sqm1;
      double d            = (kc * dgky2px - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[y2] + in;
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_bandpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void bandpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double kc        = .95 * co[k];
    double gd2k      = detail::f_g (st[d2] * kc);
    double tanhterm1 = std::tanh (-st[d1] - in - gd2k);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + in + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2          = kc * st[y2];
      double gky2         = detail::f_g (ky2);
      double dgky2        = detail::f_dg (ky2);
      double sig1         = -in - st[y1] - gky2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = in + st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = 1. - hhthsig1sqm1;
      double b            = -kc * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (kc * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[d2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_notch : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void notch (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double gd2k      = detail::f_g (st[d2] * co[k]);
    double tanhterm1 = std::tanh (-st[d1] - in - gd2k);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + in + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2          = co[k] * st[y2];
      double gky2         = detail::f_g (ky2);
      double dgky2        = detail::f_dg (ky2);
      double sig1         = -in - st[y1] - gky2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = in + st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = 1. - hhthsig1sqm1;
      double b            = -co[k] * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (co[k] * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return in - st[y2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//##############################################################################
// MS20 /K35: "MS20_nonlin_asym_tanh" in the original sources.
//##############################################################################
//----------------------------------------------------------------------------
struct ms20_asym_lowpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double gd2k = detail::f_g_asym (st[d2] * co[k]);
    double sig1 = -st[d1] + in - gd2k;
    double qq   = (sig1 < 0.) ? .6 : 1.0;

    sig1 *= qq;
    double tanhterm1 = std::tanh (sig1);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2   = co[k] * st[y2];
      double gky2  = detail::f_g_asym (ky2);
      double dgky2 = detail::f_dg_asym (ky2);
      sig1         = in - st[y1] - gky2;
      qq           = (sig1 < 0.) ? .6 : 1.0;
      sig1 *= qq;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = -qq * hhthsig1sqm1 + 1.;
      double b            = -qq * co[k] * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (co[k] * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[d2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_highpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void highpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double kc        = co[k];
    double gkd2px    = detail::f_g_asym (kc * (st[d2] + in));
    double tanhterm1 = std::tanh (-st[d1] - gkd2px);
    double tanhterm2 = std::tanh (st[d1] - st[d2] - in + gkd2px);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double kxpy2        = kc * (in + st[y2]);
      double gkxpy2       = detail::f_g_asym (kxpy2);
      double dgky2px      = detail::f_dg_asym (kxpy2);
      double sig1         = -st[y1] - gkxpy2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = -in + st[y1] - st[y2] + gkxpy2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = (thsig1sq - 1.);
      double hhthsig2sqm1 = (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = -hhthsig1sqm1 + 1.;
      double b            = -kc * hhthsig1sqm1 * dgky2px;
      double c            = hhthsig2sqm1;
      double d            = (kc * dgky2px - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1.0 / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[y2] + in;
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_bandpass : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void bandpass (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double kc   = co[k];
    double gd2k = detail::f_g_asym (st[d2] * kc);
    double sig1 = -st[d1] - in - gd2k;
    double qq   = (sig1 < 0.) ? .6 : 1.0;
    sig1 *= qq;
    double tanhterm1 = std::tanh (sig1);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + in + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2   = kc * st[y2];
      double gky2  = detail::f_g_asym (ky2);
      double dgky2 = detail::f_dg_asym (ky2);
      sig1         = -in - st[y1] - gky2;
      qq           = (sig1 < 0.) ? .6 : 1.0;
      sig1 *= qq;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = in + st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = 1. - qq * hhthsig1sqm1;
      double b            = -qq * kc * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (kc * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return st[d2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
struct ms20_asym_notch : public detail::ms20_base {
  //----------------------------------------------------------------------------
  static void notch (crange<double> c, double freq, double reso, double sr)
  {
    init (c, freq, reso, sr);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double gd2k      = detail::f_g_asym (st[d2] * co[k]);
    double tanhterm1 = std::tanh (-st[d1] - in - gd2k);
    double tanhterm2 = std::tanh (st[d1] - st[d2] + in + gd2k);

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double ky2          = co[k] * st[y2];
      double gky2         = detail::f_g_asym (ky2);
      double dgky2        = detail::f_dg_asym (ky2);
      double sig1         = -in - st[y1] - gky2;
      double thsig1       = std::tanh (sig1);
      double thsig1sq     = thsig1 * thsig1;
      double sig2         = in + st[y1] - st[y2] + gky2;
      double thsig2       = std::tanh (sig2);
      double thsig2sq     = thsig2 * thsig2;
      double hhthsig1sqm1 = co[hh] * (thsig1sq - 1.);
      double hhthsig2sqm1 = co[hh] * (thsig2sq - 1.);
      double f1           = st[y1] - st[d1] - co[hh] * (tanhterm1 + thsig1);
      double f2           = st[y2] - st[d2] - co[hh] * (tanhterm2 + thsig2);
      res                 = std::fabs (f1) + std::fabs (f2);
      double a            = 1. - hhthsig1sqm1;
      double b            = -co[k] * hhthsig1sqm1 * dgky2;
      double c            = hhthsig2sqm1;
      double d            = (co[k] * dgky2 - 1.) * hhthsig2sqm1 + 1.;
      double norm         = 1. / (a * d - b * c);
      st[y1] -= (d * f1 - b * f2) * norm;
      st[y2] -= (a * f2 - c * f1) * norm;
    }
    st[d1] = st[y1];
    st[d2] = st[y2];
    return in - st[y2];
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
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
  static void lowpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0., sr);
  }

  static void bandpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0.25, sr);
  }

  static void highpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0.5, sr);
  }

  static void notch (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0.75, sr);
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
  static void init (
    crange<double> co,
    double         freq,
    double         reso,
    double         morph,
    double         sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.

    assert (co.size() >= n_coeffs);
    freq                  = std ::min (sr < 80000. ? 14000. : 20000., freq);
    double oversampling   = (sr / 44100.);
    double norm_freq      = (freq / 20000.);
    double f0             = norm_freq / oversampling;
    co[h]                 = tan (0.5 * M_PI * f0);
    co[k]                 = 3.98 * reso;
    co[hsq]               = co[h] * co[h]; // calculated instead?
    co[kh]                = co[k] * co[h]; // calculated instead?
    co[normalizing_const] = 1.0 / (-kh + hsq + 2. * h + 1.);

    /*alpha = 20.94153124476462;
    beta = 0.057872340425531923;
    vref = log(h / (alpha * beta)) / alpha;*/

    if (morph < 0.25) {
      co[lp] = 1.0 - morph * 4;
      co[bp] = morph * 4;
      co[hp] = 0;
    }
    else if (morph < 0.5) {
      co[hp] = (morph - 0.25) * 4;
      co[bp] = 1.0 - (morph - 0.25) * 4;
      co[lp] = 0;
    }
    else if (morph < 0.75) {
      co[hp] = 1.0;
      co[bp] = -2.0 * (morph - 0.5) * 4;
      co[lp] = (morph - 0.5) * 4;
    }
    else {
      co[hp] = 1.0 - (morph - 0.75) * 4;
      co[bp] = -2.0 * (1 - (morph - 0.75) * 4);
      co[lp] = 1;
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
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double x_xn = st[x] + in;
    double v1n  = (-co[kh] * st[v1] + (co[bp] - co[hp]) * co[h] * x_xn
                  + co[hsq] * ((co[lp] - co[hp]) * x_xn - st[v1])
                  + 2. * co[h] * st[v2] + st[v1])
      * co[normalizing_const];
    double v2n
      = (co[kh] * ((co[bp] - co[hp]) * x_xn + st[v2] - 2.0 * st[v1])
         + (co[lp] - co[bp]) * co[hsq] * x_xn + (co[lp] - co[bp]) * co[h] * x_xn
         - co[hsq] * st[v2] + st[v2])
      * co[normalizing_const];
    double kp1 = co[k] + 1.0;
    double fb  = detail::f_g ((co[k] + 1.) * (co[hp] * st[x] + st[v1]));
    double s1_fixed
      = 2.0 * co[bp] * st[x] - co[hp] * st[x] - st[v1] + st[v2] + fb;
    double s2_fixed = -4.0 * co[bp] * st[x] + co[hp] * st[x]
      - 2.0 * (st[v2] + fb) + st[v1] + co[lp] * st[x];

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double fb_s1      = kp1 * (co[hp] * in + v1n);
      double fb_clipped = detail::f_g (fb_s1);
      double s1         = s1_fixed - co[hp] * in - v1n + v2n + fb_clipped;
      double s2
        = s2_fixed + co[hp] * in + co[lp] * in + v1n - 2.0 * (v2n + fb_clipped);
      double f1   = -v1n + st[v1] + co[h] * s1;
      double f2   = -v2n + st[v2] + co[h] * s2;
      res         = std::fabs (f1) + std::fabs (f2);
      double a    = co[h] * (kp1 * detail::f_dg (fb_s1) - 1.0) - 1.0;
      double b    = co[h];
      double c    = co[h] * (1.0 - 2.0 * (co[k] + 1.0) * detail::f_dg (fb_s1));
      double d    = -2.0 * co[h] - 1.0;
      double norm = 1.0 / (a * d - b * c);
      v1n -= (d * f1 - b * f2) * norm;
      v2n -= (a * f2 - c * f1) * norm;
    }
    st[x]  = in;
    st[v1] = v1n;
    st[v2] = v2n;
    return (v1n + co[hp] * in);
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
  //----------------------------------------------------------------------------
};
//##############################################################################
// Steiner1: Crazy still!!! (eval_steiner_asym).
//##############################################################################
struct steiner_2 : public detail::steiner_base {
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double x_xn = st[x] + in;
    double v1n  = (-co[kh] * st[v1] + (co[bp] - co[hp]) * co[h] * x_xn
                  + co[hsq] * ((co[lp] - co[hp]) * x_xn - st[v1])
                  + 2. * co[h] * st[v2] + st[v1])
      * co[normalizing_const];
    double v2n
      = (co[kh] * ((co[bp] - co[hp]) * x_xn + st[v2] - 2.0 * st[v1])
         + (co[lp] - co[bp]) * co[hsq] * x_xn + (co[lp] - co[bp]) * co[h] * x_xn
         - co[hsq] * st[v2] + st[v2])
      * co[normalizing_const];
    double kp1 = co[k] + 1.0;
    double fb  = detail::f_g_asym ((co[k] + 1.) * (co[hp] * st[x] + st[v1]));
    double s1_fixed
      = 2.0 * co[bp] * st[x] - co[hp] * st[x] - st[v1] + st[v2] + fb;
    double s2_fixed = -4.0 * co[bp] * st[x] + co[hp] * st[x]
      - 2.0 * (st[v2] + fb) + st[v1] + co[lp] * st[x];

    double res = INFINITY;
    for (uint i = 0; i < max_iter && epsilon < res; ++i) {
      double fb_s1      = kp1 * (co[hp] * in + v1n);
      double fb_clipped = detail::f_g_asym (fb_s1);
      double s1         = s1_fixed - co[hp] * in - v1n + v2n + fb_clipped;
      double s2
        = s2_fixed + co[hp] * in + co[lp] * in + v1n - 2.0 * (v2n + fb_clipped);
      double f1 = -v1n + st[v1] + co[h] * s1;
      double f2 = -v2n + st[v2] + co[h] * s2;
      res       = std::fabs (f1) + std::fabs (f2);
      double a  = co[h] * (kp1 * detail::f_dg_asym (fb_s1) - 1.0) - 1.0;
      double b  = co[h];
      double c
        = co[h] * (1.0 - 2.0 * (co[k] + 1.0) * detail::f_dg_asym (fb_s1));
      double d    = -2.0 * co[h] - 1.0;
      double norm = 1.0 / (a * d - b * c);
      v1n -= (d * f1 - b * f2) * norm;
      v2n -= (a * f2 - c * f1) * norm;
    }
    st[x]  = in;
    st[v1] = v1n;
    st[v2] = v2n;
    return (v1n + co[hp] * in);
  }
  //----------------------------------------------------------------------------
  // TODO SIMD
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
    vt2,
    vt2i,
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
  static void lowpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void bandpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 1, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 2, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void notch (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 3, 0., sr);
  }
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double yo = std::tanh (co[k0g] * (in + st[sg1]));
    double a  = yo;
    double yi = yo;
    double yd = co[k0s] * (yi + st[sf1]);
    double y  = yd + st[si1];
    yo        = std::tanh (co[vt2i] * y);
    double b  = yo;
    st[si1]   = yd + y;
    st[sf1]   = co[r1s] * yi - co[q0s] * yo;
    yi        = yo;
    yd        = co[k0s] * (yi + st[sf2]);
    y         = yd + st[si2];
    yo        = std::tanh (co[vt2i] * y);
    double c  = yo;
    st[si2]   = yd + y;
    st[sf2]   = co[r1s] * yi - co[q0s] * yo;
    yi        = yo;
    yd        = co[k0s] * (yi + st[sf3]);
    y         = yd + st[si3];
    yo        = std::tanh (co[vt2i] * y);
    double d  = yo;
    st[si3]   = yd + y;
    st[sf3]   = co[r1s] * yi - co[q0s] * yo;
    yi        = yo;
    yd        = co[k0s] * (yi + st[sf4]);
    y         = yd + st[si4];
    yo        = std::tanh (co[vt2i] * y);
    st[si4]   = yd + y;
    st[sf4]   = co[r1s] * yi - co[q0s] * yo;
    double yf = co[k] * y;
    st[sg1]   = co[rg1] * in + co[qg1] * yf + st[sg2];
    st[sg2]   = co[rg2] * in + co[qg2] * yf + st[sg3];
    st[sg3]   = co[rg3] * in + co[qg3] * yf + st[sg4];
    st[sg4]   = co[rg4] * in + co[qg4] * yf;

    double f1, f2;
    double fracv = co[frac];
    switch (*(u64*) &co[choice]) {
    case 0:
      f1 = y * (1. + co[k]);
      f2 = vt2 * (2. * c - 2. * b);
      break;
    case 1:
      f1 = vt2 * (2. * c - 2. * b);
      f2 = vt2 * (a - 4. * b + 6. * c - 4. * d + yo);
      break;
    case 2:
      fracv *= fracv;
      fracv *= fracv;
      f1 = vt2 * (a - 4. * b + 6. * c - 4. * d + yo);
      f2 = vt2 * (a - 4. * b + 6. * c - 4. * d);
      break;
    case 3:
      f1 = vt2 * (a - 4. * b + 6. * c - 4. * d);
      f2 = y * (1. + co[k]);
      break;
    default:
      f1 = 0;
      f2 = 0;
      break;
    }
    return (f2 * fracv + f1 * (1.0 - fracv)) * co[vt2i];
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (
    crange<double>       co_smooth,
    crange<const double> co)
  {
    assert (co_smooth.size() >= n_coeffs);
    assert (co.size() >= n_coeffs);
    // the filter type as integer can't be smoothed
    co_smooth[choice] = co[choice];
  }
  //----------------------------------------------------------------------------
  static void reset_states (crange<double> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static void init (
    crange<double> co,
    double         freq,
    double         reso,
    uint           filter_type, // choice on the original code
    double         fraction, // frac on the original code
    double         sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.
    assert (co.size() >= n_coeffs);

    // needs generous sample/rate oversampling
    freq = std ::min (sr < 176200. ? 12000. : 20000., freq);
    freq = std ::min (sr < 96000. ? 11000. : 20000., freq);
    freq = std ::min (sr < 88200. ? 6000. : 20000., freq);
    freq = std ::min (sr < 48000. ? 5500. : 20000., freq);

    double norm_freq = (freq / 20000.);
    double fc        = .5 * sr * norm_freq;
    double fs        = sr;

    co[k]    = reso * 3.9999999999999987;
    double g = std::tan (3.141592653589793 * fc / fs)
      / std::sqrt (
                 1.0 + std::sqrt (co[k])
                 - 2. * std::pow (co[k], 0.25) * 0.7071067811865476);
    co[vt2]    = 0.052;
    co[vt2i]   = 19.23076923076923;
    double p0s = 1.0 / (1.0 + g);
    co[q0s]    = 1.0 - g;
    co[r1s]    = -g;
    co[k0s]    = co[vt2] * g * p0s;
    double nmp = (1.0 - p0s);
    double gn  = nmp * nmp * nmp * nmp;
    double kgn = co[k] * gn;
    double p0g = 1.0 / (1.0 + kgn);
    co[k0g]    = -co[vt2i] * p0g;
    co[rg1]    = -4.0 * kgn;
    co[rg2]    = -6.0 * kgn;
    co[rg3]    = -4.0 * kgn;
    co[rg4]    = -1.0 * kgn;
    double tmp = p0s * (g - 1.0);
    double acc = tmp;
    co[qg1]    = -4.0 * (kgn + acc);
    acc *= tmp;
    co[qg2] = -6.0 * (kgn + acc);
    acc *= tmp;
    co[qg3] = -4.0 * (kgn + acc);
    acc *= tmp;
    co[qg4]               = -1.0 * (kgn + acc);
    *((u64*) &co[choice]) = filter_type;
    co[frac]              = fraction;
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
    vt2,
    vt2i,
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
  static void lowpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 0, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void bandpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 1, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 2, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void notch (crange<double> co, double freq, double reso, double sr)
  {
    init (co, freq, reso, 3, 0., sr);
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (
    crange<double>       co_smooth,
    crange<const double> co)
  {
    assert (co_smooth.size() >= n_coeffs);
    assert (co.size() >= n_coeffs);
    // the filter type as integer can't be smoothed
    co_smooth[choice] = co[choice];
  }
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

    double yo = std::tanh (co[k0g] * (in + st[sg1]));
    double a  = yo;
    double yi = yo;
    double yd = co[k0s] * (yi + st[sf1]);
    double y  = yd + st[si1];
    yo        = std::tanh (co[vt2i] * y);
    double b  = yo;
    st[si1]   = yd + y;
    st[sf1]   = co[r1s] * yi - co[q0s] * yo;
    yi        = yo;
    yd        = co[k0s] * (yi + st[sf2]);
    y         = yd + st[si2];
    yo        = std::tanh (co[vt2i] * y);
    st[si2]   = yd + y;
    st[sf2]   = co[r1s] * yi - co[q0s] * yo;
    double yf = co[k] * y;
    st[sg1]   = co[rg1] * in + co[qg1] * yf + st[sg2];
    st[sg2]   = co[rg2] * in + co[qg2] * yf;

    double f1, f2;
    double fracv = co[frac];
    switch (*((u64*) &co[choice])) {
    case 0:
      f1 = y * (1. + co[k]);
      f2 = vt2 * (2. * b - 2. * yo) * 8.;
      break;
    case 1:
      f1 = vt2 * (2. * b - 2. * yo) * 8.;
      f2 = vt2 * (a - b);
      break;
    case 2:
      fracv *= fracv;
      fracv *= fracv;
      f1 = vt2 * (a - b);
      f2 = -vt2 * (2. * b - 2. * yo - a);
      break;
    case 3:
      f1 = -vt2 * (2. * b - 2. * yo - a);
      f2 = y * (1. + co[k]);
      break;
    default:
      f1 = 0;
      f2 = 0;
      break;
    }
    return (f2 * fracv + f1 * (1.0 - fracv)) * co[vt2i];
  }
  //----------------------------------------------------------------------------
  static void init (
    crange<double> co,
    double         freq,
    double         reso,
    uint           filter_type, // choice on the original code
    double         fraction, // frac on the original code
    double         sr)
  {
    // The original filter freq seems to be off by 10% 22000 vs 20000kHz.
    // It also doesn't respect the frequency.
    assert (co.size() >= n_coeffs);

    // needs generous sample/rate oversampling
    freq = std ::min (sr < 176200. ? 12000. : 20000., freq);
    freq = std ::min (sr < 96000. ? 11000. : 20000., freq);
    freq = std ::min (sr < 88200. ? 6000. : 20000., freq);
    freq = std ::min (sr < 48000. ? 5500. : 20000., freq);

    double norm_freq = (freq / 20000.);
    double fc        = .5 * sr * norm_freq;
    double fs        = sr;

    co[k]    = reso * 120.;
    double g = std::tan (3.141592653589793 / fs * fc) / std::sqrt (1. + co[k]);
    co[vt2]  = 0.052;
    co[vt2i] = 19.23076923076923;
    double p0s = 1.0 / (1.0 + g);
    co[q0s]    = 1.0 - g;
    co[r1s]    = -g;
    co[k0s]    = co[vt2] * g * p0s;
    double nmp = (1.0 - p0s);
    double gn  = std::pow (nmp, 2.);
    double kgn = co[k] * gn;
    double p0g = 1.0 / (1.0 + kgn);
    co[k0g]    = -co[vt2i] * p0g;
    co[rg1]    = -2.0 * kgn;
    co[rg2]    = -1.0 * kgn;
    double tmp = p0s * (g - 1.0);
    double acc = tmp;
    co[qg1]    = -2.0 * (kgn + acc);
    acc *= tmp;
    co[qg2]               = -1.0 * (kgn + acc);
    *((u64*) &co[choice]) = filter_type;
    co[frac]              = fraction;
  }
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
}} // namespace artv::saike
