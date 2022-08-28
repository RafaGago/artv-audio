#pragma once

#include <cmath>

#include "artv-common/dsp/own/parts/traits.hpp"

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

namespace artv {

//------------------------------------------------------------------------------
// A TDF2 biquad. For audio use the trapezoidal SVF instead. The goal of this
// class is mostly to be able to extract poles and zeros from the coefficients.
//------------------------------------------------------------------------------
struct biquad {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, b2, a1, a2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum states { z1, z2, n_states };

  using rev_zeros_tag = part_tick_tag<0>;
  //----------------------------------------------------------------------------
  // RBJ
  //----------------------------------------------------------------------------
  template <class T>
  struct rbj_tag {};

  using rbj_lowpass_tag   = rbj_tag<lowpass_tag>;
  using rbj_highpass_tag  = rbj_tag<highpass_tag>;
  using rbj_bandpass_tag  = rbj_tag<bandpass_tag>;
  using rbj_notch_tag     = rbj_tag<notch_tag>;
  using rbj_peak_tag      = rbj_tag<peak_tag>;
  using rbj_allpass_tag   = rbj_tag<allpass_tag>;
  using rbj_bell_tag      = rbj_tag<bell_tag>;
  using rbj_lowshelf_tag  = rbj_tag<lowshelf_tag>;
  using rbj_highshelf_tag = rbj_tag<highshelf_tag>;

  template <class T>
  struct mvic_hq_tag {};

  using mvic_lowpass_hq_tag  = mvic_hq_tag<lowpass_tag>;
  using mvic_highpass_hq_tag = mvic_hq_tag<highpass_tag>;
  using mvic_bandpass_hq_tag = mvic_hq_tag<bandpass_tag>;
  using mvic_bell_hq_tag     = mvic_hq_tag<bell_tag>;

  template <class T>
  struct mvic_tag {};

  using mvic_lowpass_tag  = mvic_tag<lowpass_tag>;
  using mvic_highpass_tag = mvic_tag<highpass_tag>;
  using mvic_bandpass_tag = mvic_tag<bandpass_tag>;

  //----------------------------------------------------------------------------
  // Defaults
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    lowpass_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_lowpass_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    highpass_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_highpass_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    bandpass_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_bandpass_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    notch_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_notch_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    peak_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_peak_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    allpass_tag)
  {
    reset_coeffs<V> (co, freq, q, t_spl, rbj_allpass_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    bell_tag)
  {
    reset_coeffs<V> (co, freq, q, gain_db, t_spl, rbj_bell_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    lowshelf_tag)
  {
    reset_coeffs<V> (co, freq, q, gain_db, t_spl, rbj_lowshelf_tag {});
  }
  //----------------------------------------------------------------------------
  template <class V>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    highshelf_tag)
  {
    reset_coeffs<V> (co, freq, q, gain_db, t_spl, rbj_highshelf_tag {});
  }
  //----------------------------------------------------------------------------
  // thiran interpolator of order 2, wastes memory and ops on unrequired
  // coefficients...
  template <class V>
  static void reset_coeffs (
    crange<V> co,
    V         d, // 0 to 1
    thiran_tag)
  {
    using T = vec_value_type_t<V>;

    // http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt3_allpass.pdf
    d += (T) 1.403;
    co[a1] = (T) -2 * (d - (T) 2) / (d + (T) 1); // a1
    co[a2] = ((d - (T) 1) * (d - (T) 2)) / ((d + (T) 1) * (d + (T) 2)); // a2
    co[b0] = co[a2];
    co[b1] = co[a1];
    co[b2] = vec_set<V> ((T) 0);
  }
  //----------------------------------------------------------------------------
  // Robert Bristow Jhonson's
  // https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    rbj_lowpass_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V a0   = (T) 1 + alpha;
    V norm = (T) 1 / a0;

    co[b0] = co[b2] = (((T) 1 - cosw0) * (T) 0.5) * norm;
    co[b1]          = ((T) 1 - cosw0) * norm;
    co[a1]          = ((T) -2 * cosw0) * norm;
    co[a2]          = ((T) 1 - alpha) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    rbj_highpass_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V a0   = (T) 1 + alpha;
    V norm = (T) 1 / a0;

    co[b0] = co[b2] = (((T) 1 + cosw0) * (T) 0.5) * norm;
    co[b1]          = (-((T) 1 + cosw0)) * norm;
    co[a1]          = ((T) -2 * cosw0) * norm;
    co[a2]          = ((T) 1 - alpha) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    rbj_bandpass_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V a0   = (T) 1 + alpha;
    V norm = (T) 1 / a0;

    co[b0] = alpha * norm;
    co[b1] = vec_set<V> ((T) 0);
    co[b2] = -co[b0];
    co[a1] = ((T) -2 * cosw0) * norm;
    co[a2] = ((T) 1 - alpha) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    rbj_notch_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V a0   = (T) 1 + alpha;
    V norm = (T) 1 / a0;

    co[b0] = co[b2] = norm;
    co[b1] = co[a1] = ((T) -2 * cosw0) * norm;
    co[a2]          = ((T) 1 - alpha) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    rbj_allpass_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V a0   = (T) 1 + alpha;
    V norm = (T) 1 / a0;

    co[b0] = co[a2] = ((T) 1 - alpha) * norm;
    co[b1] = co[a1] = ((T) -2 * cosw0) * norm;
    co[b2]          = vec_set<V> ((T) 1); // 1 + alpha / a0
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    rbj_bell_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V alpha_div_A = alpha / A;
    V a0          = (T) 1 + alpha_div_A;
    V norm        = (T) 1 / a0;

    co[b0] = ((T) 1 + alpha * A) * norm;
    co[b1] = co[a1] = ((T) -2 * cosw0) * norm;
    co[b2]          = ((T) 1 - alpha * A) * norm;
    co[a2]          = ((T) 1 - alpha_div_A) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    rbj_lowshelf_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V g    = (T) 2 * vec_sqrt (A) * alpha;
    V a0   = (A + (T) 1) + (A - (T) 1) * cosw0 + g;
    V norm = (T) 1 / a0;

    co[b0] = A * ((A + (T) 1) - (A - (T) 1) * cosw0 + g) * norm;
    co[b1] = (T) 2. * A * ((A - (T) 1) - (A + (T) 1) * cosw0) * norm;
    co[b2] = A * ((A + (T) 1) - (A - (T) 1) * cosw0 - g) * norm;
    co[a1] = (T) -2. * ((A - (T) 1) + (A + (T) 1) * cosw0) * norm;
    co[a2] = ((A + (T) 1) + (A - (T) 1) * cosw0 - g) * norm;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    rbj_highshelf_tag)
  {
    // RBJ's
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq * t_spl;
    V cosw0 = vec_cos (w0);
    V alpha = vec_sin (w0) / ((T) 2. * q);

    V g    = (T) 2 * vec_sqrt (A) * alpha;
    V a0   = (A + (T) 1) - (A - (T) 1) * cosw0 + g;
    V norm = (T) 1 / a0;

    co[b0] = A * ((A + (T) 1) + (A - (T) 1) * cosw0 + g) * norm;
    co[b1] = (T) -2. * A * ((A - (T) 1) + (A + (T) 1) * cosw0) * norm;
    co[b2] = A * ((A + (T) 1) + (A - (T) 1) * cosw0 - g) * norm;
    co[a1] = (T) 2. * ((A - (T) 1) - (A + (T) 1) * cosw0) * norm;
    co[a2] = ((A + (T) 1) - (A - (T) 1) * cosw0 - g) * norm;
  }
  //----------------------------------------------------------------------------
  // Martin Vicanek's loose fit versions
  // https://www.mvic.de/articles/BiquadFits.pdf
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_lowpass_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);

    V r0    = (T) 1.0 + co[a1] + co[a2];
    V f02   = f0 * f0;
    V omf02 = (T) 1.0 - f02;
    V num   = ((T) 1.0 - co[a1] + co[a2]) * f02;
    V r1    = num / vec_sqrt (omf02 * omf02 + f02 / (q * q));

    co[b0] = (r0 + r1) * (T) 0.5;
    co[b1] = r0 - co[b0];
    co[b2] = vec_set<V> (0.0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_highpass_tag)
  {

    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);

    V f02   = f0 * f0;
    V omf02 = (T) 1.0 - f02;
    V num   = (T) 1.0 - co[a1] + co[a2];
    V r1    = num / vec_sqrt (omf02 * omf02 + f02 / (q * q));

    co[b0] = r1 * (T) 0.25;
    co[b1] = (T) -2 * co[b0];
    co[b2] = co[b0];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_bandpass_tag)
  {
    // At first sight, this one seems to be breaking at low Q.
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);

    V r0 = ((T) 1.0 + co[a1] + co[a2]) / ((T) M_PI * f0 * q);

    V inv_q = (T) 1. / q;
    V f02   = f0 * f0;
    V omf02 = (T) 1.0 - f02;
    V num   = ((T) 1.0 - co[a1] + co[a2]) * f0 * inv_q;
    V r1    = num / vec_sqrt (omf02 * omf02 + f02 * inv_q * inv_q);

    co[b0] = -r1 * (T) 0.5;
    co[b1] = (r0 - co[b1]) * (T) 0.5;
    co[b2] = -co[b0] - co[b1];
  }
  //----------------------------------------------------------------------------
  // Martin Vicanek's tight fit versions
  // https://www.mvic.de/articles/BiquadFits.pdf
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_lowpass_hq_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);
    auto m = get_mvic_matched (co[a1], co[a2], f0 * (T) M_PI);

    V R1 = (m.A0 * m.phi0 + m.A1 * m.phi1 + m.A2 * m.phi2) * (q * q);
    V B0 = m.A0;
    V B1 = (R1 - B0 * m.phi0) / m.phi1;
    // according to the paper the sqrt on B1 requires double precision
    static_assert (sizeof (T) >= sizeof (double), "Not enough precision.");

    co[b0] = (T) 0.5 * (m.A0_sqrt + vec_sqrt (B1));
    co[b1] = m.A0_sqrt - co[b0];
    co[b2] = vec_set<V> (0.0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_highpass_hq_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);
    auto m = get_mvic_matched (co[a1], co[a2], f0 * (T) M_PI);

    // according to the paper the sqrt on B1 requires double precision
    static_assert (sizeof (T) >= sizeof (double), "Not enough precision.");

    V num  = vec_sqrt (m.A0 * m.phi0 + m.A1 * m.phi1 + m.A2 * m.phi2) * q;
    co[b0] = num / ((T) 4 * m.phi1);
    co[b1] = (T) -2 * co[b0];
    co[b2] = co[b0];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> t_spl,
    mvic_bandpass_hq_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);
    auto m = get_mvic_matched (co[a1], co[a2], f0 * (T) M_PI);

    V R1 = m.A0 * m.phi0 + m.A1 * m.phi1 + m.A2 * m.phi2;
    V R2 = -m.A0 + m.A1 + (T) 4 * (m.phi0 - m.phi1) * m.A2;

    V B2 = (R1 - R2 * m.phi1) / ((T) 4 * m.phi1 * m.phi1);
    V B1 = R2 + (T) 4 * (m.phi1 - m.phi0) * B2;

    // according to the paper some of the sqrts on B1 require double precision
    static_assert (sizeof (T) >= sizeof (double), "Not enough precision.");

    co[b1] = (T) -0.5 * vec_sqrt (B1);
    co[b0] = (T) 0.5 * (vec_sqrt (B2 + co[b1] * co[b1]) - co[b1]);
    co[b2] = -co[b0] - co[b1];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    V                   gain_db,
    vec_value_type_t<V> t_spl,
    mvic_bell_hq_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V f0 = (T) 2. * freq * t_spl; // normalized freq 0-1
    get_impulse_invariance_poles (co[a1], co[a2], f0, q);
    auto m = get_mvic_matched (co[a1], co[a2], f0 * (T) M_PI);

    V GG = vec_exp (gain_db * (T) (1. / 20.) * (T) M_LN10);
    GG *= GG;

    V R1 = (m.A0 * m.phi0 + m.A1 * m.phi1 + m.A2 * m.phi2) * GG;
    V R2 = (-m.A0 + m.A1 + (T) 4 * (m.phi0 - m.phi1) * m.A2) * GG;

    V B0 = m.A0;
    V B2 = (R1 - R2 * m.phi1 - B0) / ((T) 4 * m.phi1 * m.phi1);
    V B1 = R2 + B0 + (T) 4 * (m.phi1 - m.phi0) * B2;

    // according to the paper some of the sqrts on B1 require double precision
    static_assert (sizeof (T) >= sizeof (double), "Not enough precision.");

    V B0sqrt = vec_sqrt (B0);
    V W      = (T) 0.5 * (B0sqrt + vec_sqrt (B1));

    co[b0] = (T) 0.5 * (W + vec_sqrt (W * W + B2));
    co[b1] = (T) 0.5 * (B0sqrt - vec_sqrt (B1));
    co[b2] = -B2 / ((T) 4 * co[b0]);
  }
  //----------------------------------------------------------------------------
  // processing functions
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> co, crange<V> st, V x)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V out  = x * co[b0] + st[z1];
    st[z1] = x * co[b1] + st[z2] - co[a1] * out;
    st[z2] = x * co[b2] - co[a2] * out;
    return out;
  }
  //----------------------------------------------------------------------------
  // Running with revesed zeros.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> co, crange<V> st, V x, rev_zeros_tag)
  {
    assert (co.size() >= b2 + 1);
    assert (st.size() >= n_states);

    V out  = x * co[b2] + st[z1];
    st[z1] = x * co[b1] + st[z2];
    st[z2] = x * co[b0];
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void get_zeros (
    crange<const V>        co,
    crange<vec_complex<V>> zeros,
    V&                     zeros_gain)
  {
    assert (zeros.size() >= 2);

    get_quadratic_roots (zeros[0], zeros[1], co[b0], co[b1], co[b2]);
    zeros_gain = co[b0];
  } //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void get_poles (crange<vec_complex<V>> poles, crange<const V> co)
  {
    assert (poles.size() >= 2);

    // ffast-math is enabled, so the compiler can decide to remove the divisions
    // by 2 * a when seeing it is set to 1.
    get_quadratic_roots (poles[0], poles[1], vec_set<V> (1.), co[a1], co[a2]);
  }
  //----------------------------------------------------------------------------
private:
  // probably this belongs elsewhere...
  //----------------------------------------------------------------------------
  // quadratic equation roots:
  // (a * x^2) + (b * x) + c = a * (x - r1) * (x - r2)
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void get_quadratic_roots (
    vec_complex<V>& r1,
    vec_complex<V>& r2,
    V               a,
    V               b,
    V               c)
  {
    using T = vec_value_type_t<V>;

    // roots = -b / 2a ± sqrt (b^2 - 4ac) / 2a

    V inv_2a = (T) 1. / ((T) 2. * a);
    r1.re    = -b * inv_2a;
    r1.im    = (b * b) - (T) 4. * a * c;

    for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
      if (r1.im[i] < 0.) {
        // complex conjugate
        r1.im[i] = sqrt (-r1.im[i]) * inv_2a[i];
        r2.im[i] = -r1.im[i];
        r2.re[i] = r1.re[i];
      }
      else {
        // real poles
        r1.im[i] = sqrt (r1.im[i]) * inv_2a[i];
        r2.re[i] = r1.re[i] - r1.im[i];
        r1.re[i] = r1.re[i] + r1.im[i];
        r1.im[i] = r2.im[i] = (T) 0.;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void get_impulse_invariance_poles (V& a1, V& a2, V f0, V q)
  {
    using T = vec_value_type_t<V>;

    V i2q = (T) 1 / ((T) 2 * q);

    V w0        = f0 * (T) M_PI; // from 0-1 to 0-PI
    V exp_mq_w0 = vec_exp (-i2q * w0); // e^(-q w0)
    a2          = exp_mq_w0 * exp_mq_w0; // e^(-2 q w0)
    get_impulse_invariance_a1 (a1, w0, i2q, exp_mq_w0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void get_impulse_invariance_a1 (V& a1, V w0, V i2q, V exp_mq_w0)
  {
    using T = vec_value_type_t<V>;

    auto cmp = (i2q <= vec_set<V> ((T) 1));
    if (vec_is_all_ones (cmp)) {
      a1 = (T) -2 * exp_mq_w0 * vec_cos (vec_sqrt ((T) 1 - i2q * i2q) * w0);
      return;
    }
    if (vec_is_all_zeros (cmp)) {
      a1 = (T) -2 * exp_mq_w0 * vec_cosh (vec_sqrt (i2q * i2q - (T) 1) * w0);
      return;
    }
    for (uint i = 0; i < vec_traits<V>().size; ++i) {
      // do individually.
      auto a1x1 = vec_set<1> (a1[i]);
      get_impulse_invariance_a1 (
        a1x1,
        vec_set<1> (w0[i]),
        vec_set<1> (i2q[i]),
        vec_set<1> (exp_mq_w0[i]));
      a1[i] = a1x1[0];
    }
  }
  //----------------------------------------------------------------------------
  template <class V>
  struct mvic_matched {
    V A0_sqrt, A0, A1, A2, phi0, phi1, phi2;
  };
  //----------------------------------------------------------------------------
  template <class V>
  static mvic_matched<V> get_mvic_matched (V a1, V a2, V w0)
  {
    using T = vec_value_type_t<V>;

    V A0_sqrt = (T) 1.0 + a1 + a2;
    V A0      = A0_sqrt * A0_sqrt;
    V A1      = (T) 1.0 - a1 + a2;
    A1 *= A1;
    V A2 = (T) -4.0 * a2;

    V phi1 = vec_sin (w0 * (T) 0.5);
    phi1 *= phi1;
    V phi0 = (T) 1.0 - phi1;
    V phi2 = (T) 4.0 * phi0 * phi1;

    return {A0_sqrt, A0, A1, A2, phi0, phi1, phi2};
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
