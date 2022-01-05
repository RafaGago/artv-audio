#pragma once

#include <cmath>

#include "artv-common/dsp/own/parts/traits.hpp"

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// A TDF2 biquad. For audio use the trapezoidal SVF instead. The goal of this
// class is to be able to extract poles and zeros from the coefficients.
//
// Coefficient calculation from RBJ's cookbook.
// TODO: use M.Vicaneks matched biquads:
// https://www.vicanek.de/articles/BiquadFits.pdf
//------------------------------------------------------------------------------

struct biquad {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, b2, a1, a2, n_coeffs };
  enum states { z1, z2, n_states };

  using rev_zeros_tag = part_tick_tag<0>;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    lowpass_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    highpass_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    bandpass_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    notch_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    allpass_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    bell_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    lowshelf_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
    vec_value_type_t<V> sr,
    highshelf_tag)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    V A     = vec_exp (gain_db * T {1. / 40.} * T {M_LN10});
    V w0    = (T) 2. * (T) M_PI * freq / sr;
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
  static void get_poles (crange<const V> co, crange<vec_complex<V>> poles)
  {
    assert (poles.size() >= 2);

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
};
//------------------------------------------------------------------------------
} // namespace artv