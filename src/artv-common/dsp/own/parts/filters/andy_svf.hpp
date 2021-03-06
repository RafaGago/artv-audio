#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <cmath>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/dsp/own/parts/traits.hpp"

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv { namespace andy {
namespace detail {
//------------------------------------------------------------------------------
template <class T>
struct svf_coeffs {
  T a1;
  T a2;
  T a3;
  T k;
};

template <class T, class U>
static svf_coeffs<T> get_main_coeffs (T freq, T q, U sr)
{
  // "U" should always be a builtin type.
  static_assert (std::is_floating_point<U>::value, "");

  svf_coeffs<T> ret;

  T g    = vec_tan ((U) M_PI * freq / sr);
  ret.k  = (U) 1.0 / q;
  ret.a1 = (U) 1.0 / ((U) 1.0 + g * (g + ret.k));
  ret.a2 = g * ret.a1;
  ret.a3 = g * ret.a2;
  return ret;
}
//------------------------------------------------------------------------------
template <class T>
struct svf_coeffs_ext : public svf_coeffs<T> {
  T A;
};

static constexpr uint bell_flag   = 1 << 0;
static constexpr uint lshelf_flag = 1 << 1;
static constexpr uint hshelf_flag = 1 << 2;

template <class T, class U>
static svf_coeffs_ext<T> get_main_coeffs (
  T    freq,
  T    q,
  T    db,
  U    sr,
  uint flags = 0)
{
  // "U" should always be a builtin type.
  static_assert (std::is_floating_point<U>::value, "");

  svf_coeffs_ext<T> ret;
  T                 g;

  if constexpr (is_vec_v<T>) {
    ret.A = vec_exp (db * (U) (1. / 40.) * (U) M_LN10);
    g     = vec_tan ((U) M_PI * freq / sr);
  }
  else {
    ret.A = exp (db * (U) (1. / 40.) * (U) M_LN10);
    g     = tan ((U) M_PI * freq / sr);
  }
  if ((flags & bell_flag)) {
    q *= ret.A;
  }
  if ((flags & lshelf_flag)) {
    g /= vec_sqrt (ret.A);
  }
  if ((flags & hshelf_flag)) {
    g *= vec_sqrt (ret.A);
  }
  ret.k  = (U) 1.0 / q;
  ret.a1 = (U) 1.0 / ((U) 1.0 + g * (g + ret.k));
  ret.a2 = g * ret.a1;
  ret.a3 = g * ret.a2;
  return ret;
}

// https://cytomic.com/files/dsp/SvfLinearTrapezoidalSin.pdf
template <class T, class U>
static svf_coeffs_ext<T> get_main_coeffs_precise (
  T    freq,
  T    q,
  T    db,
  U    sr,
  uint flags = 0)
{
  // "U" should always be a builtin type.
  static_assert (std::is_floating_point<U>::value, "");

  svf_coeffs_ext<T> ret;
  T                 s1, s2;

  T w = (U) M_PI * freq / sr;
  if constexpr (is_vec_v<T>) {
    ret.A = vec_exp (db * (U) (1. / 40.) * (U) M_LN10);
    s1    = vec_sin (w);
    // sometimes sin and cos can be returned on the same op.
    // sin(2x) = 2*cos(x)*sin(x)
    s2 = (U) 2 * vec_cos (w) * s1;
  }
  else {
    ret.A = exp (db * (U) (1. / 40.) * (U) M_LN10);
    s1    = vec_sin (w);
    // sometimes sin and cos can be returned on the same op.
    // sin(2x) = 2*cos(x)*sin(x)
    s2 = (U) 2 * vec_cos (w) * s1;
  }
  T nrm = (U) 1 / ((U) 2 * s2);

  ret.k  = (U) 1 / q;
  ret.a1 = s2 * nrm;
  ret.a2 = ((U) -2 * s1 * s1 - ret.k * s2) * nrm;
  ret.a3 = ((U) 2 * s1 * s1) * nrm;
  return ret;
}
//------------------------------------------------------------------------------
template <class T>
struct svf_tick_result {
  T v1;
  T v2;
  T ic1eq;
  T ic2eq;
};

template <class T>
static svf_tick_result<T> svf_tick (T ic1eq, T ic2eq, T a1, T a2, T a3, T v0)
{
  svf_tick_result<T> ret;

  T v3      = v0 - ic2eq;
  ret.v1    = a1 * ic1eq + a2 * v3;
  ret.v2    = ic2eq + a2 * ic1eq + a3 * v3;
  ret.ic1eq = (ret.v1 * 2.) - ic1eq;
  ret.ic2eq = (ret.v2 * 2.) - ic2eq;
  return ret;
}

//------------------------------------------------------------------------------
} // namespace detail

//------------------------------------------------------------------------------
// Enables the SVF multimode output. It is compile-time optimized depending on
// the enabled modes. Use any of the "_tag" classes as template parameters to
// enable a given mode. The outputs are returned on the "tick" function on the
// order on which the tags appear (natural order, left to right). If only one
// tag is passed the output is returned as a single value instead of an array.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template <class... Tags>
struct svf_multimode {

  // maybe do tag validation, (things won't compile anyways).
  using enabled_modes = mp11::mp_unique<mp11::mp_list<Tags...>>;

  static constexpr bool needs_k_coeff
    = mp11::mp_contains<enabled_modes, highpass_tag>::value
    || mp11::mp_contains<enabled_modes, notch_tag>::value
    || mp11::mp_contains<enabled_modes, peak_tag>::value
    || mp11::mp_contains<enabled_modes, bandpass_tag>::value
    || mp11::mp_contains<enabled_modes, allpass_tag>::value;

  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;

  static_assert (
    mp11::mp_size<enabled_modes>::value >= 1,
    "One mode needs to be enabled");

  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, k, n_coeffs = k + (uint) needs_k_coeff };
  enum coeffs_int { n_coeffs_int };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, V q, vec_value_type_t<V> sr)
  {
    assert (c.size() >= n_coeffs);

    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    if constexpr (needs_k_coeff) {
      c[k] = coeffs.k;
    }
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
  static auto tick (
    crange<const V> c, // coeffs (interleaved, SIMD aligned)
    crange<V>       s, // states (interleaved, SIMD aligned)
    V               v0)
  {
    assert (c.size() >= n_coeffs);

    if constexpr (needs_k_coeff) {
      return tick (s, v0, c[a1], c[a2], c[a3], c[k]);
    }
    else {
      return tick (s, v0, c[a1], c[a2], c[a3]);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> c, // coeffs (single set)
    crange<V>                         s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    assert (c.size() >= n_coeffs);

    if constexpr (needs_k_coeff) {
      return tick (
        s,
        v0,
        vec_set<V> (c[a1]),
        vec_set<V> (c[a2]),
        vec_set<V> (c[a3]),
        vec_set<V> (c[k]));
    }
    else {
      return tick (
        s, v0, vec_set<V> (c[a1]), vec_set<V> (c[a2]), vec_set<V> (c[a3]));
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<V>          s,
    V                  v0,
    V                  a1_v,
    V                  a2_v,
    V                  a3_v,
    [[maybe_unused]] V k_v = vec_set<V> (0.))
  {
    using T = vec_value_type_t<V>;
    assert (s.size() >= n_states);

    auto tick_r = detail::svf_tick (s[ic1eq], s[ic2eq], a1_v, a2_v, a3_v, v0);
    s[ic1eq]    = tick_r.ic1eq;
    s[ic2eq]    = tick_r.ic2eq;

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, lowpass_tag>) {
        ret[index.value] = tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, highpass_tag>) {
        ret[index.value] = v0 - k_v * tick_r.v1 - tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, bandpass_q_gain_tag>) {
        ret[index.value] = tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, bandpass_tag>) {
        ret[index.value] = tick_r.v1 * k_v;
      }
      else if constexpr (std::is_same_v<tag, notch_tag>) {
        ret[index.value] = v0 - k_v * tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, peak_tag>) {
        ret[index.value] = v0 - k_v * tick_r.v1 - (T) 2. * tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, allpass_tag>) {
        ret[index.value] = v0 - (T) 2. * k_v * tick_r.v1;
      }
    });

    if constexpr (returns_array) {
      return ret;
    }
    else {
      return ret[0];
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
using svf_lowpass  = svf_multimode<lowpass_tag>;
using svf_highpass = svf_multimode<highpass_tag>;
using svf_allpass  = svf_multimode<allpass_tag>;
using svf_bandpass = svf_multimode<bandpass_tag>;
using svf_notch    = svf_multimode<notch_tag>;
using svf_peak     = svf_multimode<peak_tag>;
//------------------------------------------------------------------------------
struct svf {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, m0, m1, m2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    lowpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = vec_set<V> ((T) 0.);
    c[m2] = vec_set<V> ((T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    highpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = -coeffs.k;
    c[m2] = vec_set<V> ((T) -1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    bandpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = coeffs.k;
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    bandpass_q_gain_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = vec_set<V> ((T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    peak_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);

    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = -coeffs.k;
    c[m2] = vec_set<V> ((T) -2.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    vec_value_type_t<V> sr,
    allpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = (T) -2. * coeffs.k;
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> sr,
    bell_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, db, sr, detail::bell_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = coeffs.k * (coeffs.A * coeffs.A - (T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> sr,
    bell_bandpass_tag)
  {
    using T     = vec_value_type_t<V>;
    auto coeffs = detail::get_main_coeffs (freq, q, db, sr, detail::bell_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 0.);
    c[m1] = coeffs.k * (coeffs.A * coeffs.A - (T) 1.);
    c[m2] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> sr,
    lowshelf_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, sr, detail::lshelf_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = vec_set<V> ((T) 1.);
    c[m1] = coeffs.k * (coeffs.A - (T) 1.);
    c[m2] = (coeffs.A * coeffs.A) - (T) 1.;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    V                   q,
    V                   db,
    vec_value_type_t<V> sr,
    highshelf_tag)
  {
    using T = vec_value_type_t<V>;
    auto coeffs
      = detail::get_main_coeffs (freq, q, db, sr, detail::hshelf_flag);

    assert (c.size() >= n_coeffs);
    c[a1] = coeffs.a1;
    c[a2] = coeffs.a2;
    c[a3] = coeffs.a3;
    c[m0] = coeffs.A * coeffs.A;
    c[m1] = coeffs.k * ((T) 1. - coeffs.A) * coeffs.A;
    c[m2] = (T) 1. - (coeffs.A * coeffs.A);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<V>                         s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V out = calc<V> (
      s[ic1eq],
      s[ic2eq],
      vec_set<V> (c[m0]),
      vec_set<V> (c[m1]),
      vec_set<V> (c[m2]),
      vec_set<V> (c[a1]),
      vec_set<V> (c[a2]),
      vec_set<V> (c[a3]),
      v0);

    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> c, // coeffs (interleaved, SIMD aligned)
    crange<V>       s, // states (interleaved, SIMD aligned)
    V               v0)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    V out = calc<V> (
      s[ic1eq], s[ic2eq], c[m0], c[m1], c[m2], c[a1], c[a2], c[a3], v0);

    return out;
  }
  //----------------------------------------------------------------------------
private:
  template <class T>
  static T calc (
    T& ic1eq_v,
    T& ic2eq_v,
    T  m0_v,
    T  m1_v,
    T  m2_v,
    T  a1_v,
    T  a2_v,
    T  a3_v,
    T  v0)
  {
    auto tick_r = detail::svf_tick (ic1eq_v, ic2eq_v, a1_v, a2_v, a3_v, v0);
    T    out    = m0_v * v0 + m1_v * tick_r.v1 + m2_v * tick_r.v2;
    ic1eq_v     = tick_r.ic1eq;
    ic2eq_v     = tick_r.ic2eq;
    return out;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

}} // namespace artv::andy
