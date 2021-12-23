#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <cmath>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/dsp/own/parts/traits.hpp"

#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

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

template <class T, class U>
static svf_coeffs_ext<T> get_main_coeffs (T freq, T q, T db, U sr, bool is_bell)
{
  // "U" should always be a builtin type.
  static_assert (std::is_floating_point<U>::value, "");

  svf_coeffs_ext<T> ret;
  T                 g;

  if constexpr (is_vec_v<T>) {
    ret.A = vec_exp (db * (U) (1. / 40.) * (U) M_LN10);
    g     = vec_tan (M_PI * freq / sr);
  }
  else {
    ret.A = exp (db * (U) (1. / 40.) * (U) M_LN10);
    g     = tan (M_PI * freq / sr);
  }
  if (is_bell) {
    q *= ret.A;
  }
  ret.k  = (U) 1.0 / q;
  ret.a1 = (U) 1.0 / ((U) 1.0 + g * (g + ret.k));
  ret.a2 = g * ret.a1;
  ret.a3 = g * ret.a2;
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
    || mp11::mp_contains<enabled_modes, allpass_tag>::value;

  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;

  static_assert (
    mp11::mp_size<enabled_modes>::value >= 1,
    "One mode needs to be enabled");

  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, k, n_coeffs = k + (uint) needs_k_coeff };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {

    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    if constexpr (needs_k_coeff) {
      vec_store (&c[k * traits.size], coeffs.k);
    }
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);

    if constexpr (needs_k_coeff) {
      return tick (
        s,
        v0,
        vec_load<V> (&c[a1 * traits.size]),
        vec_load<V> (&c[a2 * traits.size]),
        vec_load<V> (&c[a3 * traits.size]),
        vec_load<V> (&c[k * traits.size]));
    }
    else {
      return tick (
        s,
        v0,
        vec_load<V> (&c[a1 * traits.size]),
        vec_load<V> (&c[a2 * traits.size]),
        vec_load<V> (&c[a3 * traits.size]));
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0,
    single_coeff_set_tag)
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
    crange<vec_value_type_t<V>> s,
    V                           v0,
    V                           a1_v,
    V                           a2_v,
    V                           a3_v,
    [[maybe_unused]] V          k_v = vec_set<V> (0.))
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (s.size() >= traits.size * n_states);

    auto tick_r = detail::svf_tick (
      vec_load<V> (&s[ic1eq * traits.size]),
      vec_load<V> (&s[ic2eq * traits.size]),
      a1_v,
      a2_v,
      a3_v,
      v0);
    vec_store (&s[ic1eq * traits.size], tick_r.ic1eq);
    vec_store (&s[ic2eq * traits.size], tick_r.ic2eq);

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, lowpass_tag>) {
        ret[index.value] = tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, highpass_tag>) {
        ret[index.value] = v0 - k_v * tick_r.v1 - tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, bandpass_tag>) {
        ret[index.value] = tick_r.v1;
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
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 0.));
    vec_store (&c[m1 * traits.size], vec_set<V> ((T) 0.));
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], -coeffs.k);
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) -1.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr,
    peak_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= (n_coeffs * traits.size));

    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], -coeffs.k);
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) -2.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr,
    allpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, sr);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], (T) -2. * coeffs.k);
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 0.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    V                           db,
    vec_value_type_t<V>         sr,
    bell_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, db, sr, true);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], coeffs.k * (coeffs.A * coeffs.A - (T) 1.));
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 0.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    V                           db,
    vec_value_type_t<V>         sr,
    lowshelf_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, db, sr, false);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], coeffs.k * (coeffs.A - 1.));
    vec_store (&c[m2 * traits.size], (coeffs.A * coeffs.A) - 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    V                           db,
    vec_value_type_t<V>         sr,
    highshelf_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    auto           coeffs = detail::get_main_coeffs (freq, q, db, sr, false);

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a1 * traits.size], coeffs.a1);
    vec_store (&c[a2 * traits.size], coeffs.a2);
    vec_store (&c[a3 * traits.size], coeffs.a3);
    vec_store (&c[m0 * traits.size], coeffs.A * coeffs.A);
    vec_store (&c[m1 * traits.size], coeffs.k * (1. - coeffs.A) * coeffs.A);
    vec_store (&c[m2 * traits.size], 1. - (coeffs.A * coeffs.A));
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
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0,
    single_coeff_set_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    assert (s.size() >= traits.size * n_states);

    V ic1eq_v = vec_load<V> (&s[ic1eq * traits.size]);
    V ic2eq_v = vec_load<V> (&s[ic2eq * traits.size]);

    V out = calc<V> (
      ic1eq_v,
      ic2eq_v,
      vec_set<V> (c[m0]),
      vec_set<V> (c[m1]),
      vec_set<V> (c[m2]),
      vec_set<V> (c[a1]),
      vec_set<V> (c[a2]),
      vec_set<V> (c[a3]),
      v0);

    vec_store (&s[ic1eq * traits.size], ic1eq_v);
    vec_store (&s[ic2eq * traits.size], ic2eq_v);

    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);
    assert (s.size() >= traits.size * n_states);

    V ic1eq_v = vec_load<V> (&s[ic1eq * traits.size]);
    V ic2eq_v = vec_load<V> (&s[ic2eq * traits.size]);

    V out = calc<V> (
      ic1eq_v,
      ic2eq_v,
      vec_load<V> (&c[m0 * traits.size]),
      vec_load<V> (&c[m1 * traits.size]),
      vec_load<V> (&c[m2 * traits.size]),
      vec_load<V> (&c[a1 * traits.size]),
      vec_load<V> (&c[a2 * traits.size]),
      vec_load<V> (&c[a3 * traits.size]),
      v0);

    vec_store (&s[ic1eq * traits.size], ic1eq_v);
    vec_store (&s[ic2eq * traits.size], ic2eq_v);

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
