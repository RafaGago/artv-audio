#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <cmath>
#include <gcem.hpp>
#include <type_traits>

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
  T             g;

  if constexpr (is_vec_v<T>) {
    g = vec_tan (M_PI * freq / sr);
  }
  else {
    g = tan (M_PI * freq / sr);
  }
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
struct svf_multimode_tags {
  struct lowpass {};
  struct highpass {};
  struct bandpass {};
  struct notch {};
  struct peak {};
  struct allpass {};
};

template <class... Tags>
struct svf_multimode {

  using tags = svf_multimode_tags;

  // maybe do tag validation, (things won't compile anyways).
  using enabled_modes = mp11::mp_unique<mp11::mp_list<Tags...>>;

  static constexpr bool needs_k_coeff
    = mp11::mp_contains<enabled_modes, tags::highpass>::value
    || mp11::mp_contains<enabled_modes, tags::notch>::value
    || mp11::mp_contains<enabled_modes, tags::peak>::value
    || mp11::mp_contains<enabled_modes, tags::allpass>::value;

  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;

  static_assert (
    mp11::mp_size<enabled_modes>::value >= 1,
    "One mode needs to be enabled");

  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, k, n_coeffs = k + (uint) needs_k_coeff };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void init (crange<T> c, T freq, T q, T sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {

    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  static auto tick (crange<const T> c, crange<T> s, T v0)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto tick_r
      = detail::svf_tick (s[ic1eq], s[ic2eq], c[a1], c[a2], c[a3], v0);
    s[ic1eq] = tick_r.ic1eq;
    s[ic2eq] = tick_r.ic2eq;

    std::array<T, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, tags::lowpass>) {
        ret[index.value] = tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::highpass>) {
        ret[index.value] = v0 + -c[k] * tick_r.v1 + tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::bandpass>) {
        ret[index.value] = tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, tags::notch>) {
        ret[index.value] = v0 + -c[k] * tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, tags::peak>) {
        ret[index.value] = v0 + -c[k] * tick_r.v1 + (T) -2. * tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::allpass>) {
        ret[index.value] = v0 + (T) -2. * c[k] * tick_r.v1;
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
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static auto tick_simd (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);
    assert (s.size() >= traits.size * n_states);

    auto tick_r = detail::svf_tick (
      vec_load<V> (&s[ic1eq * traits.size]),
      vec_load<V> (&s[ic2eq * traits.size]),
      vec_load<V> (&c[a1 * traits.size]),
      vec_load<V> (&c[a2 * traits.size]),
      vec_load<V> (&c[a3 * traits.size]),
      v0);
    vec_store (&s[ic1eq * traits.size], tick_r.ic1eq);
    vec_store (&s[ic2eq * traits.size], tick_r.ic2eq);

    [[maybe_unused]] V k_v;
    if constexpr (needs_k_coeff) {
      vec_load (k_v, &c[k * traits.size]);
    }

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, tags::lowpass>) {
        ret[index.value] = tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::highpass>) {
        ret[index.value] = v0 + k_v * tick_r.v1 + tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::bandpass>) {
        ret[index.value] = tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, tags::notch>) {
        ret[index.value] = v0 + -k_v * tick_r.v1;
      }
      else if constexpr (std::is_same_v<tag, tags::peak>) {
        ret[index.value] = v0 + -k_v * tick_r.v1 + (T) -2. * tick_r.v2;
      }
      else if constexpr (std::is_same_v<tag, tags::allpass>) {
        ret[index.value] = v0 + (T) -2. * k_v * tick_r.v1;
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
using svf_lowpass  = svf_multimode<svf_multimode_tags::lowpass>;
using svf_highpass = svf_multimode<svf_multimode_tags::highpass>;
using svf_allpass  = svf_multimode<svf_multimode_tags::allpass>;
using svf_bandpass = svf_multimode<svf_multimode_tags::bandpass>;
using svf_notch    = svf_multimode<svf_multimode_tags::notch>;
using svf_peak     = svf_multimode<svf_multimode_tags::peak>;
//------------------------------------------------------------------------------
struct svf {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, m0, m1, m2, n_coeffs };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 0.;
    *co.m1 = 0.;
    *co.m2 = 1.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  static void bandpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 0.;
    *co.m1 = 1.;
    *co.m2 = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void highpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = -coeffs.k;
    *co.m2 = -1.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  static void notch (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = -coeffs.k;
    *co.m2 = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void peak (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = -coeffs.k;
    *co.m2 = -2.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void allpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, sr);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = -2. * coeffs.k;
    *co.m2 = 0.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void allpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  static void bell (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (interleaved_pack_offset);
    auto co     = unpack_interleaved_coeffs (c, interleaved_pack_size);
    auto coeffs = detail::get_main_coeffs (freq, q, db, sr, true);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = coeffs.k * (coeffs.A * coeffs.A - 1.);
    *co.m2 = 0.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void bell_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    V                           db,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  static void low_shelf (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (coeff_offset * coeff_distance * n_coeffs);
    auto co     = unpack_interleaved_coeffs (c, coeff_distance);
    auto coeffs = detail::get_main_coeffs (freq, q, db, sr, false);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = 1.;
    *co.m1 = coeffs.k * (coeffs.A - 1.);
    *co.m2 = (coeffs.A * coeffs.A) - 1.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void high_shelf (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c           = c.shrink_head (coeff_offset * coeff_distance * n_coeffs);
    auto co     = unpack_interleaved_coeffs (c, coeff_distance);
    auto coeffs = detail::get_main_coeffs (freq, q, db, sr, false);

    *co.a1 = coeffs.a1;
    *co.a2 = coeffs.a2;
    *co.a3 = coeffs.a3;
    *co.m0 = coeffs.A * coeffs.A;
    *co.m1 = coeffs.k * (1. - coeffs.A) * coeffs.A;
    *co.m2 = 1. - (coeffs.A * coeffs.A);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T> s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               v0)
  {
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    return calc (
      s[ic1eq], s[ic2eq], c[m0], c[m1], c[m2], c[a1], c[a2], c[a3], v0);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Sparse
  // states version.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>                              c, // coeffs
    std::array<crange<vec_value_type_t<V>>, vec_traits_t<V>::size> s, // states
    V                                                              v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    for (uint i = 0; i < traits.size; ++i) {
      assert (s[i].size() >= n_states);
      assert (s[i].size() >= n_states);
    }

    V ic1eq_v, ic2eq_v;

    for (uint i = 0; i < traits.size; ++i) {
      ic1eq_v[i] = s[i][ic1eq];
      ic2eq_v[i] = s[i][ic2eq];
    }

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

    for (uint i = 0; i < traits.size; ++i) {
      s[i][ic1eq] = ic1eq_v[i];
      s[i][ic2eq] = ic2eq_v[i];
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_aligned (
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  template <class T>
  struct coeff_ptrs {
    T* a1;
    T* a2;
    T* a3;
    T* m0;
    T* m1;
    T* m2;
  };
  //----------------------------------------------------------------------------
  template <class T>
  static coeff_ptrs<T> unpack_interleaved_coeffs (crange<T> c, uint distance)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= ((n_coeffs - 1) * distance) + 1);

    coeff_ptrs<T> ret;
    ret.a1 = &c[a1 * distance];
    ret.a2 = &c[a2 * distance];
    ret.a3 = &c[a3 * distance];
    ret.m0 = &c[m0 * distance];
    ret.m1 = &c[m1 * distance];
    ret.m2 = &c[m2 * distance];
    return ret;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

}} // namespace artv::andy
