#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/filters/onepole.hpp"
#include "artv-common/misc/simd.hpp"

//------------------------------------------------------------------------------
namespace artv {

// Classic IIR DC blocker based on a highpass
// https://ccrma.stanford.edu/~jos/filters/DC_Blocker.html
// https://www.musicdsp.org/en/latest/Filters/135-dc-filter.html?highlight=DC
struct iir_dc_blocker {
  //----------------------------------------------------------------------------
  enum coeffs { R, n_coeffs };
  enum state { x1, y1, n_states };
  //----------------------------------------------------------------------------
  // warning, if going to very low frequencies, use "double".
  template <class T>
  static void init (crange<T> c, T freq, T sr)
  {
    assert (c.size() >= n_coeffs);
    c[R] = 1. - (M_PI * 2. * freq / sr);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));

    auto Rv = 1. - (M_PI * 2. * freq / sr);
    vec_store (&c[R * traits.size], Rv);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               x)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    T y   = x - s[x1] + c[R] * s[y1];
    s[y1] = y;
    s[x1] = x;
    return y;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       s,
    V                                 x)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));
    assert (s.size() >= (traits.size * n_states));

    auto x1_ptr = &s[x1 * traits.size];
    auto y1_ptr = &s[y1 * traits.size];

    auto x1v = vec_load<V> (x1_ptr);
    auto y1v = vec_load<V> (y1_ptr);
    auto Rv  = vec_load<V> (&c[R * traits.size]);

    auto y = x - x1v + Rv * y1v;

    vec_store (x1_ptr, x);
    vec_store (y1_ptr, y);
    return y;
  }
  //----------------------------------------------------------------------------
};

// I don't know if he is the original author but he does an extreme good job
// sharing knowledge and helping. Described here:
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=545280#top
struct mystran_dc_blocker {
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs = onepole_smoother::n_coeffs;
  static constexpr uint n_states = onepole_smoother::n_states;
  //----------------------------------------------------------------------------
  // warning, if going to very low frequencies, use "double".
  template <class T>
  static void init (crange<T> c, T freq, T sr)
  {
    onepole_smoother::lowpass (c, freq, sr);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    onepole_smoother::lowpass_simd (c, freq, sr);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               x)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto prev_lp_out = s[onepole_smoother::z1];
    onepole_smoother::tick (c, s, x);
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c,
    crange<vec_value_type_t<V>>       s,
    V                                 x)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));
    assert (s.size() >= (traits.size * n_states));

    auto prev_lp_out = vec_load<V> (&s[onepole_smoother::z1 * traits.size]);
    onepole_smoother::tick_simd (c, s, x);
    return x - prev_lp_out;
  }
  //----------------------------------------------------------------------------
};
// TODO: Things to try
// -Moving average DC blocker (linear phase):
//    https://www.dsprelated.com/showthread/comp.dsp/66509-1.php
//------------------------------------------------------------------------------
} // namespace artv
