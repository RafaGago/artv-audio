#pragma once

#include <cassert>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// a slew limiter. Notice that an envelope follower can be made by modifying
// the input, and sending:
//
// tick (fabs (in)): peak follower
// tick (in * in): mean abs
// sqrt (tick (in * in)): instant RMS
//
struct slew_limiter {
  //----------------------------------------------------------------------------
  enum coeffs { attack, release, n_coeffs };
  enum state { prev, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void init (crange<T> c, T attack_sec, T release_sec, T samplerate)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (attack_sec >= 0.);
    assert (release_sec >= 0.);

    double_x2 v = {attack_sec, release_sec};
    v *= samplerate;
    v          = vec_exp ((T) -1. / v);
    c[attack]  = v[0];
    c[release] = v[1];
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (
    crange<vec_value_type_t<V>> c,
    V                           attack_sec,
    V                           release_sec,
    vec_value_type_t<V>         samplerate)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    attack_sec *= samplerate;
    attack_sec = vec_exp ((T) -1. / attack_sec);
    vec_store (&c[attack * traits.size], attack_sec);

    release_sec *= samplerate;
    release_sec = vec_exp ((T) -1. / release_sec);
    vec_store (&c[release * traits.size], release_sec);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    T in_prev = s[prev];
    T coeff   = in_prev < in ? c[attack] : c[release];
    s[prev]   = in + coeff * (in_prev - in);
    return s[prev];
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>>
                                c, // coeffs interleaved, ready to SIMD load
    crange<vec_value_type_t<V>> s, // states interleaved, ready to SIMD load
    V                           in)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (traits.size * n_coeffs));
    assert (s.size() >= (traits.size * n_states));

    V in_prev = vec_load<V> (&s[prev * traits.size]);
    V att     = vec_load<V> (&c[attack * traits.size]);
    V rel     = vec_load<V> (&c[release * traits.size]);

    auto coeffs = in_prev < in ? att : rel;
    in_prev     = in + coeffs * (in_prev - in);
    vec_store (&s[prev * traits.size], in_prev);
    return in_prev;
  }
};
//------------------------------------------------------------------------------

} // namespace artv
