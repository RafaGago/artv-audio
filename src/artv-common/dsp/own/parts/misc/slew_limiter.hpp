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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           attack_sec,
    V                           release_sec,
    vec_value_type_t<V>         samplerate)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    auto zero = vec_set<V> ((T) 0);

    V atk = attack_sec;
    atk *= samplerate;
    atk = vec_exp ((T) -1. / atk);
    atk = attack_sec != zero ? atk : zero;
    vec_store (&c[attack * traits.size], atk);

    V rel = release_sec;
    rel *= samplerate;
    rel = vec_exp ((T) -1. / rel);
    rel = release_sec != zero ? rel : zero;
    vec_store (&c[release * traits.size], rel);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
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
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>
                                c, // coeffs interleaved, ready to SIMD load
    crange<vec_value_type_t<V>> s, // states interleaved, ready to SIMD load
    V                           in)
  {
    using T               = vec_value_type_t<V>;
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
