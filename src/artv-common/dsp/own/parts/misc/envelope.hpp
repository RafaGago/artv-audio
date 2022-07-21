#pragma once

#pragma once

#include <cassert>
#include <limits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
// an envelope follower/smoother.
//
struct envelope {
  //----------------------------------------------------------------------------
  enum coeffs { time_k, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { prev, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   time_sec,
    vec_value_type_t<V> samplerate)
  {
    using T = vec_value_type_t<V>;
    assert (c.size() >= n_coeffs);

    V    k    = vec_exp ((T) -1. / (time_sec * samplerate));
    auto zero = vec_set<V> ((T) 0);
    c[time_k] = k != zero ? k : zero;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> s, V in)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    s[prev] = in + c[time_k] * (s[prev] - in);
    return s[prev];
  }
  //----------------------------------------------------------------------------
  struct rms_tag {};
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> s, V in, rms_tag)
  {
    using T = vec_value_type_t<V>;
    auto v  = tick (c, s, in * in);
    return vec_sqrt (vec_max ((T) 0, v));
  }
};
//------------------------------------------------------------------------------

} // namespace artv
