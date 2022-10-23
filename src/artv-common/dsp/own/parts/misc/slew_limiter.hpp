#pragma once

#include <cassert>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// a slew limiter. Notice that an envelope follower can be made by modifying
// the input, and sending:
//
// tick (fabs (in)): peak follower
// tick (in * in): mean abs
// sqrt (tick (in * in)): instant RMS
struct slew_limiter {
  //----------------------------------------------------------------------------
  enum coeffs { attack, release, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { prev, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   attack_sec,
    V                   release_sec,
    vec_value_type_t<V> t_spl)
  {
    using T = vec_value_type_t<V>;
    assert (c.size() >= n_coeffs);

    auto zero = vec_set<V> ((T) 0);

    V atk     = vec_exp ((T) -t_spl / attack_sec);
    c[attack] = attack_sec != zero ? atk : zero;

    V rel      = vec_exp ((T) -t_spl / release_sec);
    c[release] = release_sec != zero ? rel : zero;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (xspan<V const> c, xspan<V> s, V in)
  {
    assert (c.size() >= n_coeffs);
    assert (s.size() >= n_states);

    auto coeff = s[prev] < in ? c[attack] : c[release];
    s[prev]    = in + coeff * (s[prev] - in);
    return s[prev];
  }
};
//------------------------------------------------------------------------------

} // namespace artv
