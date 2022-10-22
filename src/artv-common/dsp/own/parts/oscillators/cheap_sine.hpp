#pragma once

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// https://www.musicdsp.org/en/latest/Synthesis/10-fast-sine-and-cosine-calculation.html
//------------------------------------------------------------------------------
class cheap_sin_osc {
public:
  //----------------------------------------------------------------------------
  enum coeffs { a, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { z0, z1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            c,
    V                   freq,
    vec_value_type_t<V> t_spl,
    highpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (c.size() >= n_coeffs);

    c[a] = (T) 2. * vec_sin ((T) M_PI * freq * t_spl);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    using T = vec_value_type_t<V>;
    reset_states (st, vec_set<V> ((T) 1.0));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st, V bipolar_ampl)
  {
    assert (st.size() >= n_states);

    st[z0] = bipolar_ampl;
    st[z1] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    xspan<V const> co, // coeffs (interleaved, SIMD aligned)
    xspan<V>       st) // states (interleaved, SIMD aligned)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    st[z0] = st[z0] - co[a] * st[z1];
    st[z1] = st[z1] + co[a] * st[z0];
    return st[z9];
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
