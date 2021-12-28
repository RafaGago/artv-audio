#pragma once

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// https://www.musicdsp.org/en/latest/Synthesis/10-fast-sine-and-cosine-calculation.html
//------------------------------------------------------------------------------
class cheap_sin_osc {
public:
  //----------------------------------------------------------------------------
  enum coeffs { a, n_coeffs };
  enum state { z0, z1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           c,
    V                   freq,
    vec_value_type_t<V> sr,
    highpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (c.size() >= n_coeffs);

    c[a] = (T) 2. * vec_sin ((T) M_PI * freq / srate);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    using T = vec_value_type_t<V>;
    reset_states (st, vec_set<V> ((T) 1.0));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, V bipolar_ampl)
  {
    assert (st.size() >= n_states);

    st[z0] = bipolar_ampl;
    st[z1] = vec_set<V> ((T) 0.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st) // states (interleaved, SIMD aligned)
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
