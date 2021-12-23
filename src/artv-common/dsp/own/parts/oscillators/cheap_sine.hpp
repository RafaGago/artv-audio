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
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));
    vec_store (&c[a * traits.size], (T) 2. * vec_sin ((T) M_PI * freq / srate));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    reset_states (st, vec_set<V> ((vec_value_type_t<V>) 1.0));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, V bipolar_ampl)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= (n_states * traits.size));

    vec_store (&st[z0 * traits.size], bipolar_ampl);
    vec_store (&st[z1 * traits.size], vec_set<V> ((T) 0.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st) // states (interleaved, SIMD aligned)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= traits.size * n_coeffs);
    assert (st.size() >= traits.size * n_states);

    V a_v  = vec_load<V> (&co[a * traits.size]);
    V z0_v = vec_load<V> (&st[z0 * traits.size]);
    V z1_v = vec_load<V> (&st[z1 * traits.size]);

    z0_v = z0_v - a_v * z1_v;
    z1_v = z1_v + a_v * z0_v;

    vec_store (&st[z0 * traits.size], z0_v);
    vec_store (&st[z1 * traits.size], z1_v);
    return z0_v;
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
