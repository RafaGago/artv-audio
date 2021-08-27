#pragma once

#include "artv-common/misc/util.hpp"

namespace artv {

// Seen on Saike TanhAA saturator
struct allpass_interpolator {
  //----------------------------------------------------------------------------
  enum coeffs { nu, n_coeffs };
  enum state { y0, d0, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (crange<vec_value_type_t<V>> c, V frac)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    auto v = ((T) 1.0 - frac) / (1.0 + frac);
    vec_store (&c[nu * traits.size], v);
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
    crange<vec_value_type_t<V>> z, // states interleaved, ready to SIMD load
    V                           in) // N inputs ready to SIMD load
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);
    assert (z.size() >= traits.size * n_states);

    T const* nu_ptr = &c[nu * traits.size];
    T*       d0_ptr = &z[d0 * traits.size];
    T*       y0_ptr = &z[y0 * traits.size];

    V nu_v = vec_load<V> (nu_ptr);
    V d0_v = vec_load<V> (d0_ptr);
    V y0_v = vec_load<V> (y0_ptr);

    y0_v = nu_v * in + d0_v - nu_v * y0_v;

    vec_store (y0_ptr, y0_v);
    vec_store (d0_ptr, in);
    return y0_v;
  }
};
//----------------------------------------------------------------------------

} // namespace artv
