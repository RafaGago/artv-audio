#pragma once

#include "artv-common/misc/util.hpp"

namespace artv {

// Seen on Saike TanhAA saturator
struct allpass_interpolator {
  //----------------------------------------------------------------------------
  enum coeffs { nu, n_coeffs };
  enum state { y0, d0, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void init (crange<T> c, T frac)
  {
    c[nu] = (1.0 - frac) / (1.0 + frac);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (crange<vec_value_type_t<V>> c, V frac)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    auto v = ((T) 1.0 - frac) / (1.0 + frac);
    vec_store (&c[nu * traits.size], v);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr T tick (
    crange<const T> c, // coeffs
    crange<T>       z, // state
    T               in)
  {
    z[y0] = c[nu] * in + z[d0] - c[nu] * z[y0];
    z[d0] = in;
    return z[y0];
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>>
                                c, // coeffs interleaved, ready to SIMD load
    crange<vec_value_type_t<V>> z, // states interleaved, ready to SIMD load
    V                           in) // N inputs ready to SIMD load
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
