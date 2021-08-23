#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct pow2_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    return vec_abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return abs (x * x * x * (T) (1. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_abs (x * x * x * (T) (1. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return abs (x) * x * x * x * (T) (1. / 12.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_abs (x) * x * x * x * (T) (1. / 12.);
  }
  //----------------------------------------------------------------------------
};
#if 1
// TODO: is the sign preservation breaking this?
//------------------------------------------------------------------------------
// This one has a trivial simplification, it doesn't require branches.
//
// The formula is:
// ((1/6) * x**3) - ((1/6) * x1**3) / x - x1
//
// It can easily be factored to:
//
// (1/6) * (x**2 + x*x1 + x1**2)
//------------------------------------------------------------------------------
class pow2_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_pow2, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void init_states (crange<T> s)
  {}
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_states_simd (crange<vec_value_type_t<V>> s)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    T x1v      = st[x1];
    T x1_pow2v = st[x1_pow2];

    T xpow2 = x * x;

    st[x1]      = x;
    st[x1_pow2] = xpow2;

    return abs (xpow2 + (x1v * x) + x1_pow2v)
      * sgn_no_zero (x, (T) -1. / 6., (T) 1. / 6.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (st.size() >= traits.size * n_states);

    T* x1v_ptr     = &st[x1 * traits.size];
    T* x1_pow2_ptr = &st[x1_pow2 * traits.size];

    V x1v      = vec_load<V> (x1v_ptr);
    V x1_pow2v = vec_load<V> (x1_pow2_ptr);

    V xpow2 = x * x;

    vec_store (x1v_ptr, x);
    vec_store (x1_pow2_ptr, xpow2);

    return vec_abs (xpow2 + (x1v * x) + x1_pow2v)
      * vec_sgn_no_zero (
             x, vec_set<V> ((T) -1. / 6.), vec_set<V> ((T) 1. / 6.));
  }
};
//------------------------------------------------------------------------------
template <uint order>
using pow2_adaa = std::conditional_t<
  order == 1,
  pow2_adaa_1,
  adaa::waveshaper<pow2_functions, order>>;
//------------------------------------------------------------------------------
#else
template <uint order>
using pow2_adaa = adaa::waveshaper<pow2_functions, order>;
#endif
} // namespace artv
