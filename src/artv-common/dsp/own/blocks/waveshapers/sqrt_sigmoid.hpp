#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
// sqrt waveshaper with first derivative AA. Based on this thread:
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=521377&sid=18fa45082ea11269f4c29c3ccf36becd
struct sqrt_sigmoid_functions {
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return x / sqrt (x * x + T (1.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    return x / vec_sqrt (x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return sqrt (x * x + T (1.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_sqrt (x * x + (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return ((x * sqrt (x * x + T (1.))) + asinh (x)) * T (0.5);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return ((x * vec_sqrt (x * x + T (1.))) + vec_asinh (x)) * T (0.5);
  }
};
//------------------------------------------------------------------------------
// This is special cased, so it doesn't have to branch on the problematic
// division:
// https://www.kvraudio.com/forum/viewtopic.php?t=521377&start=30
// message by "martinvicanek"
// (sqrt(1 + x^2) - sqrt(1+ x1^2))/(x - x1) =
//    (x + x1)/(sqrt(1 + x^2) + sqrt(1 + x1^2))

class sqrt_sigmoid_adaa_1 {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_sqrt, n_states };
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
    T x1v     = st[x1];
    T x1vsqrt = st[x1_sqrt];
    T xsqrt   = sqrt (x * x + (T) 1);

    st[x1]      = x;
    st[x1_sqrt] = xsqrt;

    return (x + x1v) / (xsqrt + x1vsqrt);
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

    T* x1v_ptr      = &st[x1 * traits.size];
    T* x1v_sqrt_ptr = &st[x1_sqrt * traits.size];

    V x1v     = vec_load<T> (x1v_ptr);
    V x1vsqrt = vec_load<T> (x1v_sqrt_ptr);

    V xsqrt = vec_sqrt (x * x + (T) 1.);

    vec_store (x1v_ptr, x);
    vec_store (x1v_sqrt_ptr, xsqrt);

    return (x + x1v) / (xsqrt + x1vsqrt);
  }
};

#if 0
template <uint order>
using sqrt_sigmoid_adaa = std::conditional_t<
  order == 1,
  sqrt_sigmoid_adaa_1,
  adaa::waveshaper<sqrt_sigmoid_functions, order>>;
#else
template <uint order>
using sqrt_sigmoid_adaa = adaa::waveshaper<sqrt_sigmoid_functions, order>;
#endif
//------------------------------------------------------------------------------
} // namespace artv
