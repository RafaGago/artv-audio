#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct pow2_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    return vec_abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_abs (x * x * x * (T) (1. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>>)
  {}
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>,
    crange<vec_value_type_t<V>> st,
    V                           x)
  {
    using T               = vec_value_type_t<V>;
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
