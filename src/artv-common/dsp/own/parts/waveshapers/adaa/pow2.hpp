#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

// clang-format on
struct pow2_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V fn (V x)
  {
    return vec_abs (x) * x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_abs (x * x * x * (T) (1. / 3.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_abs (x) * x * x * x * (T) (1. / 12.);
  }
  //----------------------------------------------------------------------------
};
#if 0
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
  enum coeffs_int { n_coeffs_int };
  enum state { x1, x1_pow2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (
    xspan<const V>,
    xspan<V> st,
    V                           x)
  {
    using T               = vec_value_type_t<V>;

    assert (st.size() >= n_states);

    T* x1v_ptr     = &st[x1];
    T* x1_pow2_ptr = &st[x1_pow2];

    V x1v      = st[x1];
    V x1_pow2v = st[x1_pow2];

    V xpow2 = x * x;

    st[x1] = x;
    st[x1_pow2] = xpow2;

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
