#pragma once

#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

// Quick and dirty. As of now just using the variant with M/L = 2. No need to
// care about a recursive/non-recursive implementation (the recursive variant is
// very efficient but has cummulative errors. Implementing as a template.
//------------------------------------------------------------------------------
namespace artv {

// Aka box/boxcar
template <uint L>
struct moving_average;

template <>
struct moving_average<2> {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum state { z1, n_states };
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
    crange<vec_value_type_t<V>> z, // state 'z1' 1 to N
    V                           in) // in' 1 to N
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (z.size() >= (traits.size * n_states));

    auto z1_ptr = &z[z1 * traits.size];
    V    ret    = vec_load<V> (z1_ptr);
    ret += in;
    ret *= (T) 0.5;
    vec_store (z1_ptr, in);
    return ret;
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} // namespace artv
