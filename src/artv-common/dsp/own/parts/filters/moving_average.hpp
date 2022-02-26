#pragma once

#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// Quick and dirty. As of now just using the variant with M/L = 2. No need to
// care about a recursive/non-recursive implementation (the recursive variant is
// very efficient but has cummulative errors). Very small sizes are used.
// Implementing as a template.
//------------------------------------------------------------------------------
template <uint L>
struct moving_average {
  static constexpr uint length = L;
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states = length - 1 };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V> z, // state 'z1' 1 to N
    V         in) // in' 1 to N
  {
    using T = vec_value_type_t<V>;
    assert (z.size() >= n_states);

    V sum  = in;
    V prev = in;
    for (uint i = 0; i < (length - 1); ++i) {
      V next = z[0];
      z[0]   = prev;
      sum += next;
      prev = next;
      z.cut_head (1);
    }
    static constexpr double coeff = (T) 1 / (T) length;
    sum *= coeff;
    return sum;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
