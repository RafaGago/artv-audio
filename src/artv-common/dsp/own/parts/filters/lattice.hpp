#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

//------------------------------------------------------------------------------
namespace artv {

// Nested Schroeder allpasses with unit delay. N is the number of nestings.
// https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
// It seems equivalent to a two-multiply lattice. No need to bother with the
// 1-multiply form I guess.
template <uint N>
struct lattice {
  //----------------------------------------------------------------------------
  static_assert (N > 0, "");
  enum coeffs { n_coeffs = N };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states = N };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, xspan<V> g)
  {
    assert (co.size() >= n_coeffs);
    assert (g.size() >= n_coeffs);
    xspan_memcpy (co, g.cut_head (n_coeffs));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (xspan<const V> co, xspan<V> st, V in)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    T fwd = in + st[0] * -co[0];
    T out = st[0] + fwd * co[0];
    for (uint i = 1; i < N; ++i) {
      fwd += st[i] * co[i];
      T bwd     = st[i] + fwd * -co[i];
      st[i - 1] = bwd;
    }
    st[N - 1] = fwd; // this is the all-pole output, sending the allpass out
    return out;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
