#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"

namespace artv {

// Seen on Saike TanhAA saturator
struct allpass_interpolator {
  //----------------------------------------------------------------------------
  enum coeffs { nu, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y0, d0, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V frac)
  {
    using T = vec_value_type_t<V>;

    assert (c.size() >= n_coeffs);

    c[nu] = ((T) 1.0 - frac) / (1.0 + frac);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V> c, crange<V> z, V in)
  {
    assert (c.size() >= n_coeffs);
    assert (z.size() >= n_states);

    z[y0] = c[nu] * in + z[d0] - c[nu] * z[y0];
    z[d0] = in;
    return z[y0];
  }
};
//----------------------------------------------------------------------------

} // namespace artv
