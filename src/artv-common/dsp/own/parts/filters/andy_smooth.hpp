#pragma once

// smoothing filter from https://cytomic.com/files/dsp/DynamicSmoothing.pdf

#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace andy {
// Algo 1
struct smoother {
  //----------------------------------------------------------------------------
  enum coeffs { g0, s, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { z1, z2, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   freq,
    V                   sensitivity,
    vec_value_type_t<V> t_spl)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    V gc   = vec_tan ((T) M_PI * freq * t_spl);
    co[g0] = (T) 2 * gc / ((T) 1 + gc);
    co[s]  = sensitivity * (T) 4;
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
  static V tick (xspan<V const> co, xspan<V> st, V in)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V z1v   = st[z1];
    V z2v   = st[z2];
    V bandz = z1v - z2v;
    V g     = vec_min (co[g0] + co[s] * vec_abs (bandz), (T) 1);
    st[z1]  = z1v + g * (in - z1v);
    st[z2]  = z2v + g * (st[z1] - z2v);
    return st[z2];
  }
  //----------------------------------------------------------------------------
};

}} // namespace artv::andy
