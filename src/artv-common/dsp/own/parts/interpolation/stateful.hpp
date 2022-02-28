#pragma once

#include <cassert>

#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/biquad.hpp"

namespace artv {

// https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolators.html

template <uint N> // order
struct thiran_interp;

template <>
struct thiran_interp<1> {
  //----------------------------------------------------------------------------
#ifdef NDEBUG
  enum coeffs {a1, n_coeffs};
#else
  enum coeffs { a1, frac, n_coeffs };
#endif
  enum coeffs_int { n_coeffs_int };
  enum state { y1, n_states };
  static constexpr uint n_points = 1;
  static constexpr uint x_offset = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V fractional)
  {
    using T     = vec_value_type_t<V>;
    V corrected = vec_max (fractional, vec_set<V> ((T) 0.000001));
    co[a1]      = (1 - corrected) / (1 + corrected);
#ifndef NDEBUG
    co[frac] = fractional;
#endif
  }
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
    crange<const V>         co,
    crange<V>               st,
    std::array<V, n_points> y,
    V                       x [[maybe_unused]])
  {
    assert (x == co[frac]);
    V ret  = co[a1] * (y[0] - st[y1]);
    st[y1] = ret;
    return ret;
  }
  //----------------------------------------------------------------------------
};

template <>
struct thiran_interp<2> {
  //----------------------------------------------------------------------------
#ifdef NDEBUG
  enum coeffs {n_coeffs = biquad::n_coeffs};
#else
  enum coeffs { frac = biquad::n_coeffs, n_coeffs };
#endif
  enum coeffs_int { n_coeffs_int };
  enum state { n_states = biquad::n_states };
  static constexpr uint n_points = 1;
  static constexpr uint x_offset = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V fractional)
  {
    biquad::reset_coeffs (co, fractional, thiran_tag {});
#ifndef NDEBUG
    co[frac] = fractional;
#endif
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    biquad::reset_states (st);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>         co,
    crange<V>               st,
    std::array<V, n_points> y,
    V                       x [[maybe_unused]])
  {
    assert (x == co[frac]);
    return biquad::tick (co, st, y[0]);
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
