#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/traits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

//------------------------------------------------------------------------------
namespace artv {

struct onepole_smoother {
  //----------------------------------------------------------------------------
  enum coeffs { b1, n_coeffs };
  enum state { z1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         srate)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    constexpr T    pi_x2  = (T) 6.283185307179586476925286766559;

    vec_store (c, vec_exp (-pi_x2 * freq / srate)); // just on coeff
  }
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
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<vec_value_type_t<V>>       z, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (z.size() >= traits.size * n_states);
    assert (c.size() >= n_coeffs);

    V a0_v = vec_set<V> (((T) 1.) - c[b1]);
    V b1_v = vec_set<V> (c[b1]);
    V z1_v = vec_load<V> (z); // 1 coeff only

    z1_v = (in * a0_v) + (z1_v * b1_v);
    vec_store (z, z1_v); // 1 coeff only
    return z1_v;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       z, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (z.size() >= traits.size * n_states);
    assert (c.size() >= n_coeffs);

    V b1v = vec_load<V> (c); // 1 coeff only
    V z1v = vec_load<V> (z); // 1 coeff only

    z1v = (in * (1. - b1v)) + (z1v * b1v);
    vec_store (z, z1v);
    return z1v;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct onepole {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, a1, n_coeffs };
  enum state { s1, n_states };

  template <class T>
  struct sub_coeffs {
    T w;
    T n;
  };
  //----------------------------------------------------------------------------
  template <class T, class U>
  static sub_coeffs<T> get_sub_coeffs (T freq, U srate)
  {
    sub_coeffs<T> r;
    if constexpr (is_vec_v<T>) {
      r.w = vec_tan ((U) M_PI * freq / srate);
    }
    else {
      r.w = tan (M_PI * freq / srate);
    }
    r.n = (U) 1. / ((U) 1. + r.w);
    return r;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    sub_coeffs<V>               wn,
    lowpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    vec_store (&co[b0 * traits.size], wn.w * wn.n);
    vec_store (&co[b1 * traits.size], wn.w * wn.n);
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    sub_coeffs<V>               wn,
    highpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    vec_store (&co[b0 * traits.size], wn.n);
    vec_store (&co[b1 * traits.size], -wn.n);
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    sub_coeffs<V>               wn,
    allpass_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    // Found empirically
    vec_store (&co[b0 * traits.size], wn.n * (wn.w - (T) 1.));
    vec_store (&co[b1 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - (T) 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate,
    lowpass_tag                 t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate,
    highpass_tag                t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate,
    allpass_tag                 t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
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
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= traits.size * n_states);

    V b0v = vec_set<V> (co[b0]);
    V b1v = vec_set<V> (co[b1]);
    V a1v = vec_set<V> (co[a1]);
    V s1v = vec_load<V> (&st[s1 * traits.size]);

    // TDF II
    auto out = in * b0v + s1v;
    s1v      = in * b1v - out * a1v;
    vec_store<V> (&st[s1 * traits.size], s1v);
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= traits.size * n_coeffs);
    assert (st.size() >= traits.size * n_states);

    V b0v = vec_load<V> (&co[b0 * traits.size]);
    V b1v = vec_load<V> (&co[b1 * traits.size]);
    V a1v = vec_load<V> (&co[a1 * traits.size]);
    V s1v = vec_load<V> (&st[s1 * traits.size]);

    // TDF II
    auto out = in * b0v + s1v;
    s1v      = in * b1v - out * a1v;
    vec_store<V> (&st[s1 * traits.size], s1v);
    return out;
  }
};
//------------------------------------------------------------------------------
} // namespace artv
