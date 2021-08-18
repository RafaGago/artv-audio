#pragma once

#include <cmath>
#include <type_traits>

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
  template <class T>
  static void lowpass (crange<T> c, T freq, T sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
    constexpr double pi_x2 = 6.283185307179586476925286766559;
    c[b1]                  = (T) exp (-pi_x2 * freq / sr);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       z, // state
    T               in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (z.size() >= n_states);
    assert (c.size() >= n_coeffs);

    z[z1] = (in * (1. - c[b1])) + (z[z1] * c[b1]);
    return z[z1];
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_aligned (
    crange<const vec_value_type_t<V>> c, // coeffs, just 'b1'
    crange<vec_value_type_t<V>>       z, // state 'z1' 1 to N
    V                                 in) // in' 1 to N
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (z.size() >= traits.size * n_states);
    assert (c.size() >= n_coeffs);

    V a0_v = vec_set<V> (((T) 1.) - c[b1]);
    V b1_v = vec_set<V> (c[b1]);
    V z1_v = vec_load<V> (z);

    z1_v = (in * a0_v) + (z1_v * b1_v);
    vec_store (z, z1_v);
    return z1_v;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// From ReEQ.
struct onepole {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, a1, n_coeffs };
  enum state { z1, z0, n_states };
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> co, double freq, double sr)
  {
    assert (co.size() >= n_coeffs);

    double w = tan (M_PI * freq / sr);
    double n = 1. / (1. + w);
    co[b0]   = w * n;
    co[b1]   = w * n;
    co[a1]   = n * (w - 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    auto w = vec_tan (M_PI * freq / sr);
    auto n = 1. / (1. + w);
    vec_store (&co[b0 * traits.size], w * n);
    vec_store (&co[b1 * traits.size], w * n);
    vec_store (&co[a1 * traits.size], n * (w - 1.));
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double sr)
  {
    assert (co.size() >= n_coeffs);

    double w = tan (M_PI * freq / sr);
    double n = 1. / (1. + w);
    co[b0]   = n;
    co[b1]   = -n;
    co[a1]   = n * (w - 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    auto w = vec_tan (M_PI * freq / sr);
    auto n = 1. / (1. + w);
    vec_store (&co[b0 * traits.size], n);
    vec_store (&co[b1 * traits.size], -n);
    vec_store (&co[a1 * traits.size], n * (w - 1.));
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void reset_states (crange<double> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    st[z1] = in * co[b0] + st[z0] * co[b1] - st[z1] * co[a1];
    st[z0] = in;
    return st[z1];
  }
  //----------------------------------------------------------------------------
  static double_x2 tick (
    crange<const double>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    double_x2                     in)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    // I don't know if it's worth to use SIMD for this actually.
    double_x2 a1_v = vec_set<2> (co[a1]);
    double_x2 b0_v = vec_set<2> (co[b0]);
    double_x2 b1_v = vec_set<2> (co[b1]);
    double_x2 z0_v = {st[0][z0], st[1][z0]};
    double_x2 z1_v = {st[0][z1], st[1][z1]};

    z1_v      = (in * b0_v) + (z0_v * b1_v) - (z1_v * a1_v);
    st[0][z1] = z1_v[0];
    st[0][z0] = in[0];
    st[1][z1] = z1_v[1];
    st[1][z0] = in[1];
    return z1_v;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= traits.size * n_coeffs);
    assert (st.size() >= traits.size * n_states);

    V b0v = vec_load<V> (&co[b0 * traits.size]);
    V b1v = vec_load<V> (&co[b1 * traits.size]);
    V a1v = vec_load<V> (&co[a1 * traits.size]);
    V z0v = vec_load<V> (&st[z0 * traits.size]);
    V z1v = vec_load<V> (&st[z1 * traits.size]);

    z1v = in * b0v + z0v * b1v - z1v * a1v;
    z0v = in;

    vec_store<V> (&st[z0 * traits.size], z0v);
    vec_store<V> (&st[z1 * traits.size], z1v);
    return z1v;
  }
};
//------------------------------------------------------------------------------
} // namespace artv
