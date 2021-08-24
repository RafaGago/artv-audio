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
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();
    constexpr T    pi_x2  = (T) 6.283185307179586476925286766559;

    vec_store (c, vec_exp (-pi_x2 * freq / sr)); // just on coeff
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
    V z1_v = vec_load<V> (z); // 1 coeff only

    z1_v = (in * a0_v) + (z1_v * b1_v);
    vec_store (z, z1_v); // 1 coeff only
    return z1_v;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       z, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
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
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> co, double freq, double sr)
  {
    assert (co.size() >= n_coeffs);

    // BLT(?) from ReEQ
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

    // BLT(?) from ReEQ
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

    // BLT(?) from ReEQ
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

    // BLT(?) from ReEQ
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

    // TDF II
    double out = in * co[b0] + st[s1];
    st[s1]     = in * co[b1] - out * co[a1];
    return out;
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
    double_x2 a1v = vec_set<2> (co[a1]);
    double_x2 b0v = vec_set<2> (co[b0]);
    double_x2 b1v = vec_set<2> (co[b1]);
    double_x2 s1v = {st[0][s1], st[1][s1]};

    // TDF II
    auto out  = in * b0v + s1v;
    s1v       = in * b1v - out * a1v;
    st[0][s1] = s1v[0];
    st[1][s1] = s1v[1];
    return out;
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
