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
  static void lowpass (crange<T> c, T freq, T srate)
  {
    static_assert (std::is_floating_point<T>::value, "");
    constexpr double pi_x2 = 6.283185307179586476925286766559;
    c[b1]                  = (T) exp (-pi_x2 * freq / srate);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    vec_value_type_t<V>         srate)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();
    constexpr T    pi_x2  = (T) 6.283185307179586476925286766559;

    vec_store (c, vec_exp (-pi_x2 * freq / srate)); // just on coeff
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
  static void lowpass (crange<double> co, sub_coeffs<double> wn)
  {
    assert (co.size() >= n_coeffs);

    // BLT(?) from ReEQ
    co[b0] = wn.w * wn.n;
    co[b1] = wn.w * wn.n;
    co[a1] = wn.n * (wn.w - 1.);
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, sub_coeffs<double> wn)
  {
    assert (co.size() >= n_coeffs);

    // BLT(?) from ReEQ
    co[b0] = wn.n;
    co[b1] = -wn.n;
    co[a1] = wn.n * (wn.w - 1.);
  }
  //----------------------------------------------------------------------------
  static void allpass (crange<double> co, sub_coeffs<double> wn)
  {
    assert (co.size() >= n_coeffs);

    // Found empirically
    co[b0] = wn.n * (wn.w - 1.);
    co[b1] = 1.;
    co[a1] = co[b0];
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (crange<vec_value_type_t<V>> co, sub_coeffs<V> wn)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    vec_store (&co[b0 * traits.size], wn.w * wn.n);
    vec_store (&co[b1 * traits.size], wn.w * wn.n);
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (crange<vec_value_type_t<V>> co, sub_coeffs<V> wn)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    vec_store (&co[b0 * traits.size], wn.n);
    vec_store (&co[b1 * traits.size], -wn.n);
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - 1.));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void allpass_simd (crange<vec_value_type_t<V>> co, sub_coeffs<V> wn)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));

    // Found empirically
    vec_store (&co[b0 * traits.size], wn.n * (wn.w - (T) 1.));
    vec_store (&co[b1 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&co[a1 * traits.size], wn.n * (wn.w - (T) 1.));
  }
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> co, double freq, double srate)
  {
    lowpass (co, get_sub_coeffs (freq, srate));
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double srate)
  {
    highpass (co, get_sub_coeffs (freq, srate));
  }
  //----------------------------------------------------------------------------
  static void allpass (crange<double> co, double freq, double srate)
  {
    allpass (co, get_sub_coeffs (freq, srate));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate)
  {
    lowpass_simd (co, get_sub_coeffs (freq, srate));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate)
  {
    highpass_simd (co, get_sub_coeffs (freq, srate));
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void allpass_simd (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         srate)
  {
    allpass_simd (co, get_sub_coeffs (freq, srate));
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
    auto out = in * b0v + s1v;
    s1v      = in * b1v - out * a1v;

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
