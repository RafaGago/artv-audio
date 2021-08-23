#pragma once

/* Fiters from: https://cytomic.com/files/dsp/SvfLinearTrapOptimised2.pdf*/

#include <cmath>
#include <gcem.hpp>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
namespace andy {

struct svf {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, m0, m1, m2, n_coeffs };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 0.;
    *co.m1   = 0.;
    *co.m2   = 1.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    V g = vec_tan (freq * ((T) M_PI) / sr);
    V k = ((T) 1.0) / q;

    V a1_v = (T) 1. / (g * (g + k) + ((T) 1.));
    vec_store (&c[a1 * traits.size], a1_v);

    V a2_v = g * a1_v;
    vec_store (&c[a2 * traits.size], a2_v);

    V a3_v = g * a2_v;
    vec_store (&c[a3 * traits.size], a3_v);

    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 0.));
    vec_store (&c[m1 * traits.size], vec_set<V> ((T) 0.));
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 1.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void bandpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 0.;
    *co.m1   = 1.;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void highpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = -1.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    V g = vec_tan ((T) M_PI * freq / sr);
    V k = (T) 1.0 / q;

    V a1_v = (T) 1.0 / ((T) 1.0 + g * (g + k));
    vec_store (&c[a1 * traits.size], a1_v);

    V a2_v = g * a1_v;
    vec_store (&c[a2 * traits.size], a2_v);

    V a3_v = g * a2_v;
    vec_store (&c[a3 * traits.size], a3_v);

    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], -k);
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) -1.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void notch (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void peak (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -k;
    *co.m2   = -2.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void allpass (
    crange<T> c,
    T         freq,
    T         q,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = -2. * k;
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void allpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    V g = vec_tan (freq * ((T) M_PI) / sr);
    V k = ((T) 1.0) / q;

    V a1_v = (T) 1.0 / (g * (g + k) + ((T) 1.));
    vec_store (&c[a1 * traits.size], a1_v);

    V a2_v = g * a1_v;
    vec_store (&c[a2 * traits.size], a2_v);

    V a3_v = g * a2_v;
    vec_store (&c[a3 * traits.size], a3_v);

    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], -2. * k);
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 0.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void bell (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      interleaved_pack_offset = 0,
    uint      interleaved_pack_size   = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double A = exp (db * (T) (1. / 40.) * (T) M_LN10);
    double g = tan (M_PI * freq / sr);
    double k = 1.0 / (q * A);
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = k * (A * A - 1.);
    *co.m2   = 0.;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void bell_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    V                           db,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    V    A = vec_exp (db * (T) (1. / 40.) * (T) M_LN10);
    auto g = vec_tan ((freq * (T) M_PI) / (T) sr);
    auto k = (T) 1. / (q * A);

    V a1_v = ((T) 1.0) / (g * (g + k) + ((T) 1.));
    vec_store (&c[a1 * traits.size], a1_v);

    V a2_v = g * a1_v;
    vec_store (&c[a2 * traits.size], a2_v);

    V a3_v = g * a2_v;
    vec_store (&c[a3 * traits.size], a3_v);

    vec_store (&c[m0 * traits.size], vec_set<V> ((T) 1.));
    vec_store (&c[m1 * traits.size], k * (A * A - (T) 1.));
    vec_store (&c[m2 * traits.size], vec_set<V> ((T) 0.));
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void low_shelf (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (coeff_offset * coeff_distance * n_coeffs);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    double A = exp (db * (T) (1. / 40.) * (T) M_LN10);
    double g = tan (M_PI * freq / sr) / sqrt (A);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = 1.;
    *co.m1   = k * (A - 1.);
    *co.m2   = (A * A) - 1.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void high_shelf (
    crange<T> c,
    T         freq,
    T         q,
    T         db,
    T         sr,
    uint      coeff_offset   = 0,
    uint      coeff_distance = 1)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (coeff_offset * coeff_distance * n_coeffs);
    auto co = unpack_interleaved_coeffs (c, coeff_distance);

    constexpr T ln10 = (T) gcem::log (10.);

    double A = exp (db * (T) (1. / 40.) * (T) M_LN10);
    double g = tan (M_PI * freq / sr) * sqrt (A);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
    *co.m0   = A * A;
    *co.m1   = k * (1. - A) * A;
    *co.m2   = 1. - (A * A);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {
    static_assert (std::is_floating_point<T>::value, "");
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void reset_states (crange<T> s)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               v0)
  {
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    return calc (
      s[ic1eq], s[ic2eq], c[m0], c[m1], c[m2], c[a1], c[a2], c[a3], v0);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Sparse
  // states version.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>>                              c, // coeffs
    std::array<crange<vec_value_type_t<V>>, vec_traits_t<V>::size> s, // states
    V                                                              v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    for (uint i = 0; i < traits.size; ++i) {
      assert (s[i].size() >= n_states);
      assert (s[i].size() >= n_states);
    }

    V ic1eq_v, ic2eq_v;

    for (uint i = 0; i < traits.size; ++i) {
      ic1eq_v[i] = s[i][ic1eq];
      ic2eq_v[i] = s[i][ic2eq];
    }

    V out = calc<V> (
      ic1eq_v,
      ic2eq_v,
      vec_set<V> (c[m0]),
      vec_set<V> (c[m1]),
      vec_set<V> (c[m2]),
      vec_set<V> (c[a1]),
      vec_set<V> (c[a2]),
      vec_set<V> (c[a3]),
      v0);

    for (uint i = 0; i < traits.size; ++i) {
      s[i][ic1eq] = ic1eq_v[i];
      s[i][ic2eq] = ic2eq_v[i];
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_aligned (
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= n_coeffs);
    assert (s.size() >= traits.size * n_states);

    V ic1eq_v = vec_load<V> (&s[ic1eq * traits.size]);
    V ic2eq_v = vec_load<V> (&s[ic2eq * traits.size]);

    V out = calc<V> (
      ic1eq_v,
      ic2eq_v,
      vec_set<V> (c[m0]),
      vec_set<V> (c[m1]),
      vec_set<V> (c[m2]),
      vec_set<V> (c[a1]),
      vec_set<V> (c[a2]),
      vec_set<V> (c[a3]),
      v0);

    vec_store (&s[ic1eq * traits.size], ic1eq_v);
    vec_store (&s[ic2eq * traits.size], ic2eq_v);

    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);
    assert (s.size() >= traits.size * n_states);

    V ic1eq_v = vec_load<V> (&s[ic1eq * traits.size]);
    V ic2eq_v = vec_load<V> (&s[ic2eq * traits.size]);

    V out = calc<V> (
      ic1eq_v,
      ic2eq_v,
      vec_load<V> (&c[m0 * traits.size]),
      vec_load<V> (&c[m1 * traits.size]),
      vec_load<V> (&c[m2 * traits.size]),
      vec_load<V> (&c[a1 * traits.size]),
      vec_load<V> (&c[a2 * traits.size]),
      vec_load<V> (&c[a3 * traits.size]),
      v0);

    vec_store (&s[ic1eq * traits.size], ic1eq_v);
    vec_store (&s[ic2eq * traits.size], ic2eq_v);

    return out;
  }
  //----------------------------------------------------------------------------
private:
  template <class T>
  static T calc (
    T& ic1eq_v,
    T& ic2eq_v,
    T  m0_v,
    T  m1_v,
    T  m2_v,
    T  a1_v,
    T  a2_v,
    T  a3_v,
    T  v0)
  {
    T v3    = v0 - ic2eq_v;
    T v1    = a1_v * ic1eq_v + a2_v * v3;
    T v2    = ic2eq_v + a2_v * ic1eq_v + a3_v * v3;
    T out   = m0_v * v0 + m1_v * v1 + m2_v * v2;
    ic1eq_v = (v1 * 2.) - ic1eq_v;
    ic2eq_v = (v2 * 2.) - ic2eq_v;
    return out;
  }
  //----------------------------------------------------------------------------
  template <class T>
  struct coeff_ptrs {
    T* a1;
    T* a2;
    T* a3;
    T* m0;
    T* m1;
    T* m2;
  };
  //----------------------------------------------------------------------------
  template <class T>
  static coeff_ptrs<T> unpack_interleaved_coeffs (crange<T> c, uint distance)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (c.size() >= ((n_coeffs - 1) * distance) + 1);

    coeff_ptrs<T> ret;
    ret.a1 = &c[a1 * distance];
    ret.a2 = &c[a2 * distance];
    ret.a3 = &c[a3 * distance];
    ret.m0 = &c[m0 * distance];
    ret.m1 = &c[m1 * distance];
    ret.m2 = &c[m2 * distance];
    return ret;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct svf_lowpass {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, a3, n_coeffs };
  enum state { ic1eq, ic2eq, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (crange<T> c, T freq, T q, T sr)
  {
    static_assert (std::is_floating_point<T>::value, "");
    c       = c.shrink_head (interleaved_pack_offset);
    auto co = unpack_interleaved_coeffs (c, interleaved_pack_size);

    double g = tan (M_PI * freq / sr);
    double k = 1.0 / q;
    *co.a1   = 1.0 / (1.0 + g * (g + k));
    *co.a2   = g * *co.a1;
    *co.a3   = g * *co.a2;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> c,
    V                           freq,
    V                           q,
    vec_value_type_t<V>         sr)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= (n_coeffs * traits.size));

    V g = vec_tan (freq * ((T) M_PI) / sr);
    V k = ((T) 1.0) / q;

    V a1_v = (T) 1. / (g * (g + k) + ((T) 1.));
    vec_store (&c[a1 * traits.size], a1_v);

    V a2_v = g * a1_v;
    vec_store (&c[a2 * traits.size], a2_v);

    V a3_v = g * a2_v;
    vec_store (&c[a3 * traits.size], a3_v);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T tick (
    crange<const T> c, // coeffs
    crange<T>       s, // state
    T               v0)
  {
    assert (s.size() >= n_states);
    assert (c.size() >= n_coeffs);

    T v3     = v0 - s[ic2eq];
    T v1     = c[a1] * s[ic1eq] + c[a2] * v3;
    T v2     = s[ic2eq] + c[a2] * s[ic1eq] + c[a3] * v3;
    s[ic1eq] = (v1 * 2.) - s[ic1eq];
    s[ic2eq] = (v2 * 2.) - s[ic2eq];
    return v2;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> c, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       s, // states (interleaved, SIMD aligned)
    V                                 v0)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (c.size() >= traits.size * n_coeffs);
    assert (s.size() >= traits.size * n_states);

    V ic1eq_v = vec_load<V> (&s[ic1eq * traits.size]);
    V ic2eq_v = vec_load<V> (&s[ic2eq * traits.size]);
    V a1_v    = vec_load<V> (&c[a1 * traits.size]);
    V a2_v    = vec_load<V> (&c[a2 * traits.size]);
    V a3_v    = vec_load<V> (&c[a3 * traits.size]);

    V v3    = v0 - ic2eq_v;
    V v1    = a1_v * ic1eq_v + a2_v * v3;
    V v2    = ic2eq_v + a2_v * ic1eq_v + a3_v * v3;
    ic1eq_v = (v1 * 2.) - ic1eq_v;
    ic2eq_v = (v2 * 2.) - ic2eq_v;

    vec_store (&s[ic1eq * traits.size], ic1eq_v);
    vec_store (&s[ic2eq * traits.size], ic2eq_v);

    return v2;
  }
  //----------------------------------------------------------------------------
}
} // namespace artv::andy
