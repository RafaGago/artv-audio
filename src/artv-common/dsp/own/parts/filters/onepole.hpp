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
  enum coeffs_int { n_coeffs_int };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c, V freq, vec_value_type_t<V> srate)
  {
    using T              = vec_value_type_t<V>;
    constexpr auto pi_x2 = (T) (2. * M_PI);

    c[b1] = vec_exp (-pi_x2 * freq / srate);
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
    crange<const vec_value_type_t<V>> c, // coeffs (1 set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    using T = vec_value_type_t<V>;

    assert (st.size() >= n_states);
    assert (c.size() >= n_coeffs);

    V a0_v = vec_set<V> (((T) 1.) - c[b1]);
    V b1_v = vec_set<V> (c[b1]);

    st[y1] = (in * a0_v) + (st[y1] * b1_v);
    return st[y1];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> c, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
  {
    assert (st.size() >= n_states);
    assert (c.size() >= n_coeffs);

    st[y1] = (in * (1. - c[b1])) + (st[y1] * c[b1]);
    return st[y1];
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// Chapter 3.10 THE ART OF VA FILTER DESIGN
//
// Vadim Zavalishin
// https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/VAFilterDesign_2.1.0.pdf
//------------------------------------------------------------------------------
template <class... Tags>
struct onepole {
  //----------------------------------------------------------------------------
  using enabled_modes                 = mp11::mp_unique<mp11::mp_list<Tags...>>;
  static constexpr bool returns_array = mp11::mp_size<enabled_modes>::value > 1;
  //----------------------------------------------------------------------------
  enum coeffs { G, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { s, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V freq, vec_value_type_t<V> srate)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    T inv_sr_div_2 = (T) 0.5 / (srate);
    V wd           = (T) (M_PI * 2.) * freq;
    V wa           = (T) 2. * srate * vec_tan (wd * inv_sr_div_2);
    V g            = wa * inv_sr_div_2;
    co[G]          = g / ((T) 1.0 + g);
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
  static auto tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               x)
  {
    assert (co.size() >= n_coeffs);
    return tick (co[G], st, x);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 x)
  {
    assert (co.size() >= n_coeffs);
    return tick (vec_set<V> (co[G]), st, x);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (V G_v, crange<V> st, V x)
  {
    constexpr auto traits = vec_traits<V>();

    V lp, hp;

    V v   = (x - st[s]) * G_v;
    lp    = v + st[s];
    st[s] = lp + v;

    std::array<V, mp11::mp_size<enabled_modes>::value> ret;

    mp_foreach_idx (enabled_modes {}, [&] (auto index, auto mode) {
      using tag = decltype (mode);

      if constexpr (std::is_same_v<tag, lowpass_tag>) {
        ret[index.value] = lp;
      }
      else if constexpr (std::is_same_v<tag, highpass_tag>) {
        hp               = x - lp;
        ret[index.value] = hp;
      }
      else if constexpr (std::is_same_v<tag, allpass_tag>) {
        // The optimizer should trivially detect hp being calculated twice when
        // both allpass and higpass are enabled. Avoiding throwing templates to
        // the problem.
        hp               = x - lp;
        ret[index.value] = lp - hp;
      }
    });

    if constexpr (returns_array) {
      return ret;
    }
    else {
      return ret[0];
    }
  }
};
//------------------------------------------------------------------------------
using onepole_lowpass  = onepole<lowpass_tag>;
using onepole_highpass = onepole<highpass_tag>;
using onepole_allpass  = onepole<allpass_tag>;
//------------------------------------------------------------------------------
// When possible prfer the TPT variant, it only has 1 coeff and 1 state.
struct onepole_tdf2 {
  //----------------------------------------------------------------------------
  enum coeffs { b0, b1, a1, n_coeffs };
  enum coeffs_int { n_coeffs_int };
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
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, lowpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    co[b0] = wn.w * wn.n;
    co[b1] = co[b0];
    co[a1] = wn.n * (wn.w - (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, highpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    co[b0] = wn.n;
    co[b1] = -wn.n;
    co[a1] = wn.n * (wn.w - (T) 1.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, sub_coeffs<V> wn, allpass_tag)
  {
    using T = vec_value_type_t<V>;
    assert (co.size() >= n_coeffs);

    // Found empirically
    co[b0] = wn.n * (wn.w - (T) 1.);
    co[b1] = vec_set<V> ((T) 1.);
    co[a1] = co[b0];
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> srate,
    lowpass_tag         t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> srate,
    highpass_tag        t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co,
    V                   freq,
    vec_value_type_t<V> srate,
    allpass_tag         t)
  {
    reset_coeffs (co, get_sub_coeffs (freq, srate), t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V> co,
    V         d, // 0 to 1
    thiran_tag)
  {
    using T = vec_value_type_t<V>;

    // according to this. Table 3.2:
    // http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt3_allpass.pdf
    // The point minimizing the average error for D0 (D0 < D < D0 +1) is 0.418
    // but it broke when modulating.

    d += (T) 1; // d = delta
    d      = vec_max ((T) 1.001, d); // finite precision correction
    co[a1] = ((T) 1 - d) / ((T) 1 + d);
    co[b0] = co[a1];
    co[b1] = vec_set<V> ((T) 0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V a1v, V b0v, V b1v, raw_tag)
  {
    co[a1] = a1v;
    co[b0] = b0v;
    co[b1] = b1v;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    memset (st.data(), 0, sizeof (V) * n_states);
  }
  //----------------------------------------------------------------------------
  // 1 set of coeffs, N outs. (E.g. stereo filter using double). Interleaved
  // states version
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (1 set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    V b0v = vec_set<V> (co[b0]);
    V b1v = vec_set<V> (co[b1]);
    V a1v = vec_set<V> (co[a1]);

    auto out = in * b0v + st[s1];
    st[s1]   = in * b1v - out * a1v;
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs);
    assert (st.size() >= n_states);

    auto out = in * co[b0] + st[s1];
    st[s1]   = in * co[b1] - out * co[a1];
    return out;
  }
};
//------------------------------------------------------------------------------
} // namespace artv
