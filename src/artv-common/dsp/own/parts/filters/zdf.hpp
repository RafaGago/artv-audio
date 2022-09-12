#pragma once
// routines for solving feedbacks with zero delay

#include <cassert>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv { namespace zdf {

// tag to select overloads giving back the G and S zdf coefficients
// See "The ART ov VA filter design" (Vadim Zabalishin) for what G and S are.
struct gs_coeffs_tag {};

template <class T>
struct response {
  T G;
  T S;
};

static constexpr uint n_gs_coeffs = 2;
//------------------------------------------------------------------------------
// Accumulates the linear response of a number of elements connected in series.
// G_S contains a list of G and S values
//
// See "The ART ov VA filter design" (Vadim Zabalishin) 5.3
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
struct response<V> combine_response (crange<const V> G_S) {
  using T = vec_value_type_t<V>;
  assert ((G_S.size() % n_gs_coeffs) == 0);

  response<V> ret;
  ret.G = vec_set<V> ((T) 1);
  ret.S = vec_set<V> ((T) 0);

  for (uint i = 0; i < G_S.size(); i += n_gs_coeffs) {
    ret.G *= G_S[i];
    ret.S *= G_S[i];
    ret.S += G_S[i + 1];
  }
  return ret;
}
//------------------------------------------------------------------------------
// topologies
struct linear_tag {};
struct sqrt_sig_after_fb_juction_tag {};
struct sqrt_sig_before_fb_juction_tag {};
struct sqrt_sig_before_fb_juction_pp_tag {};

// numeric methods

// See Urs Heckmann and Mystran's comments:
//  https://www.kvraudio.com/forum/viewtopic.php?start=375&t=350246&sid=eb35460aab9d1bda3038aa78018b5182
struct lin_mystran_2_tag {};

// TO Add/test
// struct lin_mystran_tag {};
// struct newton_rhapson_tag {};
// struct euler_tag {};
// struct heun_tag {};
// struct runge_kutta_tag {};

//------------------------------------------------------------------------------
template <class Topology_tag, class Numerical_method_tag>
struct feedback;
//------------------------------------------------------------------------------
// Solve zdf feedback loop. Linear case. "k" is the feedback ratio.
// See "The ART ov VA filter design" (Vadim Zabalishin) 5.3
// Topology:
//          __________
//       u  |        |
// x -------| G*_x+S |------ y
//    |-    |________|     |
//    |      ________      |
//    |     |        |     |
//    ------|  * k   |------
//          |________|
//
// where:
// u = x - k * (G * u + S)

template <class Any>
struct feedback<linear_tag, Any> {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (crange<const V>, crange<V>, V in, response<V> resp, V k)
  {
    using T = vec_value_type_t<V>;
    return (in - k * resp.S) / (k * resp.G + (T) 1);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// global ZDF feedback loop with "x/sqrt(x**2*h + 1)" sigmoid with variable
// hardness (h parameter) solved/approximated with 2 rounds Mystran's
// linearization.
//
// Topology:
//          ________  __________
//       u |        | |        |
// x ------|  f(x)  |-| G*_x+S |------ y
//    |-   |________| |________|     |
//    |          ________            |
//    |         |        |           |
//    ----------|  * k   |-  ---------
//              |________|
//
// where:
// u = x - k * (G * linramp * u + S)

template <>
class feedback<sqrt_sig_after_fb_juction_tag, lin_mystran_2_tag> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    st[lin] = vec_set<V> (1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> resp,
    V           k,
    V           h)
  {
    using T = vec_value_type_t<V>;
    // 1st round
    auto scaled = resp;
    scaled.G *= st[lin];
    V u = feedback<linear_tag, void>::tick<V> ({}, {}, scaled, k, in);
    if (vec_is_all_ones (h == vec_set<V> (0))) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    st[lin] = (T) 1 / vec_sqrt (u * u * h + (T) 1); // sqrt_sigmoid (u) / u
    // 2 st round (compensation)
    scaled = resp;
    scaled.G *= st[lin];
    u += feedback<linear_tag, void>::tick<V> ({}, {}, scaled, k, in);
    u *= (T) 0.5;
    st[lin] = (T) 1 / vec_sqrt (u * u * h + (T) 1); // sqrt_sigmoid (u) / u
    return u;
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// global ZDF feedback loop with "x/sqrt(x**2*h + 1)" sigmoid with variable
// hardness (h parameter) solved/approximated with 2 rounds Mystran's
// linearization.
//
// Topology:
//         _________
//      u |        |
// x -----| G*_x+S |-------------- y
//    |-  |________|             |
//    |    ________    ________  |
//    |   |        |  |        | |
//    ----|  f(x)  |--|  * k   |--
//        |________|  |________|
//
// where:
// u = x - linramp * k * (Gu + S)

template <>
class feedback<sqrt_sig_before_fb_juction_tag, lin_mystran_2_tag> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    st[lin] = vec_set<V> (1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> resp,
    V           k,
    V           h)
  {
    using T = vec_value_type_t<V>;
    // 1st round
    auto scaled = resp;
    scaled.G *= st[lin];
    scaled.S *= st[lin];
    V u = feedback<linear_tag, void>::tick<V> ({}, {}, scaled, k, in);
    if (vec_is_all_ones (h == vec_set<V> (0))) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    V sig_in = k * (resp.G * u + resp.S);
    // sqrt_sigmoid (u) / u
    st[lin] = (T) 1 / vec_sqrt (sig_in * sig_in * h + (T) 1);
    // 2 st round (compensation)
    scaled = resp;
    scaled.G *= st[lin];
    scaled.S *= st[lin];
    u += feedback<linear_tag, void>::tick<V> ({}, {}, scaled, k, in);
    u *= (T) 0.5;
    sig_in = k * (resp.G * u + resp.S);
    // sqrt_sigmoid (u) / u
    st[lin] = (T) 1 / vec_sqrt (sig_in * sig_in * h + (T) 1);
    return u;
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
// global ZDF feedback loop with "x/sqrt(x**2*h + 1)" sigmoid with variable
// hardness (h parameter) solved/approximated with 2 rounds Mystran's
// linearization. This one allows additional filtering (e.g DC blocker) after
// the nonlinearity.
//
// Topology:
//         __________
//      u |          |
// x -----| G1*_x+S1 |---------------------------- y
//    |-  |__________|                          |
//    |    __________   ________    ________    |
//    |   |          |  |        |  |        |  |
//    ----| G2*_x+S2 |--|  f(x)  |--|  * k   |--
//        |__________|  |________|  |________|
//
// where:
// u = x - (G2 (linramp * k * (G1 u + S1)) + S2)
// u = (-G2 S1 k linramp - S2 + in) / (G1 G2 k linramp + 1)

template <>
class feedback<sqrt_sig_before_fb_juction_pp_tag, lin_mystran_2_tag> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> c)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    assert (st.size() >= n_states);
    st[lin] = vec_set<V> (1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> r1,
    response<V> r2,
    V           k,
    V           h)
  {
    using T = vec_value_type_t<V>;
    // 1st round
    V u = -r2.G * r1.S * k * st[lin] - r2.S + in;
    u /= r1.G * r2.G * k * st[lin] + (T) 1;

    if (vec_is_all_ones (h == vec_set<V> (0))) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    V sig_in = k * (r1.G * u + r1.S);
    // sqrt_sigmoid (u) / u
    st[lin] = (T) 1 / vec_sqrt (sig_in * sig_in * h + (T) 1);
    // 2 st round (compensation)
    V u2 = -r2.G * r1.S * k * st[lin] - r2.S + in;
    u2 /= r1.G * r2.G * k * st[lin] + (T) 1;
    u      = (u + u2) * (T) 0.5;
    sig_in = k * (r1.G * u + r1.S);
    // sqrt_sigmoid (u) / u
    st[lin] = (T) 1 / vec_sqrt (sig_in * sig_in * h + (T) 1);
    return u;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

}} // namespace artv::zdf
