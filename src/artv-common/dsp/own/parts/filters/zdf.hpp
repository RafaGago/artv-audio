#pragma once
// routines for solving feedbacks with zero delay

#include <cassert>
#include <limits>
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
struct nonlin_post_fb_node_tag {};
struct nonlin_pre_fb_node_tag {};
struct lin_pre_fb_node_nonlin_after_tag {};

// numeric methods

// See Urs Heckmann and Mystran's comments:
//  https://www.kvraudio.com/forum/viewtopic.php?start=375&t=350246&sid=eb35460aab9d1bda3038aa78018b5182
template <uint n_iters> // 1 or 2
struct lin_mystran_tag {};

// TO Add/test
// struct lin_mystran_tag {};
// struct newton_rhapson_tag {};
// struct euler_tag {};
// struct heun_tag {};
// struct runge_kutta_tag {};
// etc...

//------------------------------------------------------------------------------
// - Topology_tag, one of the topology tags above
// - Numerical_method_tag one of the numerical method tags above (except for
//   linear_tag)
// - Nonlin: Nonlinearity class, e.g. one of the classes on
//  "src/artv-common/dsp/own/parts/waveshapers/sigmoid.hpp" (except for
//  linear_tag)
template <
  class Topology_tag,
  class Numerical_method_tag = void,
  class Nonlin               = void>
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

template <class Any1, class Any2>
struct feedback<linear_tag, Any1, Any2> {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { n_states };
  //----------------------------------------------------------------------------
  using nonlin = Any2;
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

template <class Nonlin, uint n_iters>
class feedback<nonlin_post_fb_node_tag, lin_mystran_tag<n_iters>, Nonlin> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  using nonlin = Nonlin;
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
  template <class V, class... Ts, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> resp,
    V           k,
    Ts&&... nonlin_args)
  {
    static_assert (n_iters == 1 || n_iters == 2);
    using T = vec_value_type_t<V>;
    // 1st round
    auto scaled = resp;
    scaled.G *= st[lin];
    V u = feedback<linear_tag>::tick<V> ({}, {}, scaled, k, in);
    if (nonlin::is_linear (std::forward<Ts> (nonlin_args)...)) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    st[lin] = nonlin::tick_div_in (u, std::forward<Ts> (nonlin_args)...);
    if constexpr (n_iters == 1) {
      return u;
    }
    // 2 st round (compensation)
    scaled = resp;
    scaled.G *= st[lin];
    u += feedback<linear_tag>::tick<V> ({}, {}, scaled, k, in);
    u *= (T) 0.5;
    st[lin] = nonlin::tick_div_in (u, std::forward<Ts> (nonlin_args)...);
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

template <class Nonlin, uint n_iters>
class feedback<nonlin_pre_fb_node_tag, lin_mystran_tag<n_iters>, Nonlin> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  using nonlin = Nonlin;
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
  template <class V, class... Ts, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> resp,
    V           k,
    Ts&&... nonlin_args)
  {
    static_assert (n_iters == 1 || n_iters == 2);
    using T = vec_value_type_t<V>;
    // 1st round
    auto scaled = resp;
    scaled.G *= st[lin];
    scaled.S *= st[lin];
    V u = feedback<linear_tag>::tick<V> ({}, {}, scaled, k, in);
    if (nonlin::is_linear (std::forward<Ts> (nonlin_args)...)) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    V sig_in = k * (resp.G * u + resp.S);
    st[lin]  = nonlin::tick_div_in (sig_in, std::forward<Ts> (nonlin_args)...);
    if constexpr (n_iters == 1) {
      return u;
    }
    // 2 st round (compensation)
    scaled = resp;
    scaled.G *= st[lin];
    scaled.S *= st[lin];
    u += feedback<linear_tag>::tick<V> ({}, {}, scaled, k, in);
    u *= (T) 0.5;
    sig_in  = k * (resp.G * u + resp.S);
    st[lin] = nonlin::tick_div_in (sig_in, std::forward<Ts> (nonlin_args)...);
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

template <class Nonlin, uint n_iters>
class feedback<
  lin_pre_fb_node_nonlin_after_tag,
  lin_mystran_tag<n_iters>,
  Nonlin> {
public:
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { lin, n_states };
  //----------------------------------------------------------------------------
  using nonlin = Nonlin;
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
  template <class V, class... Ts, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>,
    crange<V>   st,
    V           in,
    response<V> r1,
    response<V> r2,
    V           k,
    Ts&&... nonlin_args)
  {
    static_assert (n_iters == 1 || n_iters == 2);
    using T = vec_value_type_t<V>;
    // 1st round
    V u = -r2.G * r1.S * k * st[lin] - r2.S + in;
    u /= r1.G * r2.G * k * st[lin] + (T) 1;

    if (nonlin::is_linear (std::forward<Ts> (nonlin_args)...)) {
      // linear, fast-path
      st[lin] = vec_set<V> (1);
      return u;
    }
    V sig_in = k * (r1.G * u + r1.S);
    st[lin]  = nonlin::tick_div_in (sig_in, std::forward<Ts> (nonlin_args)...);
    if constexpr (n_iters == 1) {
      return u;
    }
    // 2 st round (compensation)
    V u2 = -r2.G * r1.S * k * st[lin] - r2.S + in;
    u2 /= r1.G * r2.G * k * st[lin] + (T) 1;
    u       = (u + u2) * (T) 0.5;
    sig_in  = k * (r1.G * u + r1.S);
    st[lin] = nonlin::tick_div_in (sig_in, std::forward<Ts> (nonlin_args)...);
    return u;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

}} // namespace artv::zdf
