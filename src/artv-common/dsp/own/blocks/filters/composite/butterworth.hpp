#pragma once

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/blocks/filters/andy_svf.hpp"
#include "artv-common/dsp/own/blocks/filters/onepole.hpp"

namespace artv {

class butterworth_any_order {
public:
  static constexpr uint max_order = 12;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * andy::svf::n_coeffs + n_onepole * onepole::n_coeffs;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_states_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * andy::svf::n_states + n_onepole * onepole::n_states;
  }
  //----------------------------------------------------------------------------
  // Making these public serves a dual purpose:
  // -To allow the function of a higher order filter to dispatch at runtiome,
  //  so it can run lower order initializations of the coefficients.
  //- To make it easy for the compiler to reduce template bloat, as they don't
  //  depend on template parameters.
  //----------------------------------------------------------------------------
  static void init (
    crange<double> co, // coeffs
    double         freq,
    double         sr,
    uint           order,
    bool           is_lowpass)
  {
    assert (co.size() >= n_coeffs_for_order (order));

    auto q_list = get_q_list (order);
    if (order & 1) {
      if (is_lowpass) {
        onepole::lowpass (co, freq, sr);
      }
      else {
        onepole::highpass (co, freq, sr);
      }
      co = co.shrink_head (onepole::n_coeffs);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      if (is_lowpass) {
        andy::svf::lowpass (co, freq, q_list[i], sr);
      }
      else {
        andy::svf::highpass (co, freq, q_list[i], sr);
      }
      co = co.shrink_head (andy::svf::n_coeffs);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void init_simd (
    crange<vec_value_type_t<V>> co, // coeffs interleaved
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order,
    bool                        is_lowpass)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs_for_order (order) * traits.size));

    auto q_list = get_q_list (order);
    if (order & 1) {
      if (is_lowpass) {
        onepole::lowpass_simd (co, freq, sr);
      }
      else {
        onepole::highpass_simd (co, freq, sr);
      }
      co = co.shrink_head (onepole::n_coeffs * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      if (is_lowpass) {
        andy::svf::lowpass_simd (co, freq, vec_set<V> (q_list[i]), sr);
      }
      else {
        andy::svf::highpass_simd (co, freq, vec_set<V> (q_list[i]), sr);
      }
      co = co.shrink_head (andy::svf::n_coeffs * traits.size);
    }
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               v0,
    uint                 order)
  {
    assert (co.size() >= n_coeffs_for_order (order));
    assert (st.size() >= n_states_for_order (order));

    if (order & 1) {
      v0 = onepole::tick (co, st, v0);
      co = co.shrink_head (onepole::n_coeffs);
      st = st.shrink_head (onepole::n_states);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      v0 = andy::svf::tick (co, st, v0);
      co = co.shrink_head (andy::svf::n_coeffs);
      st = st.shrink_head (andy::svf::n_states);
    }
    return v0;
  }
  //----------------------------------------------------------------------------
  static double_x2 tick (
    crange<double const>          co, // coeffs (unaligned, uninterleaved)
    std::array<crange<double>, 2> st, // state (unaligned, uninterleaved)
    double_x2                     v0,
    uint                          order)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs_for_order (order));
    assert (st[0].size() >= n_states_for_order (order));
    assert (st[1].size() >= n_states_for_order (order));

    auto out = v0;

    if (order & 1) {
      out   = onepole::tick (co, st, out);
      co    = co.shrink_head (onepole::n_coeffs);
      st[0] = st[0].shrink_head (onepole::n_states);
      st[1] = st[1].shrink_head (onepole::n_states);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out   = andy::svf::tick (co, st, out);
      co    = co.shrink_head (andy::svf::n_coeffs);
      st[0] = st[0].shrink_head (andy::svf::n_states);
      st[1] = st[1].shrink_head (andy::svf::n_states);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 v0,
    uint                              order)
  {
    using T = vec_value_type_t<V>;
    static_assert (std::is_floating_point<T>::value, "");
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = v0;

    if (order & 1) {
      out = onepole::tick_simd (co, st, out);
      co  = co.shrink_head (onepole::n_coeffs * traits.size);
      st  = st.shrink_head (onepole::n_states * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out = andy::svf::tick_simd (co, st, out);
      co  = co.shrink_head (andy::svf::n_coeffs * traits.size);
      st  = st.shrink_head (andy::svf::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static crange<const double> get_q_list (uint order)
  {
    // https://www.earlevel.com/main/2016/09/29/cascading-filters/
    // one pole values (0.5) removed from the list.
    static const double q2[] = {0.70710678};
    static const double q3[] = {1.};
    static const double q4[] = {0.54119610, 1.3065630};
    static const double q5[] = {0.61803399, 1.6180340};
    static const double q6[] = {0.51763809, 0.70710678, 1.9318517};
    static const double q7[] = {0.55495813, 0.80193774, 2.2469796};
    static const double q8[] = {0.50979558, 0.60134489, 0.89997622, 2.5629154};
    static const double q9[] = {0.53208889, 0.65270364, 1.0000000, 2.8793852};
    static const double q10[]
      = {0.50623256, 0.56116312, 0.70710678, 1.1013446, 3.1962266};
    static const double q11[]
      = {0.52110856, 0.59435114, 0.76352112, 1.2036156, 3.5133371};
    static const double q12[]
      = {0.50431448, 0.54119610, 0.63023621, 0.82133982, 1.3065630, 3.8306488};

    switch (order) {
    case 1:
      return {};
    case 2:
      return q2;
    case 3:
      return q3;
    case 4:
      return q4;
    case 5:
      return q5;
    case 6:
      return q6;
    case 7:
      return q7;
    case 8:
      return q8;
    case 9:
      return q9;
    case 10:
      return q10;
    case 11:
      return q11;
    case 12:
      return q12;
    default:
      assert (false);
      return {};
    }
  }
};
//------------------------------------------------------------------------------
template <uint N>
class butterworth {
public:
  static constexpr uint order = N;
  static constexpr uint n_states
    = butterworth_any_order::n_states_for_order (order);
  static constexpr uint n_coeffs
    = butterworth_any_order::n_coeffs_for_order (order);

  static_assert (
    order > 0 && order <= butterworth_any_order::max_order
    && "Unsupported filter order.");
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void lowpass (crange<double> co, double freq, double sr)
  {
    butterworth_any_order::init (co, freq, sr, order, true);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void lowpass_simd (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    butterworth_any_order::init_simd (co, freq, sr, order, true);
  }
  //----------------------------------------------------------------------------
  static void highpass (crange<double> co, double freq, double sr)
  {
    butterworth_any_order::init (co, freq, sr, order, false);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static void highpass_simd (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    butterworth_any_order::init_simd (co, freq, sr, order, false);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               v0)
  {
    return butterworth_any_order::tick (co, st, v0, order);
  }
  //----------------------------------------------------------------------------
  static double_x2 tick (
    crange<double const>          co, // coeffs (unaligned, uninterleaved)
    std::array<crange<double>, 2> st, // coeffs (state, uninterleaved)
    double_x2                     v0s)
  {
    return butterworth_any_order::tick (co, st, v0s, order);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V tick_simd (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 v0s)
  {
    return butterworth_any_order::tick_simd (co, st, v0s, order);
  }
};
//------------------------------------------------------------------------------
} // namespace artv
