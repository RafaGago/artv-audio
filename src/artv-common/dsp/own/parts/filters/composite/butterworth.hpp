#pragma once

#include <array>
#include <cmath>

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/traits.hpp"

namespace artv {

struct butterworth_2p_cascade_q_list {
  //----------------------------------------------------------------------------
  static constexpr uint size (uint order) { return order / 2; };
  //----------------------------------------------------------------------------
  template <uint Order>
  static constexpr uint size_v = size (Order);
  //----------------------------------------------------------------------------
  template <uint Max_order>
  using array = std::array<double, size_v<Max_order>>;
  //----------------------------------------------------------------------------
  static crange<double> get (crange<double> res_mem, uint order)
  {
    // odd orders have a single pole at the front.
    uint sz = size (order);
    assert (res_mem.size() >= sz);

    long double step  = M_PI / (double) order;
    long double angle = order & 1 ? step : step * 0.5;

    for (uint i = 0; i < sz; ++i, angle += step) {
      res_mem[i] = 1. / (2. * cos (angle));
    }
    return {res_mem.data(), sz};
  }
  //----------------------------------------------------------------------------
  template <uint Order>
  static array<Order> get()
  {
    array<Order> ret;
    get (ret, Order);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Order>
  static constexpr array<Order> cget()
  {
    array<Order> ret {};

    long double step  = M_PI / (long double) Order;
    long double angle = Order & 1 ? step : step * 0.5;

    for (uint i = 0; i < ret.size(); ++i, angle += step) {
      ret[i] = 1. / (2. * gcem::cos (angle));
    }
    return ret;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// Butterworth as a cascade of Andy SVF's
class butterworth_any_order {
public:
  static constexpr uint max_order = 16;
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
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> co, // coeffs interleaved
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order,
    bool                        is_lowpass)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs_for_order (order) * traits.size));
    assert (order <= max_order);

    constexpr auto max_size = butterworth_2p_cascade_q_list::size (max_order);
    std::array<double, max_size> q_list_mem;
    auto q_list = butterworth_2p_cascade_q_list::get (q_list_mem, order);

    if (order & 1) {
      if (is_lowpass) {
        onepole::init (co, freq, sr, lowpass_tag {});
      }
      else {
        onepole::init (co, freq, sr, highpass_tag {});
      }
      co = co.shrink_head (onepole::n_coeffs * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      if (is_lowpass) {
        andy::svf::init (co, freq, vec_set<V> (q_list[i]), sr, lowpass_tag {});
      }
      else {
        andy::svf::init (co, freq, vec_set<V> (q_list[i]), sr, highpass_tag {});
      }
      co = co.shrink_head (andy::svf::n_coeffs * traits.size);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states_for_order (order);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order,
    single_coeff_set_tag)
  {

    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs_for_order (order));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    if (order & 1) {
      out = onepole::tick (co, st, out, single_coeff_set_tag {});
      co  = co.shrink_head (onepole::n_coeffs);
      st  = st.shrink_head (onepole::n_states * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out = andy::svf::tick (co, st, out, single_coeff_set_tag {});
      co  = co.shrink_head (andy::svf::n_coeffs);
      st  = st.shrink_head (andy::svf::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    if (order & 1) {
      out = onepole::tick (co, st, out);
      co  = co.shrink_head (onepole::n_coeffs * traits.size);
      st  = st.shrink_head (onepole::n_states * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out = andy::svf::tick (co, st, out);
      co  = co.shrink_head (andy::svf::n_coeffs * traits.size);
      st  = st.shrink_head (andy::svf::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
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
  static void fix_unsmoothable_coeffs (crange<double>, crange<const double>) {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    butterworth_any_order::init (co, freq, sr, order, true);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr,
    highpass_tag)
  {
    butterworth_any_order::init (co, freq, sr, order, false);
  }
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
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag              t)
  {
    return butterworth_any_order::tick (co, st, in, order, t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    return butterworth_any_order::tick (co, st, in, order);
  }
};
//------------------------------------------------------------------------------
// The lowpass-only variant uses 3 coefficients less. The tick function saves
// 3 muls and two adds for each cascaded SVF.
class butterworth_lowpass_any_order {
public:
  static constexpr uint max_order = butterworth_any_order::max_order;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * andy::svf_lowpass::n_coeffs
      + n_onepole * onepole::n_coeffs;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_states_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * andy::svf_lowpass::n_states
      + n_onepole * onepole::n_states;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void init (
    crange<vec_value_type_t<V>> co, // coeffs interleaved
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs_for_order (order) * traits.size));
    assert (order <= max_order);

    constexpr auto max_size = butterworth_2p_cascade_q_list::size (max_order);
    std::array<double, max_size> q_list_mem;
    auto q_list = butterworth_2p_cascade_q_list::get (q_list_mem, order);

    if (order & 1) {
      onepole::init (co, freq, sr, lowpass_tag {});
      co = co.shrink_head (onepole::n_coeffs * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      andy::svf_lowpass::init (co, freq, vec_set<V> (q_list[i]), sr);
      co = co.shrink_head (andy::svf_lowpass::n_coeffs * traits.size);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * n_states_for_order (order);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    if (order & 1) {
      out = onepole::tick (co, st, out);
      co  = co.shrink_head (onepole::n_coeffs * traits.size);
      st  = st.shrink_head (onepole::n_states * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out = andy::svf_lowpass::tick (co, st, out);
      co  = co.shrink_head (andy::svf_lowpass::n_coeffs * traits.size);
      st  = st.shrink_head (andy::svf_lowpass::n_states * traits.size);
    }
    return out;
  }
#if 0 // I was ready to introduce these premature optimizations, delaying...
  //----------------------------------------------------------------------------
  // for block processing optimization, so one pass of the onepole is done on
  // all the samples, then a pass of the cascades.
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick_onepole (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    if (order & 1) {
      out = onepole::tick (co, st, out);
      co  = co.shrink_head (onepole::n_coeffs * traits.size);
      st  = st.shrink_head (onepole::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  // for block processing optimization, so one pass of the onepole is done on
  // all the samples, then a pass of the cascades.
  // N sets of coeffs, N outs calculated at once.
  template <class V, :enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick_two_pole_cascade (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    co = co.shrink_head (onepole::n_coeffs * traits.size * (order & 1));
    st = st.shrink_head (onepole::n_states * traits.size * (order & 1));

    for (uint i = 0; i < (order / 2); ++i) {
      out = andy::svf_lowpass::tick (co, st, out);
      co  = co.shrink_head (andy::svf_lowpass::n_coeffs * traits.size);
      st  = st.shrink_head (andy::svf_lowpass::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
#endif
};

//------------------------------------------------------------------------------
} // namespace artv
