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

// A variable order cascade of Cytomic SVF's 2-poles and TDF2 1-pole (odd
// orders).
template <class Filter_mode_tag>
class filter_cascade_any_order {
public:
  //----------------------------------------------------------------------------
  using mode_tag                  = Filter_mode_tag;
  using svf_type                  = andy::svf_multimode<mode_tag>;
  using onepole_type              = onepole<mode_tag>;
  static constexpr uint max_order = 16;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs_for_order (uint order)
  {
    uint n_onepole = order % 2;
    uint n_twopole = order / 2;
    return n_twopole * svf_type::n_coeffs + n_onepole * onepole_type::n_coeffs;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_states_for_order (uint order)
  {
    uint n_onepole = order % 2;
    uint n_twopole = order / 2;
    return n_twopole * svf_type::n_states + n_onepole * onepole_type::n_states;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>                   co, // coeffs interleaved
    V                           freq,
    crange<vec_value_type_t<V>> q_list, // one Q for each cascaded SVF
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    uint           n_svfs = (order / 2);

    assert (co.size() >= (n_coeffs_for_order (order)));
    assert (order <= max_order);
    assert (q_list.size() >= n_svfs);

    if (order & 1) {
      onepole_type::reset_coeffs (co, freq, sr);
      co = co.shrink_head (onepole_type::n_coeffs);
    }
    for (uint i = 0; i < n_svfs; ++i) {
      svf_type::reset_coeffs (co, freq, vec_set<V> (q_list[i]), sr);
      co = co.shrink_head (svf_type::n_coeffs);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint order)
  {
    uint numstates = n_states_for_order (order);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    return tick_impl (co, st, in, order);
  }
  //----------------------------------------------------------------------------
  // N sets of coeffs, N outs calculated at once.
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in,
    uint            order)
  {
    return tick_impl (co, st, in, order);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class Co, class V>
  static V tick_impl (
    crange<const Co> co, // coeffs (interleaved or single set)
    crange<V>        st, // states (interleaved, SIMD aligned)
    V                in,
    uint             order)
  {
    assert (co.size() >= n_coeffs_for_order (order));
    assert (st.size() >= n_states_for_order (order));

    auto out = in;

    if (order & 1) {
      out = onepole_type::tick (co, st, out);
      co  = co.shrink_head (onepole_type::n_coeffs);
      st  = st.shrink_head (onepole_type::n_states);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      out = svf_type::tick (co, st, out);
      co  = co.shrink_head (svf_type::n_coeffs);
      st  = st.shrink_head (svf_type::n_states);
    }
    return out;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
using lowpass_cascade_any_order  = filter_cascade_any_order<lowpass_tag>;
using highpass_cascade_any_order = filter_cascade_any_order<highpass_tag>;
using allpass_cascade_any_order  = filter_cascade_any_order<allpass_tag>;
//------------------------------------------------------------------------------
// A fixed order cascade of Cytomic SVF's and TDF2 onepole (odd orders).
template <uint N, class Filter_mode_tag>
class filter_cascade {
public:
  using mode_tag                 = Filter_mode_tag;
  using cascade_type             = filter_cascade_any_order<mode_tag>;
  static constexpr uint order    = N;
  static constexpr uint n_states = cascade_type::n_states_for_order (order);
  static constexpr uint n_coeffs = cascade_type::n_coeffs_for_order (order);

  static_assert (
    order > 0 && order <= cascade_type::max_order
    && "Unsupported filter order.");

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>                   co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    crange<vec_value_type_t<V>> q_list, // one Q for each cascaded SVF
    vec_value_type_t<V>         sr)
  {
    cascade_type::reset_coeffs (co, freq, q_list, sr, order);
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
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    return cascade_type::tick (co, st, in, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
  {
    return cascade_type::tick (co, st, in, order);
  }
};
//------------------------------------------------------------------------------
template <uint N>
using lowpass_cascade = filter_cascade<N, lowpass_tag>;
template <uint N>
using highpass_cascade = filter_cascade<N, highpass_tag>;
template <uint N>
using allpass_cascade = filter_cascade<N, allpass_tag>;
//------------------------------------------------------------------------------
} // namespace artv
