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

// A variable order cascade of Cytomic SVF's and TDF2 onepole (odd orders).
template <class Filter_mode_tag>
class filter_cascade_any_order {
public:
  //----------------------------------------------------------------------------
  using mode_tag                  = Filter_mode_tag;
  using svf_type                  = andy::svf_multimode<mode_tag>;
  static constexpr uint max_order = 16;
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * svf_type::n_coeffs + n_onepole * onepole::n_coeffs;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_states_for_order (uint order)
  {
    uint n_onepole = (order % 2);
    uint n_twopole = order - n_onepole;
    return n_twopole * svf_type::n_states + n_onepole * onepole::n_states;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs interleaved
    V                           freq,
    crange<vec_value_type_t<V>> q_list, // one Q for each cascaded SVF
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();
    uint           n_svfs = (order / 2);

    assert (co.size() >= (n_coeffs_for_order (order) * traits.size));
    assert (order <= max_order);
    assert (q_list.size() >= n_svfs);

    if (order & 1) {
      onepole::reset_coeffs (co, freq, sr, mode_tag {});
      co = co.shrink_head (onepole::n_coeffs * traits.size);
    }
    for (uint i = 0; i < n_svfs; ++i) {
      svf_type::reset_coeffs (co, freq, vec_set<V> (q_list[i]), sr);
      co = co.shrink_head (svf_type::n_coeffs * traits.size);
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
    single_coeff_set_tag              t)
  {
    return tick_impl<single_coeff_set_tag> (co, st, in, order);
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
    return tick_impl (co, st, in, order);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class Tag = void, class V>
  static V tick_impl (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    using T                        = vec_value_type_t<V>;
    constexpr auto traits          = vec_traits<V>();
    constexpr bool is_single_coeff = !std::is_same_v<Tag, void>;
    constexpr auto coeff_vec_size  = is_single_coeff ? 1 : traits.size;

    assert (co.size() >= (coeff_vec_size * n_coeffs_for_order (order)));
    assert (st.size() >= (traits.size * n_states_for_order (order)));

    auto out = in;

    if (order & 1) {
      if constexpr (!is_single_coeff) {
        out = onepole::tick (co, st, out);
      }
      else {
        out = onepole::tick (co, st, out, Tag {});
      }
      co = co.shrink_head (onepole::n_coeffs * coeff_vec_size);
      st = st.shrink_head (onepole::n_states * traits.size);
    }
    for (uint i = 0; i < (order / 2); ++i) {
      if constexpr (!is_single_coeff) {
        out = svf_type::tick (co, st, out);
      }
      else {
        out = svf_type::tick (co, st, out, Tag {});
      }
      co = co.shrink_head (svf_type::n_coeffs * coeff_vec_size);
      st = st.shrink_head (svf_type::n_states * traits.size);
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
  static void fix_unsmoothable_coeffs (crange<double>, crange<const double>) {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    crange<vec_value_type_t<V>> q_list, // one Q for each cascaded SVF
    vec_value_type_t<V>         sr)
  {
    cascade_type::reset_coeffs (co, freq, q_list, sr, order);
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
    return cascade_type::tick (co, st, in, order, t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
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
