#pragma once

#include <array>
#include <cmath>

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/composite/cascade.hpp"
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
// Variable order Butterworth as a cascade of onepole TDF2's and Andy SVF's
template <class Filter_mode_tag>
class butterworth_any_order {
public:
  using mode_tag     = Filter_mode_tag;
  using cascade_type = filter_cascade_any_order<mode_tag>;

  static constexpr auto max_order = cascade_type::max_order;

  // clang-format off
  static_assert (
    std::is_same_v<mode_tag, highpass_tag> ||
    std::is_same_v<mode_tag, lowpass_tag>,
    "");
  // clang-format on
  //----------------------------------------------------------------------------
  static constexpr uint n_coeffs_for_order (uint order)
  {
    return cascade_type::n_coeffs_for_order (order);
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_states_for_order (uint order)
  {
    return cascade_type::n_states_for_order (order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs interleaved
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T                 = vec_value_type_t<V>;
    constexpr auto traits   = vec_traits<V>();
    constexpr auto max_size = butterworth_2p_cascade_q_list::size (max_order);

    std::array<double, max_size> q_list_mem {};
    auto q_list = butterworth_2p_cascade_q_list::get (q_list_mem, order);
    if constexpr (std::is_same_v<T, double>) {
      cascade_type::reset_coeffs (co, freq, q_list, sr, order);
    }
    else {
      std::array<T, max_size> q_list_cast_mem;
      for (uint i = 0; i < q_list.size(); ++i) {
        q_list_cast_mem[i] = (T) q_list[i];
      }
      auto q_list_cast = make_crange (q_list_cast_mem.data(), q_list.size());
      cascade_type::reset_coeffs (co, freq, q_list_cast, sr, order);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>       dst,
    crange<vec_value_type_t<const V>> src)
  {
    cascade_type::fix_unsmoothable_coeffs (dst, src);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint order)
  {
    cascade_type::reset_states (st, order);
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
    return cascade_type::tick (co, st, in, order, t);
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
    return cascade_type::tick (co, st, in, order);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
using butterworth_lowpass_any_order  = butterworth_any_order<lowpass_tag>;
using butterworth_highpass_any_order = butterworth_any_order<highpass_tag>;
using butterworth_allpass_any_order  = butterworth_any_order<allpass_tag>;
//------------------------------------------------------------------------------
// Fixed order Butterworth as a cascade of onepole TDF2's and Andy SVF's
template <uint N, class Filter_mode_tag>
class butterworth {
public:
  using mode_tag         = Filter_mode_tag;
  using butterworth_type = butterworth_any_order<mode_tag>;

  static constexpr uint order    = N;
  static constexpr uint n_states = butterworth_type::n_states_for_order (order);
  static constexpr uint n_coeffs = butterworth_type::n_coeffs_for_order (order);
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr,
    lowpass_tag)
  {
    butterworth_type::reset_coeffs (co, freq, sr, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st)
  {
    butterworth_type::reset_states (st, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>       dst,
    crange<vec_value_type_t<const V>> src)
  {
    butterworth_type::fix_unsmoothable_coeffs (dst, src);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in,
    single_coeff_set_tag              t)
  {
    return butterworth_type::tick (co, st, in, order, t);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    crange<vec_value_type_t<V>>       st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    return butterworth_type::tick (co, st, in, order);
  }
};
//------------------------------------------------------------------------------
template <uint N>
using butterworth_lowpass = butterworth<N, lowpass_tag>;
template <uint N>
using butterworth_highpass = butterworth<N, highpass_tag>;
template <uint N>
using butterworth_allpass = butterworth<N, allpass_tag>;
//------------------------------------------------------------------------------
} // namespace artv
