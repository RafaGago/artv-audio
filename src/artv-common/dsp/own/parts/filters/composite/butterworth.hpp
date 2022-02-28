#pragma once

#include <array>
#include <cmath>

#include <gcem.hpp>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

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
// Variable order Butterworth as a cascade of onepole TPT's and Andy SVF's
template <class Filter_mode_tag, uint MaxOrder = 16>
class butterworth_any_order {
public:
  static constexpr auto max_order = MaxOrder;

  using mode_tag     = Filter_mode_tag;
  using cascade_type = filter_cascade_any_order<mode_tag, max_order>;

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
  static constexpr uint n_coeffs     = cascade_type::n_coeffs;
  static constexpr uint n_coeffs_int = cascade_type::n_coeffs_int;
  static constexpr uint n_states     = cascade_type::n_states;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co, // coeffs interleaved
    V                   freq,
    vec_value_type_t<V> sr,
    uint                order)
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
  static void reset_states (crange<V> st, uint order)
  {
    cascade_type::reset_states (st, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in,
    uint                              order)
  {
    return cascade_type::tick (co, st, in, order);
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
    return cascade_type::tick (co, st, in, order);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
using butterworth_lowpass_any_order  = butterworth_any_order<lowpass_tag>;
using butterworth_highpass_any_order = butterworth_any_order<highpass_tag>;
using butterworth_allpass_any_order  = butterworth_any_order<allpass_tag>;
//------------------------------------------------------------------------------
// Fixed order Butterworth as a cascade of onepole TPT's and Andy SVF's
template <uint N, class Filter_mode_tag>
class butterworth {
public:
  using mode_tag         = Filter_mode_tag;
  using butterworth_type = butterworth_any_order<mode_tag, N>;

  static constexpr uint order        = N;
  static constexpr uint n_states     = butterworth_type::n_states;
  static constexpr uint n_coeffs_int = butterworth_type::n_coeffs_int;
  static constexpr uint n_coeffs     = butterworth_type::n_coeffs;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>           co, // coeffs (interleaved, SIMD aligned)
    V                   freq,
    vec_value_type_t<V> sr,
    lowpass_tag)
  {
    butterworth_type::reset_coeffs (co, freq, sr, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    butterworth_type::reset_states (st, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co, // coeffs (single set)
    crange<V>                         st, // states (interleaved, SIMD aligned)
    V                                 in)
  {
    return butterworth_type::tick (co, st, in, order);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co, // coeffs (interleaved, SIMD aligned)
    crange<V>       st, // states (interleaved, SIMD aligned)
    V               in)
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
struct butterworth_lp_complex {
  // returns first the poles:
  //  - central (if odd order)
  //  - positive imaginary
  //  - negative imaginary / conjugates
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void poles (
    crange<vec_complex<V>> poles,
    V                      freq,
    vec_value_type_t<V>    srate,
    uint                   order)
  {
    /*
     This function places the poles of a normalized filter on the S-plane.

                 (s - z1) (s - z2) ... ( s - zn)
      H(s) = K * -------------------------------
                 (s - p1) (s - p2) ... ( s - pn)


     Then does the bilinear transform to the z plane.

                (z - z1) (z - z2) ... ( z - zn)
     H(z) = K * -------------------------------
                (z - p1) (z - p2) ... ( z - pn)

    It just returns the poles. The gain and the zeros are to be obtained
    separately.
    */

    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (poles.size() >= order);

    V w         = freq / srate;
    V w_prewarp = vec_tan ((T) M_PI * w);

    T poles_angle = (T) M_PI / (T) order;
    T degrees;

    uint idx = 0;
    if (order & 1) {
      degrees = M_PI_2 + poles_angle;
    }
    else {
      degrees = M_PI_2 + (poles_angle * (T) 0.5);
    }
    // positive Im
    for (uint i = 0; i < (order / 2); ++i) {
      auto pole = vec_polar<V> (1., degrees);
      pole *= w_prewarp; // freq transform
      poles[idx] = (1. + pole) / (1. - pole); // BLT LP
      degrees += poles_angle; // float addition will add error, but not a lot.
      ++idx;
    }
    // central
    if (order & 1) {
      auto pole = vec_complex<V> {-1.};
      pole *= w_prewarp; // freq transform
      poles[idx] = (1. + pole) / (1. - pole); // BLT LP
      ++idx;
    }
    // conjugates
    for (uint i = 0; i < (order / 2); ++i) {
      poles[idx] = vec_conj (poles[idx - (order / 2)]);
      ++idx;
    }
  }
  //----------------------------------------------------------------------------
  // just a formality. All zeroes are real at -1.
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void zeros (
    crange<vec_complex<V>> zeros,
    V,
    vec_value_type_t<V>,
    uint order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (zeros.size() >= order);

    for (uint i = 0; i < order; ++i) {
      zeros[i] = vec_complex<V> {(T) -1., (T) 0.};
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V gain (const crange<vec_complex<V>> poles, uint order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    /*
    This gets poles on the Z-plane (the zeros are implied):

               (z - z1) (z - z2) ... ( z - zn)
    H(z) = K * -------------------------------
               (z - p1) (z - p2) ... ( z - pn)


    Substituting z by 1 (evaluate at DC), the gain for 1 (normalized) and the
    zeros by -1, that would leave:

                      2^(n zeros)
    1 = K * -------------------------------
             (1 - p1) (1 - p2) ... ( 1 - pn)


        (1 - p1) (1 - p2) ... ( 1 - pn)
    K = -------------------------------
                 2^(n zeros)
    */

    assert (poles.size() >= order);

    auto den = (T) (1u << order);
    auto num = vec_set<V> ((T) 1.);

    // handling conjugates.
    uint i = 0;
    for (; i < (order / 2); ++i) {
      num *= vec_real (vec_conj_mul (poles[i]));
    }
    // handling real pole if present.
    if (order & 1) {
      num *= (T) 1 - vec_real (poles[i]);
    }

    V ret = num / den;
    return ret;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
