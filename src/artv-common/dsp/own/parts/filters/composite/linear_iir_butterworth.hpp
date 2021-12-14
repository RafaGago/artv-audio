#pragma once

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros_reversed.hpp"

// TODO: it didn't behave as expected with the trapezoidal one and two pole
// filters (Andy's SVF and Vadim's onepole).
namespace artv {
//------------------------------------------------------------------------------
struct linear_iir_butterworth_order_2_lowpass {
  //----------------------------------------------------------------------------
  static constexpr uint n_complex = 2;
  static constexpr uint order     = 2;
  //----------------------------------------------------------------------------
  enum coeffs {
    gain = t_rev_rpole_rzero::n_coeffs + rpole_rzero::n_coeffs,
    n_coeffs
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return rpole_rzero::n_states + t_rev_rpole_rzero::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_latency (uint n_stages) { return 1 << n_stages; }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           freq,
    vec_value_type_t<V>         sr)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    vec_complex<V> poles, zeros;
    butterworth_lp_complex::poles (make_crange (poles), freq, sr, 1);
    butterworth_lp_complex::zeros (make_crange (zeros), freq, sr, 1);
    V gain_v = butterworth_lp_complex::gain (make_crange (poles), 1);

    auto gain_ptr = &co[gain * traits.size]; // gain is last
    vec_store (gain_ptr, gain_v);

    t_rev_rpole_rzero::reset_coeffs (co, vec_real (poles), vec_real (zeros));
    co = co.shrink_head (t_rev_rpole_rzero::n_coeffs * traits.size);

    rpole_rzero::reset_coeffs (co, vec_real (poles), vec_real (zeros));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>,
    uint stages)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 v,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (n_coeffs * traits.size));
    assert (st.size() >= (get_n_states (n_stages) * traits.size));
    assert (n_stages);

    V gain_v = vec_load<V> (&co[gain * traits.size]);
    V out    = v;

    out = t_rev_rpole_rzero::tick (co, st, out, n_stages, sample_idx);
    out *= gain_v;
    co = co.shrink_head (t_rev_rpole_rzero::n_coeffs * traits.size);
    st = st.shrink_head (
      t_rev_rpole_rzero::get_n_states (n_stages) * traits.size);

    out = rpole_rzero::tick (co, st, out);
    out *= gain_v;

    return out;
  }
  //----------------------------------------------------------------------------
};

// just multiples of 4.
struct linear_iir_butterworth_2pole_cascade_lowpass {
  //----------------------------------------------------------------------------
  static constexpr uint n_complex = 2;
  static constexpr uint min_order = 4;
  static constexpr uint max_order = 8;
  enum coeffs {
    gain,
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_coeffs (uint order)
  {
    uint n_2poles = order / 4;
    uint fwd      = fwd_2pole::n_coeffs * n_2poles;
    uint rev      = rev_2pole::n_coeffs * n_2poles;
    return fwd + rev + 1; // + 1: gain
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint order, uint n_stages)
  {
    uint n_2poles = order / 4;
    uint fwd      = fwd_2pole::n_states * n_2poles;
    uint rev      = rev_2pole::get_n_states (n_stages) * n_2poles;
    return fwd + rev + 1; // + 1: gain
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_latency (uint order, uint n_stages)
  {
    uint n_2poles = order / 4;
    return ((1 << n_stages) + 1) * n_2poles;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert ((order % 4) == 0);

    std::array<vec_complex<V>, 1>             zeros;
    std::array<vec_complex<V>, max_order / 2> poles;

    butterworth_lp_complex::zeros (make_crange (zeros), freq, sr, 1);
    butterworth_lp_complex::poles (make_crange (poles), freq, sr, order / 2);
    V gain_v = butterworth_lp_complex::gain (make_crange (poles), order / 2);

    vec_store (&co[gain * traits.size], gain_v);
    co = co.shrink_head (traits.size); // advance the gain coefficient

    for (uint i = 0; i < (order / 4); ++i) {
      rev_2pole::reset_coeffs (co, poles[i], vec_real (zeros[0]));
      co = co.shrink_head (rev_2pole::n_coeffs * traits.size);
    }
    for (uint i = 0; i < (order / 4); ++i) {
      fwd_2pole::reset_coeffs (co, poles[i], vec_real (zeros[0]));
      co = co.shrink_head (fwd_2pole::n_coeffs * traits.size);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>       dst,
    crange<vec_value_type_t<const V>> src,
    uint                              order,
    uint                              stages)
  {
    // TODO: it is actually known to be a no-op. Left for when required.
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (
    crange<vec_value_type_t<V>> st,
    uint                        order,
    uint                        stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (order, stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 v,
    uint                              order,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (get_n_coeffs (order) * traits.size));
    assert (st.size() >= (get_n_states (order, n_stages) * traits.size));
    assert (n_stages);
    assert ((order % 4) == 0);

    V out    = v;
    V gain_v = vec_load<V> (&co[gain * traits.size]);
    co       = co.shrink_head (traits.size);

    for (uint i = 0; i < (order / 4); ++i) {
      out = rev_2pole::tick (co, st, out, n_stages, sample_idx);
      co  = co.shrink_head (rev_2pole::n_coeffs * traits.size);
      st  = st.shrink_head (rev_2pole::get_n_states (n_stages) * traits.size);
    }
    out *= gain_v;

    for (uint i = 0; i < (order / 4); ++i) {
      out = fwd_2pole::tick (co, st, out);
      co  = co.shrink_head (fwd_2pole::n_coeffs * traits.size);
      st  = st.shrink_head (fwd_2pole::n_states * traits.size);
    }
    out *= gain_v;
    return out;
  }
  //----------------------------------------------------------------------------
private:
  using rev_2pole = t_rev_ccpole_pair_rzero_pair;
  using fwd_2pole = ccpole_pair_rzero_pair;
};
//------------------------------------------------------------------------------
struct linear_iir_butterworth_lowpass_any_order {
  //----------------------------------------------------------------------------
  static constexpr uint max_order
    = linear_iir_butterworth_2pole_cascade_lowpass::max_order;
  //----------------------------------------------------------------------------
  static constexpr uint get_n_coeffs (uint order)
  {
    return (order == 2) ? two::n_coeffs : any::get_n_coeffs (order);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint order, uint n_stages)
  {
    return (order == 2) ? two::get_n_states (n_stages)
                        : any::get_n_states (order, n_stages);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_latency (uint order, uint n_stages)
  {
    return (order == 2) ? two::get_latency (n_stages)
                        : any::get_latency (order, n_stages);
  }
  //----------------------------------------------------------------------------
  static uint get_n_stages (float freq, float samplerate, float snr_db)
  {
    vec_complex<double_x2> pole;
    // worst case, 1 pole.
    butterworth_lp_complex::poles (
      make_crange (pole), vec_set<double_x2> (freq), (double) samplerate, 1);
    return get_reversed_pole_n_stages (
      {vec_real (pole)[0], vec_imag (pole)[0]}, snr_db);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co, // coeffs (interleaved, SIMD aligned)
    V                           freq,
    vec_value_type_t<V>         sr,
    uint                        order)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      two::reset_coeffs (co, freq, sr);
    }
    else {
      any::reset_coeffs (co, freq, sr, order);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>       dst,
    crange<vec_value_type_t<const V>> src,
    uint                              order,
    uint                              n_stages)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      two::fix_unsmoothable_coeffs (dst, src, n_stages);
    }
    else {
      any::fix_unsmoothable_coeffs (dst, src, order, n_stages);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static auto tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 v,
    uint                              order,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      return two::tick (co, st, v, n_stages, sample_idx);
    }
    else {
      return any::tick (co, st, v, order, n_stages, sample_idx);
    }
  }
  //----------------------------------------------------------------------------
private:
  using two = linear_iir_butterworth_order_2_lowpass;
  using any = linear_iir_butterworth_2pole_cascade_lowpass;
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
