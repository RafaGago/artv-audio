// Linear-phase IIR butterworth based on a time reversed IIR
#pragma once

#include <cmath>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_complex.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/filters/onepole.hpp"

#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros_reversed.hpp"

// TODO: it didn't behave as expected with the trapezoidal one and two pole
// filters (Andy's SVF and Vadim's onepole).
namespace artv {
//------------------------------------------------------------------------------
class linear_iir_butterworth_order_2_lowpass {
private:
  using rev_1pole = t_rev_rpole_rzero;
  using fwd_1pole = onepole_lowpass;

public:
  //----------------------------------------------------------------------------
  static constexpr uint n_complex = 2;
  static constexpr uint order     = 2;
  //----------------------------------------------------------------------------
  enum coeffs { gain = fwd_1pole::n_coeffs + rev_1pole::n_coeffs, n_coeffs };
  enum coeffs_int {
    n_coeffs_int = fwd_1pole::n_coeffs_int + rev_1pole::n_coeffs_int,
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return fwd_1pole::n_states + rev_1pole::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_latency (uint n_stages) { return 1 << n_stages; }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V freq, vec_value_type_t<V> t_spl)
  {
    vec_complex<V> poles, zeros;
    butterworth_lp_complex::poles (xspan {&poles, 1}, freq, t_spl, 1);
    butterworth_lp_complex::zeros (xspan {&zeros, 1}, freq, t_spl, 1);
    // gain is last, storing before moving "co"
    co[gain] = butterworth_lp_complex::gain (xspan {&poles, 1}, 1);

    rev_1pole::reset_coeffs (co, vec_real (poles), vec_real (zeros));
    co.cut_head (rev_1pole::n_coeffs);

    fwd_1pole::reset_coeffs (co, freq, t_spl);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st, uint stages)
  {
    uint numstates = get_n_states (stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static auto tick (
    xspan<V const> co,
    xspan<V>       st,
    V              v,
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    assert (co.size() >= (n_coeffs));
    assert (st.size() >= (get_n_states (n_stages)));
    assert (n_stages);

    V gain_v = co[gain];
    V out    = v;

    out = rev_1pole::tick (co, st, out, n_stages, sample_idx);
    out *= gain_v;
    co.cut_head (rev_1pole::n_coeffs);
    st.cut_head (rev_1pole::get_n_states (n_stages));

    out = fwd_1pole::tick (co, st, out);
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void tick (
    xspan<V const> co,
    xspan<V>       st,
    xspan<V>       io, // ins on call, outs when returning
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages);

    V gain_v = co[gain];

    rev_1pole::tick (co, st, io, n_stages, sample_idx);
    co.cut_head (rev_1pole::n_coeffs);
    st.cut_head (rev_1pole::get_n_states (n_stages));

    for (uint i = 0; i < io.size(); ++i) {
      io[i] *= gain_v; // gain from previous stage
      io[i] = fwd_1pole::tick (co, st, io[i]);
    }
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
    uint fwd      = fwd_2pole::n_coeffs_for_order (order / 2);
    uint rev      = rev_2pole::n_coeffs * n_2poles;
    return fwd + rev + 1; // + 1: gain
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint order, uint n_stages)
  {
    uint n_2poles = order / 4;
    uint fwd      = fwd_2pole::n_states_for_order (order / 2);
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
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co,
    V                   freq,
    vec_value_type_t<V> t_spl,
    uint                order)
  {
    assert ((order % 4) == 0);

    std::array<vec_complex<V>, 1>             zeros;
    std::array<vec_complex<V>, max_order / 2> poles;

    butterworth_lp_complex::zeros (xspan {zeros}, freq, t_spl, 1);
    butterworth_lp_complex::poles (xspan {poles}, freq, t_spl, order / 2);

    co[gain] = butterworth_lp_complex::gain (xspan {poles}, order / 2);
    co.cut_head (1); // advance the gain coefficient

    for (uint i = 0; i < (order / 4); ++i) {
      rev_2pole::reset_coeffs (co, poles[i], vec_real (zeros[0]));
      co.cut_head (rev_2pole::n_coeffs);
    }
    fwd_2pole::reset_coeffs (co, freq, t_spl, order / 2);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st, uint order, uint stages)
  {
    uint numstates = get_n_states (order, stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static auto tick (
    xspan<V const> co,
    xspan<V>       st,
    V              v,
    uint           order,
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    assert (co.size() >= get_n_coeffs (order));
    assert (st.size() >= get_n_states (order, n_stages));
    assert (n_stages);
    assert ((order % 4) == 0);

    V out    = v;
    V gain_v = co[gain];
    co.cut_head (1);

    for (uint i = 0; i < (order / 4); ++i) {
      out = rev_2pole::tick (co, st, out, n_stages, sample_idx);
      co.cut_head (rev_2pole::n_coeffs);
      st.cut_head (rev_2pole::get_n_states (n_stages));
    }
    out *= gain_v;
    return fwd_2pole::tick (co, st, out, order / 2);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void tick (
    xspan<V const> co,
    xspan<V>       st,
    xspan<V>       io, // ins on call, outs when returning
    uint           order,
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= get_n_coeffs (order));
    assert (st.size() >= get_n_states (order, n_stages));
    assert (n_stages);
    assert ((order % 4) == 0);

    V gain_v = co[gain];
    co.cut_head (1);

    for (uint i = 0; i < (order / 4); ++i) {
      rev_2pole::tick (co, st, io, n_stages, sample_idx);
      co.cut_head (rev_2pole::n_coeffs);
      st.cut_head (rev_2pole::get_n_states (n_stages));
    }
    for (uint i = 0; i < io.size(); ++i) {
      io[i] *= gain_v;
      io[i] = fwd_2pole::tick (co, st, io[i], order / 2);
    }
  }
  //----------------------------------------------------------------------------
private:
  using rev_2pole = t_rev_ccpole_pair_rzero_eq_pair;
  using fwd_2pole = butterworth_lowpass_any_order;
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
  static uint get_n_stages (float freq, float t_spl, float snr_db)
  {
    vec_complex<f64_x2> pole;
    // worst case, 1 pole.
    butterworth_lp_complex::poles (
      xspan {&pole, 1}, vec_set<f64_x2> (freq), (double) t_spl, 1);
    return get_reversed_pole_n_stages (
      {vec_real (pole)[0], vec_imag (pole)[0]}, snr_db);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (
    xspan<V>            co, // coeffs (interleaved, SIMD aligned)
    V                   freq,
    vec_value_type_t<V> t_spl,
    uint                order)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      two::reset_coeffs<V> (co, freq, t_spl);
    }
    else {
      any::reset_coeffs<V> (co, freq, t_spl, order);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static auto tick (
    xspan<V const> co,
    xspan<V>       st,
    V              v,
    uint           order,
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      return two::tick<V> (co, st, v, n_stages, sample_idx);
    }
    else {
      return any::tick<V> (co, st, v, order, n_stages, sample_idx);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void tick (
    xspan<V const> co,
    xspan<V>       st,
    xspan<V>       io,
    uint           order,
    uint           n_stages,
    uint           sample_idx) // sample counter (external)
  {
    assert ((order % 2) == 0);
    assert (order <= max_order);
    if (order == 2) {
      two::tick<V> (co, st, io, n_stages, sample_idx);
    }
    else {
      any::tick<V> (co, st, io, order, n_stages, sample_idx);
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
