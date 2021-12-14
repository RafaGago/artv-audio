#pragma once

// a dummy plugin to experiment with

#include "artv-common/dsp/own/classes/linphase_iir_crossover.hpp"
#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/composite/linear_iir_butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros_reversed.hpp"

#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

#if 1
// variable order butterworth linear phase, wrapper class
//------------------------------------------------------------------------------
class experiments {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 2;
  //----------------------------------------------------------------------------
  struct param_a_tag {};

  void set (param_a_tag, float v) { _a = v; }

  static constexpr auto get_parameter (param_a_tag)
  {
    return float_param ("A", 61., 20000., 600., 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v) { _b = v; }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("B", 61., 20000, 1500., 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_c_tag {};

  void set (param_c_tag, float v) { _c = v; }

  static constexpr auto get_parameter (param_c_tag)
  {
    return float_param ("C", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_d_tag {};

  void set (param_d_tag, float v) { _d = v; }

  static constexpr auto get_parameter (param_d_tag)
  {
    return float_param ("D", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<param_a_tag, param_b_tag, param_c_tag, param_d_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    memset (this, 0, sizeof *this);
    _plugcontext = &pc;

    std::array<float, 3> min_freq {60., 200., 600.};
    _crossv.reset (pc.get_sample_rate(), min_freq, 120.);

    _crossv.set_crossover_point (0, 600., 600., order);
    _crossv.set_crossover_point (1, 2500., 2500., order);
    _crossv.set_crossover_point (2, 5000., 5000., order);
    pc.set_delay_compensation (_crossv.get_latency());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {

      auto vals  = _crossv.tick (double_x2 {ins[0][i], ins[1][i]});
      outs[0][i] = vals[0][0] + /*vals[1][0] + vals[2][0]*/ +vals[3][0];
      outs[1][i] = vals[0][1] + /*vals[1][1] + vals[2][1]*/ +vals[3][0];
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  float           _a, _b, _c, _d;

  static constexpr auto     order = 4;
  linphase_iir_crossover<3> _crossv;
};
//------------------------------------------------------------------------------
#endif
#if 0
// variable order butterworth linear phase
//------------------------------------------------------------------------------
class experiments {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct param_a_tag {};

  void set (param_a_tag, float v) { _a = v; }

  static constexpr auto get_parameter (param_a_tag)
  {
    return float_param ("A", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v) { _b = v; }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("B", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_c_tag {};

  void set (param_c_tag, float v) { _c = v; }

  static constexpr auto get_parameter (param_c_tag)
  {
    return float_param ("C", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_d_tag {};

  void set (param_d_tag, float v) { _d = v; }

  static constexpr auto get_parameter (param_d_tag)
  {
    return float_param ("D", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<param_a_tag, param_b_tag, param_c_tag, param_d_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    memset (this, 0, sizeof *this);
    _plugcontext = &pc;

    lp_type::reset_coeffs (
      make_crange (_coeffs).cast (double {}),
      vec_set<double_x2> (660.),
      pc.get_sample_rate(),
      order);
    // pc.set_delay_compensation (lp_type::get_latency (order, n_stages));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {
      double_x2 in {double_x2 {ins[0][i], ins[1][i]}};

      auto pos    = _sample_idx % lp_type::get_latency (order, n_stages);
      auto prev   = _delay[pos];
      _delay[pos] = in;

      auto lp = lp_type::tick (
        make_crange (_coeffs).cast (double {}),
        make_crange (_states).cast (double {}),
        in,
        order,
        n_stages,
        _sample_idx);

      // auto hp = prev - lp; // HP

      outs[0][i] = lp[0];
      outs[1][i] = lp[1];
      ++_sample_idx;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  uint            _sample_idx;
  float           _a, _b, _c, _d;

  using lp_type = linear_iir_butterworth_lowpass_any_order;

  static constexpr auto max_order
    = lp_type::linear_iir_butterworth_lowpass_any_order::max_order;

  static constexpr auto order = 2;

  static constexpr uint n_stages  = 12;
  static constexpr uint max_delay = ((1 << n_stages) + 1) * 2;

  std::array<double_x2, lp_type::get_n_coeffs (max_order) * 2> _coeffs;
  std::array<double_x2, lp_type::get_n_states (max_order, n_stages) * 2>
                                       _states;
  std::array<double_x2, max_delay * 2> _delay;
};
//------------------------------------------------------------------------------
#endif
#if 0
// 24 db butterworth linear phase
//------------------------------------------------------------------------------
class experiments {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct param_a_tag {};

  void set (param_a_tag, float v) { _a = v; }

  static constexpr auto get_parameter (param_a_tag)
  {
    return float_param ("A", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v) { _b = v; }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("B", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_c_tag {};

  void set (param_c_tag, float v) { _c = v; }

  static constexpr auto get_parameter (param_c_tag)
  {
    return float_param ("C", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_d_tag {};

  void set (param_d_tag, float v) { _d = v; }

  static constexpr auto get_parameter (param_d_tag)
  {
    return float_param ("D", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<param_a_tag, param_b_tag, param_c_tag, param_d_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    memset (this, 0, sizeof *this);
    _plugcontext = &pc;

    std::array<vec_complex<double_x2>, 2> poles;
    std::array<vec_complex<double_x2>, 2> zeros;

    butterworth_lp_complex::poles (
      make_crange (poles),
      vec_set<double_x2> (660),
      (double) pc.get_sample_rate(),
      2);

    butterworth_lp_complex::zeros (
      make_crange (zeros),
      vec_set<double_x2> (660),
      (double) pc.get_sample_rate(),
      2);

    _gain = butterworth_lp_complex::gain (make_crange (poles), 2);

    t_rev_ccpole_pair_rzero_pair::reset_coeffs (
      make_crange (co_rev_poles).cast (double {}),
      poles[0],
      vec_real (zeros[0]));

    ccpole_pair_rzero_pair::reset_coeffs (
      make_crange (co_poles).cast (double {}), poles[0], vec_real (zeros[0]));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {

      double_x2 out {ins[0][i], ins[1][i]};

      auto pos    = _sample_idx % _delay.size();
      auto prev   = _delay[pos];
      _delay[pos] = out;

      out = t_rev_ccpole_pair_rzero_pair::tick (
        make_crange (co_rev_poles).cast (double {}),
        make_crange (st_rev_poles).cast (double {}),
        out,
        n_stages,
        _sample_idx);
      out *= _gain;

      out = ccpole_pair_rzero_pair::tick (
        make_crange (co_poles).cast (double {}),
        make_crange (st_poles).cast (double {}),
        out);
      out *= _gain;

      out = prev - out;

      outs[0][i] = out[0];
      outs[1][i] = out[1];
      ++_sample_idx;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  double_x2       _gain;
  uint            _sample_idx;
  float           _a, _b, _c, _d;

  static constexpr uint n_stages = 12;
  static constexpr uint delay    = (1 << n_stages) + 1;

  std::array<double_x2, ccpole_pair_rzero_pair::n_coeffs>       co_poles;
  std::array<double_x2, t_rev_ccpole_pair_rzero_pair::n_coeffs> co_rev_poles;

  std::array<double_x2, ccpole_pair_rzero_pair::n_states> st_poles;
  std::array<double_x2, t_rev_ccpole_pair_rzero_pair::get_n_states (n_stages)>
    st_rev_poles;

  std::array<double_x2, delay> _delay;
};
//------------------------------------------------------------------------------
#endif
#if 0
// 12 db butterworth linear phase
//------------------------------------------------------------------------------
class experiments {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct param_a_tag {};

  void set (param_a_tag, float v) { _a = v; }

  static constexpr auto get_parameter (param_a_tag)
  {
    return float_param ("A", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v) { _b = v; }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("B", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_c_tag {};

  void set (param_c_tag, float v) { _c = v; }

  static constexpr auto get_parameter (param_c_tag)
  {
    return float_param ("C", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct param_d_tag {};

  void set (param_d_tag, float v) { _d = v; }

  static constexpr auto get_parameter (param_d_tag)
  {
    return float_param ("D", 0., 1., 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<param_a_tag, param_b_tag, param_c_tag, param_d_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    memset (this, 0, sizeof *this);
    _plugcontext = &pc;

    std::array<vec_complex<double_x2>, 1> poles;
    std::array<vec_complex<double_x2>, 1> zeros;

    butterworth_lp_complex::poles (
      make_crange (poles),
      vec_set<double_x2> (660),
      (double) pc.get_sample_rate(),
      1);

    butterworth_lp_complex::zeros (
      make_crange (zeros),
      vec_set<double_x2> (660),
      (double) pc.get_sample_rate(),
      1);

    _gain = butterworth_lp_complex::gain (make_crange (poles), 1);

    t_rev_rpole_rzero::reset_coeffs (
      make_crange (co_rev_poles).cast (double {}),
      vec_real (poles[0]),
      vec_real (zeros[0]));

    rpole_rzero::reset_coeffs (
      make_crange (co_poles).cast (double {}),
      vec_real (poles[0]),
      vec_real (zeros[0]));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {

      double_x2 out {ins[0][i], ins[1][i]};

      auto pos    = _sample_idx % _delay.size();
      auto prev   = _delay[pos];
      _delay[pos] = out;

      out = t_rev_rpole_rzero::tick (
        make_crange (co_rev_poles).cast (double {}),
        make_crange (st_rev_poles).cast (double {}),
        out,
        n_stages,
        _sample_idx);
      out *= _gain;

      out = rpole_rzero::tick (
        make_crange (co_poles).cast (double {}),
        make_crange (st_poles).cast (double {}),
        out);
      out *= _gain;

      out = prev - out;

      outs[0][i] = out[0];
      outs[1][i] = out[1];
      ++_sample_idx;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  double_x2       _gain;
  uint            _sample_idx;
  float           _a, _b, _c, _d;

  static constexpr uint n_stages = 12;
  static constexpr uint delay    = (1 << n_stages);

  std::array<double_x2, rpole_rzero::n_coeffs>       co_poles;
  std::array<double_x2, t_rev_rpole_rzero::n_coeffs> co_rev_poles;

  std::array<double_x2, rpole_rzero::n_states> st_poles;
  std::array<double_x2, t_rev_rpole_rzero::get_n_states (n_stages)>
    st_rev_poles;

  std::array<double_x2, delay> _delay;
};

#endif
}; // namespace artv
