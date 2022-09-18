#pragma once

#include <cstdio>
// a dummy dsp process to experiment with

#include "artv-common/dsp/own/classes/block_resampler.hpp"

#include "artv-common/dsp/own/classes/linphase_iir_crossover.hpp"

#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/composite/linear_iir_butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros.hpp"
#include "artv-common/dsp/own/parts/filters/poles_zeros_reversed.hpp"

#include "artv-common/dsp/own/classes/circular_queue.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
#if 1
// block_resampler test
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

  void set (param_a_tag, float v)
  {
    if (_a == v) {
      return;
    }
    _a = v;
    reset_src();
  }

  static constexpr auto get_parameter (param_a_tag)
  {
    return float_param ("Ratio", 0.25, 2., 1., 0.25);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v)
  {
    if (_b == v) {
      return;
    }
    _b = v;
    reset_src();
  }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("Minphase", 0., 1., 0., 1.);
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
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    _resampler.process (outs, ins, block_samples, [=] (auto block) {
      process_resampled_block (block);
    });
  }
  //----------------------------------------------------------------------------
private:
  void process_resampled_block (xspan<vec<float, 2>> io)
  {
    // process the block here...
  }
  //----------------------------------------------------------------------------
  void reset_src()
  {
    uint  src_srate = _plugcontext->get_sample_rate();
    uint  tgt_srate = (uint) (_a * (float) src_srate);
    float fc        = 0.46f;
    auto  taps      = 128;
    float fc1       = fc * (float) tgt_srate;
    float fc2       = fc * (float) std::min (tgt_srate, src_srate);

    _resampler.reset (
      tgt_srate,
      src_srate,
      fc1,
      fc2,
      taps,
      taps,
      120.f,
      _b == 1.f,
      16,
      6 * 1024);

    printf ("ratio: %f\n", (float) tgt_srate / (float) src_srate);
    printf ("src_rate(Hz): %u\n", src_srate);
    printf ("tgt_rate(Hz): %u\n", tgt_srate);
    printf ("fc1(Hz): %f\n", fc1);
    printf ("fc2(Hz): %f\n\n", fc2);
  }
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  float           _a, _b, _c, _d;

  static constexpr uint n_channels = 2;

  block_resampler<float, n_channels> _resampler;
};
#endif
#if 0
// 0.5x minphase resampler
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
    return float_param ("Ratio", 0.25, 2., 1., 0.25);
  }
  //----------------------------------------------------------------------------
  struct param_b_tag {};

  void set (param_b_tag, float v) { _b = v; }

  static constexpr auto get_parameter (param_b_tag)
  {
    return float_param ("Minphase", 0., 1., 0., 1.);
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
    std::vector<float> kernel;
    kernel.resize (128);
    kaiser_lp_kernel_2<float> (
      kernel,
      (1. / (float) ratio) * 0.48f,
      kaiser_beta_estimate (120.f),
      0,
      true);
    _decimator.reset (kernel, ratio);
    _interpolator.reset (kernel, ratio);
    _in_bf.clear();
    _in_bf.resize ((ratio - 1) * n_channels);
    _out_bf_mem.clear();
    _out_bf_mem.resize (pow2_round_ceil (ratio * n_channels));
    _out_bf.reset (_out_bf_mem);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    uint n_out = 0;
    for (uint i = 0; i < block_samples; ++i) {
      while (_out_bf.size() && n_out < i) {
        outs[0][n_out] = _out_bf.pop();
        outs[1][n_out] = _out_bf.pop();
        ++n_out;
      }
      _in_bf.push_back (ins[0][i]);
      _in_bf.push_back (ins[1][i]);
      if (_in_bf.size() == (ratio * n_channels)) {
        auto dspl = _decimator.tick (_in_bf);
        _in_bf.clear();
        std::array<float, ratio * n_channels> uspls;
        _interpolator.tick (uspls, dspl);
        _out_bf.push (uspls);
      }
    }
    while (_out_bf.size() && n_out < block_samples) {
      outs[0][n_out] = _out_bf.pop();
      outs[1][n_out] = _out_bf.pop();
      ++n_out;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  plugin_context* _plugcontext;
  float           _a, _b, _c, _d;

  static constexpr uint n_channels = 2;
  static constexpr uint ratio      = 2;

  fir_decimator<float, 2>           _decimator;
  fir_interpolator<float, 2>        _interpolator;
  std::vector<float>                _in_bf;
  std::vector<float>                _out_bf_mem;
  static_pow2_circular_queue<float> _out_bf;
};
#endif
#if 0
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
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
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
      xspan {_coeffs}.cast (double {}),
      vec_set<double_x2> (660.),
      pc.get_sample_rate(),
      order);
    // pc.set_delay_compensation (lp_type::get_latency (order, n_stages));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {
      double_x2 in {double_x2 {ins[0][i], ins[1][i]}};

      auto pos    = _sample_idx % lp_type::get_latency (order, n_stages);
      auto prev   = _delay[pos];
      _delay[pos] = in;

      auto lp = lp_type::tick (
        xspan {_coeffs}.cast (double {}),
        xspan {_states}.cast (double {}),
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

    double t_spl =  1. / (double) pc.get_sample_rate();

    butterworth_lp_complex::poles (
      xspan {poles},
      vec_set<double_x2> (660),
      t_spl,
      2);

    butterworth_lp_complex::zeros (
      xspan {zeros},
      vec_set<double_x2> (660),
      t_spl,
      2);

    _gain = butterworth_lp_complex::gain (xspan {poles}, 2);

    t_rev_ccpole_pair_rzero_eq_pair::reset_coeffs (
      xspan {co_rev_poles}.cast (double {}),
      poles[0],
      vec_real (zeros[0]));

    ccpole_pair_rzero_pair::reset_coeffs (
      xspan {co_poles}.cast (double {}), poles[0], vec_real (zeros[0]));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {

      double_x2 out {ins[0][i], ins[1][i]};

      auto pos    = _sample_idx % _delay.size();
      auto prev   = _delay[pos];
      _delay[pos] = out;

      out = t_rev_ccpole_pair_rzero_eq_pair::tick (
        xspan {co_rev_poles}.cast (double {}),
        xspan {st_rev_poles}.cast (double {}),
        out,
        n_stages,
        _sample_idx);
      out *= _gain;

      out = ccpole_pair_rzero_pair::tick (
        xspan {co_poles}.cast (double {}),
        xspan {st_poles}.cast (double {}),
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
  std::array<double_x2, t_rev_ccpole_pair_rzero_eq_pair::n_coeffs> co_rev_poles;

  std::array<double_x2, ccpole_pair_rzero_pair::n_states> st_poles;
  std::array<double_x2, t_rev_ccpole_pair_rzero_eq_pair::get_n_states (n_stages)>
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

    double t_spl =  1. / (double) pc.get_sample_rate();

    butterworth_lp_complex::poles (
      xspan {poles},
      vec_set<double_x2> (660),
      t_spl,
      1);

    butterworth_lp_complex::zeros (
      xspan {zeros},
      vec_set<double_x2> (660),
      t_spl,
      1);

    _gain = butterworth_lp_complex::gain (xspan {poles}, 1);

    t_rev_rpole_rzero::reset_coeffs (
      xspan {co_rev_poles}.cast (double {}),
      vec_real (poles[0]),
      vec_real (zeros[0]));

    rpole_rzero::reset_coeffs (
      xspan {co_poles}.cast (double {}),
      vec_real (poles[0]),
      vec_real (zeros[0]));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    for (uint i = 0; i < block_samples; ++i) {

      double_x2 out {ins[0][i], ins[1][i]};

      auto pos    = _sample_idx % _delay.size();
      auto prev   = _delay[pos];
      _delay[pos] = out;

      out = t_rev_rpole_rzero::tick (
        xspan {co_rev_poles}.cast (double {}),
        xspan {st_rev_poles}.cast (double {}),
        out,
        n_stages,
        _sample_idx);
      out *= _gain;

      out = rpole_rzero::tick (
        xspan {co_poles}.cast (double {}),
        xspan {st_poles}.cast (double {}),
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
