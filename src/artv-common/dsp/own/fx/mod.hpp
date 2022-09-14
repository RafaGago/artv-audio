#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/biquad.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/zdf.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/own/parts/waveshapers/sigmoid.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#define MOD_DBG_DENORMALS 0
#if MOD_DBG_DENORMALS
#include <fenv.h>
#endif
namespace artv {
//------------------------------------------------------------------------------
class mod {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::modulation;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct order_tag {};
  void set (order_tag, float v) { _param_smooth.target().order = v * 0.01f; }

  static constexpr auto get_parameter (order_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_rate_tag {};
  void set (lfo_rate_tag, float v)
  {
    _param_smooth.target().lfo_rate = v;
    _param.lfo_off                  = (v == 100.f);
  }

  static constexpr auto get_parameter (lfo_rate_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_time_base_tag {};
  void set (lfo_time_base_tag, int v) { _param.lfo_time_base = v; }

  struct lfo_time_base {
    enum { free, quarter_beat, two_beats };
  };

  static constexpr auto get_parameter (lfo_time_base_tag)
  {
    return choice_param (
      0, make_cstr_array ("Free", "1/4 beat=10%", "2 beats=10%"), 16);
  }
  //----------------------------------------------------------------------------
  struct lfo_depth_tag {};
  void set (lfo_depth_tag, float v)
  {
    v *= 0.01;
    _param_smooth.target().lfo_depth = v;
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_wave_tag {};
  void set (lfo_wave_tag, int v) { _param.lfo_wave = v; }

  struct lfo_wave {
    enum { sine, triangle, sample_hold, noise, trapezoid, square, saw };
  };

  static constexpr auto get_parameter (lfo_wave_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Sine", "Triangle", "S&H", "Noise", "Trapezoid", "Square", "Saw"),
      10);
  }
  //----------------------------------------------------------------------------
  struct lfo_stereo_tag {};
  void set (lfo_stereo_tag, float v) { _param_smooth.target().lfo_stereo = v; }

  static constexpr auto get_parameter (lfo_stereo_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct center_tag {};
  void set (center_tag, float v) { _param_smooth.target().center = v * 0.01; }

  static constexpr auto get_parameter (center_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct a_tag {};
  void set (a_tag, float v) { _param_smooth.target().a = v * 0.01; }

  static constexpr auto get_parameter (a_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01;
    v = sqrt (sqrt (abs (v))); // something cheaper?
    v = neg ? -v : v;

    _param_smooth.target().feedback = v * 0.98f;
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", -100., 100., 50., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_locut_tag {};
  void set (feedback_locut_tag, float v)
  {
    v *= 0.01;
    if (v != _param.feedback_locut) {
      _param.feedback_locut = v;
      _fb_filters.reset_coeffs<fb_filter::locut> (
        vec_set<4> (4.f + 620.f * v * v),
        vec_set<4> ((float) M_SQRT1_2),
        _t_spl);
    }
  }

  static constexpr auto get_parameter (feedback_locut_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_hicut_tag {};
  void set (feedback_hicut_tag, float v)
  {
    v *= 0.01;
    if (v != _param.feedback_hicut) {
      _param.feedback_hicut = v;
      _fb_filters.reset_coeffs<fb_filter::hicut> (
        vec_set<4> (21000.f - 19000.f * v * v),
        vec_set<4> ((float) M_SQRT1_2),
        _t_spl);
    }
  }

  static constexpr auto get_parameter (feedback_hicut_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_drive_tag {};
  void set (feedback_drive_tag, float v)
  {
    v *= 0.01f;
    v                                     = v * v * v * v;
    _param_smooth.target().feedback_drive = v * 1500.f;
  }

  static constexpr auto get_parameter (feedback_drive_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_curve_tag {};
  void set (feedback_curve_tag, float v)
  {
    _param_smooth.target().feedback_curve = v *= 0.01f;
  }

  static constexpr auto get_parameter (feedback_curve_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct b_tag {};
  void set (b_tag, float v) { _param_smooth.target().b = v * 0.01f; }

  static constexpr auto get_parameter (b_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct detune_tag {};
  void set (detune_tag, float v) { _param_smooth.target().detune = v * 0.01f; }

  static constexpr auto get_parameter (detune_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct depth_tag {};
  void set (depth_tag, float v)
  {
    v *= 0.01f;
    _param_smooth.target().depth = v;
  }

  static constexpr auto get_parameter (depth_tag)
  {
    return float_param ("%", 0., 100, 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    // TODO: initializations?
    _param.mode = v;
  }

  struct mode {
    enum { zdf, schroeder };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Phaser (ZDF)", "Phaser (Schroeder)"), 16);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _srate               = pc.get_sample_rate();
    _t_spl               = 1.f / _srate;
    _beat_hz             = pc.get_play_state().bpm * (1.f / 60.f);
    _n_processed_samples = 0;
    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _lfo_t_spl           = _t_spl * (_control_rate_mask + 1);
    _1spl_fb             = decltype (_1spl_fb) {};

    _phaser.reset_states_cascade();
    _fb_filters.reset_states_cascade();
    _feedback.reset_states_cascade();
    _onepole.reset_states_cascade();
    _param_smooth.reset (_t_spl, 10.f);

    // ensure that all parameter's new values will be different
    memset (&_param, -1, sizeof _param);
    memset (&_param_smooth.target(), -1, sizeof _param_smooth.target());

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _param_smooth.set_all_from_target();

    reset_mem();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
#if MOD_DBG_DENORMALS
    feenableexcept (FE_INVALID);
#endif
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    switch (_param.mode) {
    case mode::zdf:
      tick_phaser<T> (outs, ins, samples);
      break;
    case mode::schroeder:
      tick_schroeder<T> (outs, ins, samples);
      break;
    default:
      assert (false);
      break;
    }
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    order_tag,
    lfo_rate_tag,
    lfo_time_base_tag,
    lfo_depth_tag,
    lfo_wave_tag,
    lfo_stereo_tag,
    center_tag,
    a_tag,
    feedback_tag,
    feedback_locut_tag,
    feedback_hicut_tag,
    feedback_drive_tag,
    feedback_curve_tag,
    detune_tag,
    b_tag,
    depth_tag,
    mode_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr uint max_phaser_stages = 16;
  static constexpr uint max_scho_stages   = 12;
  static constexpr uint max_scho_delay_ms = 30;
  static constexpr uint n_channels        = 2;
  static constexpr uint n_ap_channels     = 4;
  //----------------------------------------------------------------------------
  struct all_parameters;
  using lfo_type       = lfo<n_channels>;
  using lfo_value_type = lfo_type::value_type;
  //----------------------------------------------------------------------------
  lfo_value_type run_lfo (all_parameters const& pars)
  {
    float hz;
    switch (pars.lfo_time_base) {
    case lfo_time_base::free: {
      float lfo_rate = pars.lfo_rate * 0.01f;
      lfo_rate       = 1.f - lfo_rate;
      lfo_rate *= lfo_rate * lfo_rate;
      hz = lfo_rate * 20.f; // 0 to 20 Hz
    } break;
    case lfo_time_base::quarter_beat:
      // a quarter beat for each 10%
      hz = (4.f * _beat_hz) / (1.f + pars.lfo_rate * 0.1f);
      break;
    case lfo_time_base::two_beats:
      // two beats for each 10%
      hz = (1 / 2.f * _beat_hz) / (1.f + pars.lfo_rate * 0.1f);
      break;
    default:
      assert (false);
      break;
    }

    _lfos.set_freq (vec_set<2> (hz), _lfo_t_spl);
    auto stereo_ph
      = phase<1> {vec_set<1> (pars.lfo_stereo), phase<1>::degrees {}}.get_raw (
        0);
    if (!pars.lfo_off) {
      // updating stereo phase diff.
      auto ph = _lfos.get_phase();
      ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (ph);
    }
    else {
      auto start_ph = phase<n_channels> {
        vec_set<n_channels> (0.f), phase<n_channels>::degrees {}};
      start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (start_ph);
    }

    lfo_value_type lfov;
    switch (pars.lfo_wave) {
    case lfo_wave::sine:
      lfov = _lfos.tick_sine();
      break;
    case lfo_wave::triangle:
      lfov = _lfos.tick_triangle();
      break;
    case lfo_wave::sample_hold:
      lfov = _lfos.tick_sample_hold();
      break;
    case lfo_wave::noise:
      lfov = _lfos.tick_filt_sample_and_hold();
      break;
    case lfo_wave::trapezoid:
      lfov = _lfos.tick_trapezoid (vec_set<2> (0.3f));
      break;
    case lfo_wave::square:
      lfov = _lfos.tick_square();
      break;
    case lfo_wave::saw:
      lfov = _lfos.tick_saw();
      break;
    }
    return lfov;
  }
  //----------------------------------------------------------------------------
  void run_phaser_mod (
    all_parameters const& pars,
    uint                  n_stages,
    float                 detune_ratio,
    lfo_value_type        lfo)
  {
    float modlim = std::min (pars.center, 0.5f);
    modlim       = std::min (1.f - pars.center, modlim);

    lfo_value_type mod    = lfo * pars.lfo_depth * modlim;
    lfo_value_type center = pars.center + mod;
    // detune by a ramp
    lfo_value_type detuned = (1.f - abs (pars.detune * detune_ratio)) * center;
    lfo_value_type diff    = (center - detuned) / (float) (n_stages * 2);
    lfo_value_type start   = (pars.detune > 0.f) ? detuned : center;
    diff                   = (pars.detune > 0.f) ? diff : -diff;
    float_x4 curr          = vec_cat (start, start + diff);
    float_x4 add           = vec_cat (diff, diff) * 2.f;

    for (uint i = 0; i < n_stages; ++i) {
      _mod[i] = vec_to_array (curr);
      curr += add;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_phaser (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {
      _param_smooth.tick();
      // access parameters without caring if they are smoothed or not.
      all_parameters pars;
      *((smoothed_parameters*) &pars)   = _param_smooth.get();
      *((unsmoothed_parameters*) &pars) = _param;

      if ((_n_processed_samples & _control_rate_mask) == 0) {
        _n_stages = 1 + (uint) (pars.order * (max_phaser_stages - 1));
        run_phaser_mod (pars, _n_stages, 0.07f, run_lfo (pars));
        for (uint s = 0; s < _n_stages; ++s) {
          float_x4 v = vec_from_array (_mod[s]);
          // cheap'ish parametric curve
          constexpr float k = 0.8f; // TODO curve: k = parameter
          auto            f = ((-1.f + v) / (1.f - k * v)) + 1.f;
          // auto f = v * v;
          auto q = f;
          f      = ((f * 0.993f) + 0.007f);
          f *= float_x4 {21200.f, 21200.f, 20000.f, 20000.f};
          // A lot of freq-dependant weight on the Q when detuning
          q = (1.f + 7.f * q * abs (pars.detune));
          q *= pars.a * pars.a * float_x4 {3.4f, 3.4f, 3.8f, 3.8f};
          q += 0.05f;
          _phaser.reset_coeffs_on_idx (s, f, q, _t_spl);
        }
        float_x4 f {450.f, 460.f, 670.f, 680.f};
        f *= vec_from_array (_mod[0]);
        _onepole.reset_coeffs<0> (f, _t_spl);
      }

      float_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      auto     onep_out = _onepole.tick<0> (wet);
      wet               = vec_shuffle (wet, onep_out, 0, 1, 4, 5);

      // Run allpass cascades with own feedback loop and saturation
      std::array<float_x4, max_phaser_stages * zdf::n_gs_coeffs> G_S_mem;

      // obtain G and S for the phaser
      auto gs = make_crange (G_S_mem);
      for (uint s = 0; s < _n_stages; ++s) {
        _phaser.tick_on_idx (
          s, gs.cut_head (zdf::n_gs_coeffs), zdf::gs_coeffs_tag {});
      }
      auto aps_resp = zdf::combine_response<float_x4> (
        make_crange (G_S_mem.data(), _n_stages * zdf::n_gs_coeffs));
      // obtain G and S for the filters
      gs = make_crange (G_S_mem);
      _fb_filters.tick<fb_filter::locut> (
        gs.cut_head (zdf::n_gs_coeffs), zdf::gs_coeffs_tag {});
      _fb_filters.tick<fb_filter::hicut> (
        gs.cut_head (zdf::n_gs_coeffs), zdf::gs_coeffs_tag {});
      auto filt_resp = zdf::combine_response<float_x4> (
        make_crange (G_S_mem.data(), fb_filter::count * zdf::n_gs_coeffs));

      // Get ZDF feedback
      wet = _feedback.tick (
        wet,
        aps_resp,
        filt_resp,
        vec_set<4> (pars.feedback),
        vec_set<4> (pars.feedback_drive),
        vec_set<4> (pars.feedback_curve));

      // Run regular filter processing
      assert (_n_stages > 0);
      for (uint s = 0; s < _n_stages; ++s) {
        wet = _phaser.tick_on_idx (s, wet);
      }
      // just running the shelves to update the states.
      auto fbv = wet * pars.feedback;
      fbv      = zdf_type::nonlin::tick (
        fbv,
        vec_set<4> (pars.feedback_drive),
        vec_set<4> (pars.feedback_curve));
      _fb_filters.tick_cascade (fbv);

      auto out   = mix (pars, wet, ins[0][i], ins[0][1]);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_schroeder (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {
      _param_smooth.tick();
      // access parameters without caring if they are smoother or not.
      all_parameters pars;
      *((smoothed_parameters*) &pars)   = _param_smooth.get();
      *((unsmoothed_parameters*) &pars) = _param;

      if ((_n_processed_samples & _control_rate_mask) == 0) {
        _n_stages = 1 + (uint) (pars.order * (max_scho_stages - 1));

        pars.center = 1.f - pars.center; // reverse range lf to hf
        pars.center /= _n_stages;
        run_phaser_mod (pars, _n_stages, 0.25f, run_lfo (pars));

        for (uint s = 0; s < _n_stages; ++s) {
          constexpr float sec_factor = (max_scho_delay_ms * 0.001f);
          float_x4        v          = vec_from_array (_mod[s]);
          float           t          = v[0];
          t                          = _srate * (t * t * sec_factor);
          _scho_del_spls.target()[s] = t;
          constexpr float min_g      = 0.0001f;
          constexpr float g_factor   = 0.15 - min_g;
          _scho_g[s]
            = min_g + pars.a * pars.a * 0.8f + v * v * g_factor * pars.a;
        }
        for (uint s = _n_stages; s < max_scho_stages; ++s) {
          _scho_del_spls.target()[s] = _scho_del_spls.target()[s - 1];
        }
        float_x4 fv {1460.f, 1460.f, 1670.f, 1780.f};
        auto     f = vec_from_array (_mod[_n_stages / 2]);
        f *= f;
        f *= (1.f + pars.a) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      _scho_del_spls.tick();

      float_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      std::array<float_x4, max_scho_stages> to_push {};
      auto&                                 del_spls = _scho_del_spls.get();
      for (uint s = 0; s < _n_stages; ++s) {
        auto n_spls = std::clamp (
          del_spls[s],
          (float) _scho.min_delay_spls(),
          (float) _scho.max_delay_spls());
        auto yn    = _scho.get (n_spls, s);
        auto r     = allpass_fn::tick<float_x4> (wet, yn, _scho_g[s]);
        wet        = r.out;
        to_push[s] = r.to_push;
      }
      _scho.push (to_push);
      auto fb_val = zdf_type::nonlin::tick (
        wet,
        vec_set<4> (pars.feedback_drive),
        vec_set<4> (pars.feedback_curve));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);

      auto out   = mix (pars, wet, ins[0][i], ins[1][i]);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  double_x2 mix (all_parameters const& p, float_x4 wet, T dry_l, T dry_r)
  {
    float dry_gain = 1.f - p.depth;
    float wet_gain = p.depth;
    wet_gain /= (get_fb_gain (p.feedback, p.feedback_drive));

    double_x2 wetdbl   = {(double) wet[0], (double) wet[1]};
    double_x2 parallel = {(double) wet[2], (double) wet[3]};
    double_x2 dry      = {(double) dry_l, (double) dry_r};
    wetdbl *= 1.f - abs (p.b);
    wetdbl += p.b * parallel;
    return wetdbl * wet_gain + dry * dry_gain;
  }
  //----------------------------------------------------------------------------
  // hardness makes a sqrt sigmmoid
  float get_fb_gain (float fb, float hardness)
  {
    float lim = zdf_type::nonlin::limit_inf<float_x1> (
      make_vec (hardness), make_vec (0.f))[0];
    return 1.f + std::min (lim, abs (fb));
  }
  //----------------------------------------------------------------------------
  void reset_mem()
  {
    uint schroeder_size
      = std::ceil (max_scho_delay_ms * 0.001f * _srate) + _scho.min_size_spls();
    schroeder_size = _scho.n_required_elems (schroeder_size, max_scho_stages);
    _mem.clear();
    _mem.resize (schroeder_size);
    _scho.reset (_mem, max_scho_stages);
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   lfo_wave;
    u32   lfo_time_base;
    u32   mode;
    float feedback_locut;
    float feedback_hicut;
    bool  lfo_off; // The lfo is smoothed, this signals if the lfo is off
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float lfo_rate;
    float lfo_depth;
    float lfo_stereo;
    float center;
    float a;
    float b;
    float depth; // aka mix (make explicit?)
    float feedback;
    float feedback_drive;
    float feedback_curve;
    float detune;
    float order;
  };
  //----------------------------------------------------------------------------
  struct all_parameters : public unsmoothed_parameters,
                          public smoothed_parameters {};
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  using allpass_type = andy::svf_zdf_allpass;
  using filters_list = mp_list<andy::svf_zdf_highpass, andy::svf_zdf_lowpass>;
  using zdf_type     = zdf::feedback<
    zdf::lin_pre_fb_node_nonlin_after_tag,
    zdf::lin_mystran_tag<2>,
    sigmoid::mystran<6, true>>;
  struct fb_filter {
    enum { locut, hicut, count };
  };

  part_classes<mp_list<zdf_type>, float_x4>                     _feedback;
  part_class_array<allpass_type, float_x4, max_phaser_stages>   _phaser;
  part_classes<filters_list, float_x4>                          _fb_filters;
  part_classes<mp_list<onepole_allpass>, float_x4>              _onepole;
  interpolated_delay_line<float_x4, linear_interp, true, false> _scho;
  value_smoother<float, std::array<float, max_scho_stages>>     _scho_del_spls;
  std::array<float_x4, max_scho_stages>                         _scho_g;

  alignas (sse_bytes) array2d<float, n_ap_channels, max_phaser_stages> _mod;

  std::vector<float_x4, overaligned_allocator<float_x4, 128>> _mem;
  lfo_type _lfos; // 0 = L, 1 = R
  float_x4 _1spl_fb;
  uint     _n_processed_samples;
  uint     _control_rate_mask;
  uint     _n_stages;
  float    _lfo_t_spl;
  float    _t_spl;
  float    _srate;
  float    _lp_smooth_coeff;
  float    _beat_hz;
};
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
#if 0
//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::modulation;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct stages_tag {};
  void set (stages_tag, int v) { _params.unsmoothed.n_stages = v + 1; }

  static constexpr auto get_parameter (stages_tag)
  {
    // TODO: before release, not so many choices make a difference, reduce to
    // e.g. 2, 4, 6, 8, 12, 16, 24, 32
    return choice_param (
      6,
      make_cstr_array (
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16"),
      32);
  }
  //----------------------------------------------------------------------------
  struct stages_mode_tag {};
  void set (stages_mode_tag, int v)
  {
    bool lin                       = v < m_total_modes;
    _params.unsmoothed.lin_lfo_mod = lin;
    _params.unsmoothed.mode        = v - (!lin * m_total_modes);
  }

  static constexpr auto get_parameter (stages_mode_tag)
  {
    return choice_param (
      8,
      make_cstr_array (
        "Lin Mod Exp F",
        "Lin Mod Exp F Spread",
        "Lin Mod Exp F Alternate",
        "Lin Mod Exp F Stereo",
        "Lin Mod Lin F",
        "Lin Mod Lin F Spread",
        "Lin Mod Lin F Alternate",
        "Lin Mod Lin F Stereo",
        "Exp Mod Exp F",
        "Exp Mod Exp F Spread",
        "Exp Mod Exp F Alternate",
        "Exp Mod Exp F Stereo",
        "Exp Mod Lin F",
        "Exp Mod Lin F Spread",
        "Exp Mod Lin F Alternate",
        "Exp Mod Lin F Stereo"),
      16);
  }

  enum stage_mode {
    m_exp,
    m_exp_spread,
    m_exp_alternate,
    m_exp_stereo,
    m_lin,
    m_lin_spread,
    m_lin_alternate,
    m_lin_stereo,
    m_total_modes,
  };

  static constexpr bool mode_is_exponential (uint v) { return v < m_lin; }
  //----------------------------------------------------------------------------
  struct lfo_rate_tag {};
  void set (lfo_rate_tag, float v)
  {
    v                              = midi_note_to_hz (v);
    _params.unsmoothed.lfo_hz_user = v;
    refresh_lfo_hz();
  }

  static constexpr auto get_parameter (lfo_rate_tag)
  {
    return frequency_parameter_from_zero (15.0, 0.66);
  }
  //----------------------------------------------------------------------------
  struct lfo_rate_sync_tag {};
  void set (lfo_rate_sync_tag, float v)
  {
    _params.unsmoothed.lfo_eights = v;
    refresh_lfo_hz();
  }

  static constexpr auto get_parameter (lfo_rate_sync_tag)
  {
    // NOTE: this has no 1/8 displayed because it breaks the display pane on
    // mix-matrix because of some (Juce?) bug.
    return float_param ("eights", 0., 32., 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct lfo_depth_tag {};
  void set (lfo_depth_tag, float v)
  {
    v *= 0.01;
    _params_smooth.target().lfo_depth = v;
  }

  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_wave_tag {};
  void set (lfo_wave_tag, int v) { _params.unsmoothed.lfo_wave = v; }

  static constexpr auto get_parameter (lfo_wave_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Sine", "Triangle", "S&H", "Noise", "Trapezoid", "Square", "Saw"),
      10);
  }
  //----------------------------------------------------------------------------
  struct lfo_start_phase_tag {};
  void set (lfo_start_phase_tag, float v)
  {
    _params_smooth.target().lfo_start_phase = v;
  }

  static constexpr auto get_parameter (lfo_start_phase_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lfo_stereo_tag {};
  void set (lfo_stereo_tag, float v)
  {
    _params_smooth.target().lfo_stereo = v;
  }

  static constexpr auto get_parameter (lfo_stereo_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct low_freq_tag {};
  void set (low_freq_tag, float v)
  {
    v                                   = midi_note_to_hz (v);
    _params_smooth.target().freq_lo = v;
  }

  static constexpr auto get_parameter (low_freq_tag)
  {
    return frequency_parameter (min_ap_freq, 10000., 225.);
  }
  //----------------------------------------------------------------------------
  struct high_freq_tag {};
  void set (high_freq_tag, float v)
  {
    v                                   = midi_note_to_hz (v);
    _params_smooth.target().freq_hi = v;
  }

  static constexpr auto get_parameter (high_freq_tag)
  {
    return frequency_parameter (min_ap_freq, 10000., 3300.);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_feedback = 0.995;
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01 * max_feedback;
    v = sqrt (abs (v));
    v = neg ? -v : v;

    _params_smooth.target().feedback = v;
    scale_feedbacks (
      _params_smooth.target().scaled_feedback,
      _params_smooth.target().scaled_delay_feedback,
      _params_smooth.target().feedback,
      _params.unsmoothed.delay_feedback);
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_hp_tag {};
  void set (feedback_hp_tag, float v)
  {
    _params_smooth.target().feedback_hp = v * 0.01f;
  }

  static constexpr auto get_parameter (feedback_hp_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_lp_tag {};
  void set (feedback_lp_tag, float v)
  {
    _params_smooth.target().feedback_lp = v * 0.01f;
  }

  static constexpr auto get_parameter (feedback_lp_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_drv_tag {};
  void set (feedback_drv_tag, float v)
  {
    v *= 0.01f;
    v = v * v * v;
    v *= 100.f;
    _params_smooth.target().feedback_drv = v;
  }

  static constexpr auto get_parameter (feedback_drv_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct q_tag {};
  void set (q_tag, float v) { _params_smooth.target().q = v; }

  static constexpr auto get_parameter (q_tag)
  {
    return float_param ("", 0.02, 3.5, 0.18, 0.001, 0.4);
  }
  //----------------------------------------------------------------------------
  struct parallel_mix_tag {};
  void set (parallel_mix_tag, float v)
  {
    _params_smooth.target().parallel_mix = v * 0.005f;
  }

  static constexpr auto get_parameter (parallel_mix_tag)
  {
    return float_param ("%", -100., 100, 100., 0.1);
  }
  //----------------------------------------------------------------------------
  struct topology_tag {};
  void set (topology_tag, int v)
  {
    if (_params.unsmoothed.topology == v) {
      return;
    }
    switch (v) {
    case 0:
    case 1:
      _params.unsmoothed.topology = v;
      break;
    case 2:
      _params.unsmoothed.topology = t_1_pole;
      break;
    case 3:
      _params.unsmoothed.topology = t_2_pole;
      break;
    case 4:
      _params.unsmoothed.topology = t_3_pole;
      break;
    case 5:
      _params.unsmoothed.topology = t_schroeder;
      break;
    case 6:
      _params.unsmoothed.topology = t_comb;
      break;
    case 7:
      _params.unsmoothed.topology = t_hybrid_1;
      break;
    case 8:
      _params.unsmoothed.topology = t_hybrid_2;
      break;
    }

    _n_processed_samples = 0; // force coeffient recalculation
  }

  static constexpr auto get_parameter (topology_tag)
  {
    return choice_param (
      3,
      make_cstr_array (
        "Phaser 1 pole (legacy)",
        "Phaser 2 pole (legacy)",
        "Phaser 1 pole",
        "Phaser 2 pole",
        "Phaser 3 pole",
        "Serial Schroeder AP",
        "Parallel Feedforward Combs",
        "Hybrid 1",
        "Hybrid 2"),
      16);
  }

  // don't change the order, it is abused as an optimization to save comparsions
  enum topologies {
    t_1_pole_legacy,
    t_2_pole_legacy,
    t_3_pole,
    t_hybrid_1,
    t_hybrid_2,
    t_2_pole,
    t_1_pole,
    t_schroeder,
    t_comb,
    t_first_non_legacy             = t_3_pole,
    t_first_non_legacy_non_2_pole  = t_1_pole,
    t_first_non_legacy_non_allpass = t_schroeder,
  };
  //----------------------------------------------------------------------------
  struct delay_feedback_tag {};
  void set (delay_feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01 * 0.96;
    v = sqrt (abs (v));
    v = neg ? -v : v;

    _params.unsmoothed.delay_feedback = v;
    scale_feedbacks (
      _params_smooth.target().scaled_feedback,
      _params_smooth.target().scaled_delay_feedback,
      _params_smooth.target().feedback,
      _params.unsmoothed.delay_feedback);
  }

  static constexpr auto get_parameter (delay_feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.01);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_delay_sec = 0.0015f;

  struct delay_time_tag {};
  void set (delay_time_tag, float v)
  {
    v *= 0.01f;
    v *= max_delay_sec * _srate;
    _params_smooth.target().delay_spls = v;
  }

  static constexpr auto get_parameter (delay_time_tag)
  {
    return float_param ("%", 0., 100, 33., 0.001);
  }
  //----------------------------------------------------------------------------
  static constexpr float max_delay_factor    = 4.f;
  static constexpr float ln_max_delay_factor = M_LN2 * 2.f;

  struct delay_lfo_tag {};
  void set (delay_lfo_tag, float v)
  {
    v *= 0.01f;
    _params_smooth.target().delay_lfo = v;
  }

  static constexpr auto get_parameter (delay_lfo_tag)
  {
    return float_param ("%", -100., 100, 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mix_tag {};
  void set (mix_tag, float v)
  {
    v *= 0.01f;
    _params_smooth.target().mix = v;
  }

  static constexpr auto get_parameter (mix_tag)
  {
    return float_param ("%", 0., 100, 100., 0.001);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _srate       = pc.get_sample_rate();
    _t_spl       = 1.f / _srate;
    _allpass1p.reset_states_cascade();
    _allpass2p.reset_states_cascade();
    _allpass1p_zdf.reset_states_cascade();
    _zdf_ap.reset_states_cascade();
    _zdf_lshelf.reset_states_cascade();
    _zdf_hshelf.reset_states_cascade();
    _feedback_shelf.reset_states_cascade();
    _feedback_delay_shelf.reset_states_cascade();
    _lfos.reset();

    _params_smooth.target().freq_lo
      = get_parameter (low_freq_tag {}).defaultv;
    _params_smooth.target().freq_hi
      = get_parameter (high_freq_tag {}).defaultv;
    _params_smooth.target().lfo_stereo
      = get_parameter (lfo_stereo_tag {}).defaultv;
    _params_smooth.target().lfo_depth
      = get_parameter (lfo_depth_tag {}).defaultv;
    _params_smooth.target().lfo_start_phase
      = get_parameter (lfo_start_phase_tag {}).defaultv;
    _params_smooth.target().feedback
      = get_parameter (feedback_tag {}).defaultv;
    _params_smooth.target().q          = get_parameter (q_tag {}).defaultv;
    _params_smooth.target().delay_spls = 0.f;
    _params_smooth.target().scaled_delay_feedback = 0.f;
    _params_smooth.target().scaled_feedback       = 0.f;
    _params_smooth.target().delay_lfo             = 0.f;
    _params_smooth.target().lfo_last              = 0.f;
    _params_smooth.target().feedback_hp
      = get_parameter (feedback_hp_tag {}).defaultv;
    _params_smooth.target().feedback_hp
      = get_parameter (feedback_lp_tag {}).defaultv;
    _params_smooth.target().feedback_drv
      = get_parameter (feedback_drv_tag {}).defaultv;
    _params_smooth.target().parallel_mix
      = get_parameter (parallel_mix_tag {}).defaultv;
    _params_smooth.target().mix = get_parameter (mix_tag {}).defaultv;

    _params.unsmoothed.lfo_hz_user = get_parameter (lfo_rate_tag {}).defaultv;
    _params.unsmoothed.lfo_eights
      = get_parameter (lfo_rate_sync_tag {}).defaultv;
    _params.unsmoothed.lfo_wave = get_parameter (lfo_wave_tag {}).defaultv;
    _params.unsmoothed.n_stages = get_parameter (stages_tag {}).defaultv;
    _params.unsmoothed.mode     = get_parameter (stages_mode_tag {}).defaultv;

    refresh_lfo_hz();

    memset (&_params.smooth_state, 0, sizeof _params.smooth_state);

    using x1_type = vec<decltype (_lp_smooth_coeff), 1>;

    onepole_smoother::reset_coeffs (
      make_crange (_lp_smooth_coeff).cast (x1_type {}),
      vec_set<x1_type> (1. / 0.01),
      _t_spl);

    _dc_blocker.reset_states();
    _dc_blocker.reset_coeffs (vec_set<4> (2.f), _t_spl);

    _n_processed_samples = 0;

    // Sample rates 44100 multiples update every 362.811us
    uint sr_order      = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask = lsb_mask<uint> (sr_order);
    _feedback_samples  = vec_set<double_x2> (0.);
    _lfo_t_spl         = _t_spl * (_control_rate_mask + 1);

    // setting up delay
    reset_memory_parts();
    _delay.set_resync_delta (10.0);

    // initialize windowed sinc table
    _ff_comb.reset_interpolator (0, false, 0.45f, 140.f);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    if (_params.unsmoothed.topology >= t_first_non_qlegacy) {
      process_zdf (outs, ins, samples);
    }
    else {
      process_legacy (outs, ins, samples);
    }
#if 0
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {

      constexpr uint vec_size = vec_traits<float_x4>().size;

      // parameter smoothing (could be done at a subblock size).
      auto smooth_iter = _params.smooth_target.arr.size() / vec_size;
      for (uint j = 0; j < smooth_iter; ++j) {
        // HACKish encapsulation violation: this one isn't assigning the result,
        // but getting the -1 state from the one pole filter directly (copied by
        // value). See "static_assert" on the code below. 1-pole digital filters
        // are NEVER going to change so this is IMO acceptable.
        auto in
          = make_crange (_params.smooth_target.arr, vec_size, vec_size * j);
        onepole_smoother::tick (
          make_crange (_lp_smooth_coeff),
          make_crange (_params.smooth_state.arr, vec_size, vec_size * j)
            .cast (float_x4 {}),
          vec_load<float_x4> (in));
      }

      all_parameters pars;
      static_assert (
        onepole_smoother::n_states == 1,
        "For the assignment below to work only the Z-1 state has to exist");
      *((smoothed_parameters*) &pars)   = _params.smooth_state.value;
      *((unsmoothed_parameters*) &pars) = _params.unsmoothed;

      // control rate/mdoulation block -----------------------------------------
      if ((_n_processed_samples & _control_rate_mask) == 0) {
        run_control_rate (pars);
      }
      // sample-wise processing block (unbreakable, it has 1 sample delays) ----
      float_x4 fb = {
        (float) _feedback_samples[0],
        (float) _feedback_samples[1],
        (float) _feedback_samples[0],
        (float) _feedback_samples[1]};
      // regular processing
      float_x4 out {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      out += fb * pars.scaled_feedback;
      double_x2 feedforward {out[0], out[1]};

      auto n_spls = pars.delay_spls;
      n_spls *= exp (ln_max_delay_factor * pars.lfo_last * pars.delay_lfo);
      auto delayed = _delay.get (n_spls, 0);
      out += delayed * pars.scaled_delay_feedback;

      if (pars.topology == t_schroeder || pars.topology == t_comb) {
        // smooth delay times calculated at control-rate
        auto smooth_iter = _stagesdl_spls.target.raw.size() / vec_size;
        for (uint j = 0; j < smooth_iter; ++j) {
          // HACKish encapsulation violation: this one isn't assigning the
          // result, but getting the -1 state from the one pole filter directly
          // (copied by value). See "static_assert" on the code below. 1-pole
          // digital filters are NEVER going to change so this is IMO
          // acceptable.
          auto in
            = make_crange (_stagesdl_spls.target.raw, vec_size, vec_size * j);
          onepole_smoother::tick (
            make_crange (_lp_smooth_coeff),
            make_crange (_stagesdl_spls.state.raw, vec_size, vec_size * j)
              .cast (float_x4 {}),
            vec_load<float_x4> (in));
        }
      }
      if (pars.topology == t_2_pole) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }
      }
      else if (pars.topology == t_3_pole) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass1p.tick_on_idx (g, out);
        }
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }
      }
      else if (pars.topology == t_1_pole) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass1p.tick_on_idx (g, out);
        }
      }
      else if (pars.topology == t_schroeder) {
        constexpr float q_scale = 0.99f / get_parameter (q_tag {}).max;

        std::array<float_x1, n_ap_channels>              acum {};
        std::array<float_x1, n_ap_channels * max_stages> to_push {};
        acum = vec1_array_wrap (vec_to_array (out));

        float_x1 gain {pars.q * q_scale};
        for (uint g = 0; g < pars.n_stages; ++g) {
          for (uint c = 0; c < n_ap_channels; ++c) {
            auto n_spls = std::clamp (
              _stagesdl_spls.state.spls[g][c],
              (float) _stagesdl.min_delay_spls(),
              (float) _stagesdl.max_delay_spls());
            uint dl_channel = g * n_ap_channels + c;
            auto yn         = _stagesdl.get (n_spls, dl_channel);
            auto r          = allpass_fn::tick<float_x1> (acum[c], yn, gain);
            acum[c]         = r.out;
            to_push[dl_channel] = r.to_push;
          }
        }
        _stagesdl.push (to_push);
        out = vec_from_array (vec1_array_unwrap (acum));
      }
      else if (pars.topology == t_comb) {
        constexpr float q_scale = 0.99f / get_parameter (q_tag {}).max;
        std::array<float, n_ap_channels> acum {};
        float                            spice_mix {pars.q * q_scale};
        spice_mix *= spice_mix;

        auto     spice = _allpass1p.tick_on_idx (0, out);
        float_x4 spice_mul {0.f, 0.f, spice_mix, spice_mix};
        out          = spice * spice_mul + out * (1.f - spice_mul);
        auto to_push = vec1_array_wrap (vec_to_array (out));
        _ff_comb.push (to_push);

        auto pan_start = make_array (
          0.999f, 0.95f, 0.98f, -0.942f, 0.999f, -0.921f, 0.98f, 0.952f);

        for (uint c = 0; c < n_ap_channels; ++c) {
          auto pans = pan_start;
          for (uint g = 0; g < pars.n_stages; ++g) {
            auto n_spls = std::clamp (
              _stagesdl_spls.state.spls[g][c],
              (float) _ff_comb.min_delay_spls(),
              (float) _ff_comb.max_delay_spls());
            auto  spl  = _ff_comb.get (n_spls, c);
            auto& panc = pans[(c + g) % pans.size()];
            acum[c] += spl[0] * panc;
            panc *= -panc;
          }
        }
        out = vec_from_array (acum);
        out /= (float) pars.n_stages;
      }
      else if (pars.topology == t_hybrid_1) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }

        auto            stage   = pars.n_stages - 1;
        constexpr float q_scale = 0.19f / get_parameter (q_tag {}).max;

        std::array<float_x1, n_ap_channels>              acum {};
        std::array<float_x1, n_ap_channels * max_stages> to_push {};
        acum = vec1_array_wrap (vec_to_array (out));

        float_x1 gain {0.8f + pars.q * q_scale};
        for (uint c = 0; c < n_ap_channels; ++c) {
          auto n_spls = std::clamp (
            _stagesdl_spls.state.spls[stage][c],
            (float) _stagesdl.min_delay_spls(),
            (float) _stagesdl.max_delay_spls());
          uint dl_channel     = stage * n_ap_channels + c;
          auto yn             = _stagesdl.get (n_spls, dl_channel);
          auto r              = allpass_fn::tick<float_x1> (acum[c], yn, gain);
          acum[c]             = r.out;
          to_push[dl_channel] = r.to_push;
        }
        _stagesdl.push (to_push);
        out = vec_from_array (vec1_array_unwrap (acum));
      }
      else if (pars.topology == t_hybrid_2) {
        std::array<float, n_ap_channels> acum {};
        auto to_push = vec1_array_wrap (vec_to_array (out));
        _ff_comb.push (to_push);

        for (uint c = 0; c < n_ap_channels; ++c) {
          auto n_spls = std::clamp (
            _stagesdl_spls.state.spls[0][c],
            (float) _ff_comb.min_delay_spls(),
            (float) _ff_comb.max_delay_spls());
          auto spl = _ff_comb.get (n_spls, c);
          acum[c]  = spl[0];
        }

        constexpr float q_scale = 0.3f / get_parameter (q_tag {}).max;
        float           gain {0.3f + pars.q * q_scale};
        out += vec_from_array (acum) * gain;
        out /= 1.f + gain;
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }
      }
      out         = _dc_blocker.tick (out);
      auto del_fb = _feedback_delay_shelf.tick_on_idx (k_shelf_lo, out);
      del_fb      = _feedback_delay_shelf.tick_on_idx (k_shelf_hi, out);
      _delay.push (make_crange (del_fb));

      double_x2 outx2    = {(double) out[0], (double) out[1]};
      double_x2 parallel = {(double) out[2], (double) out[3]};

      outx2 *= 0.5f + (0.5f - abs (pars.parallel_mix));
      outx2 += pars.parallel_mix * parallel;

      auto main_fb = outx2;
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_lo, main_fb);
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_hi, main_fb);
      main_fb = main_fb / vec_sqrt (main_fb * main_fb * pars.feedback_drv + 1.);
      _feedback_samples = main_fb;

      // feedforward outside the feedback loop
      if (pars.topology > t_2_pole_legacy && pars.topology != t_comb) {
        outx2 += feedforward * pars.scaled_feedback;
        outx2 *= 0.5f;
      }
      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];
    }
#endif
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    stages_tag,
    stages_mode_tag,
    lfo_rate_tag,
    lfo_depth_tag,
    lfo_start_phase_tag,
    lfo_stereo_tag,
    lfo_wave_tag,
    low_freq_tag,
    high_freq_tag,
    feedback_tag,
    q_tag,
    parallel_mix_tag,
    delay_feedback_tag,
    feedback_hp_tag,
    feedback_drv_tag,
    delay_time_tag,
    delay_lfo_tag>;
  //----------------------------------------------------------------------------
private:
  static constexpr uint  max_stages    = 16;
  static constexpr float min_ap_freq   = 20.f;
  static constexpr uint  n_channels    = 2;
  static constexpr uint  n_ap_channels = 4;
  //----------------------------------------------------------------------------
  template <class T>
  void process_legacy (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (_params.unsmoothed.topology < t_first_non_legacy);
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {
      // parameter smoothing (could be done on bigger blocks).
      all_parameters pars;
      smooth_parameters (pars);

      // control rate/mdoulation block -----------------------------------------
      if ((_n_processed_samples & _control_rate_mask) == 0) {
        run_control_rate (pars);
      }
      // sample-wise processing block (unbreakable, it has 1 sample delays) ----
      float_x4 fb = {
        (float) _feedback_samples[0],
        (float) _feedback_samples[1],
        (float) _feedback_samples[0],
        (float) _feedback_samples[1]};
      // regular processing
      float_x4 out {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      out += fb * pars.scaled_feedback;
      double_x2 feedforward {out[0], out[1]};

      auto n_spls = pars.delay_spls;
      n_spls *= exp (ln_max_delay_factor * pars.lfo_last * pars.delay_lfo);
      auto delayed = _delay.get (n_spls, 0);
      out += delayed * pars.scaled_delay_feedback;

      if (pars.topology == t_2_pole_legacy) {
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass2p.tick_on_idx (g, out);
        }
      }
      else {
        assert (pars.topology == t_1_pole_legacy);
        for (uint g = 0; g < pars.n_stages; ++g) {
          out = _allpass1p.tick_on_idx (g, out);
        }
      }
      out         = _dc_blocker.tick (out);
      auto del_fb = _feedback_delay_shelf.tick_on_idx (k_shelf_lo, out);
      del_fb      = _feedback_delay_shelf.tick_on_idx (k_shelf_hi, out);
      _delay.push (make_crange (del_fb));

      double_x2 outx2    = {(double) out[0], (double) out[1]};
      double_x2 parallel = {(double) out[2], (double) out[3]};

      outx2 *= 0.5f + (0.5f - abs (pars.parallel_mix));
      outx2 += pars.parallel_mix * parallel;

      auto main_fb = outx2;
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_lo, main_fb);
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_hi, main_fb);
      main_fb = main_fb / vec_sqrt (main_fb * main_fb * pars.feedback_drv + 1.);
      _feedback_samples = main_fb;

      // feedforward outside the feedback loop
      if (pars.topology > t_2_pole_legacy && pars.topology != t_comb) {
        outx2 += feedforward * pars.scaled_feedback;
        outx2 *= 0.5f;
      }
      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_zdf (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (_params.unsmoothed.topology >= t_first_non_legacy);
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {
      // parameter smoothing (could be done on bigger blocks).
      all_parameters pars;
      smooth_parameters (pars);

      // control rate/mdoulation block -----------------------------------------
      if ((_n_processed_samples & _control_rate_mask) == 0) {
        run_control_rate (pars);
      }

      float_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};

      float dry_gain = 1.f - pars.mix;
      float wet_gain = pars.mix;

      if (pars.topology < t_first_non_legacy_non_allpass) {
        // Run allpass cascades with own feedback loop and saturation
        bool has_twopoles = pars.topology < t_first_non_legacy_non_2_pole;
        bool has_onepoles
          = pars.topology == t_1_pole || pars.topology == t_3_pole;

        constexpr uint n_fb_filters = 2;
        // G and S zdf coeffs
        std::array<float_x4, (max_stages * 2 + n_fb_filters) * 2> G_S;
        auto gs_it = make_crange (G_S);
        if (has_onepoles) {
          for (uint s = 0; s < pars.n_stages; ++s) {
            _allpass1p_zdf.tick_on_idx (
              s, gs_it.cut_head (2), zdf::gs_coeffs_tag {});
          }
        }
        if (has_twopoles) {
          for (uint s = 0; s < pars.n_stages; ++s) {
            _zdf_ap.tick_on_idx (
              s, gs_it.cut_head (2), zdf::gs_coeffs_tag {});
          }
        }
        _zdf_lshelf.tick (gs_it.cut_head (2), zdf::gs_coeffs_tag {});
        _zdf_hshelf.tick (gs_it.cut_head (2), zdf::gs_coeffs_tag {});

        auto n_stages
          = pars.n_stages * ((uint) has_onepoles + (uint) has_twopoles);
        n_stages += n_fb_filters;
        auto  gs_final = make_crange (G_S.data(), n_stages * 2);
        auto  resp     = zdf::get_response_series<float_x4> (gs_final);
        float hardness = pars.feedback_drv * pars.feedback_drv * 15.f;
        wet            = _zero_feedback.solve (
          resp, vec_set<4> (pars.feedback), vec_set<4> (hardness), wet);

        if (has_onepoles) {
          for (uint s = 0; s < pars.n_stages; ++s) {
            wet = _allpass1p_zdf.tick_on_idx (s, wet);
          }
        }
        if (has_twopoles) {
          for (uint s = 0; s < pars.n_stages; ++s) {
            wet = _zdf_ap.tick_on_idx (s, wet);
          }
        }
        // just running the shelves to update the states.
        _zdf_hshelf.tick (_zdf_lshelf.tick (wet));
        // limit with x=inf of sqrt sigmoid
        float lim = 1.f / sqrt (hardness != 0.f ? hardness : 1e-30f);
        // gain of the zdf part
        float zdf_gain = 1.f + std::min (lim, abs (pars.feedback));
        // gain compensate so the other feedback loops in parallel don't blow
        // up
        wet /= zdf_gain;
      }

      double_x2 wetx2    = {(double) wet[0], (double) wet[1]};
      double_x2 parallel = {(double) wet[2], (double) wet[3]};
      double_x2 dry      = {(double) ins[0][i], (double) ins[1][i]};
      wetx2 *= 1.f - abs (pars.parallel_mix);
      wetx2 += pars.parallel_mix * parallel;
      auto outx2 = wetx2 * wet_gain + dry * dry_gain;

      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];

#if 0
      // 1 sample delay for the delay-based parts ------------------------------
      float_x4 fb = {
        (float) _feedback_samples[0],
        (float) _feedback_samples[1],
        (float) _feedback_samples[0],
        (float) _feedback_samples[1]};

      wet += fb * pars.scaled_feedback;
      dry_gain_comp *= 1.f + pars.scaled_feedback;

      auto n_spls = pars.delay_spls;
      n_spls *= exp (ln_max_delay_factor * pars.lfo_last * pars.delay_lfo);
      auto delayed = _delay.get (n_spls, 0);
      wet += delayed * pars.scaled_delay_feedback;

      if (pars.topology >= t_first_non_legacy_non_allpass) {
        // smooth delay times calculated at control-rate
        auto smooth_iter = _stagesdl_spls.target.raw.size() / vec_size;
        for (uint j = 0; j < smooth_iter; ++j) {
          // HACKish encapsulation violation: this one isn't assigning the
          // result, but getting the -1 state from the one pole filter directly
          // (copied by value). See "static_assert" on the code below. 1-pole
          // digital filters are NEVER going to change so this is IMO
          // acceptable.
          auto in
            = make_crange (_stagesdl_spls.target.raw, vec_size, vec_size * j);
          onepole_smoother::tick (
            make_crange (_lp_smooth_coeff),
            make_crange (_stagesdl_spls.state.raw, vec_size, vec_size * j)
              .cast (float_x4 {}),
            vec_load<float_x4> (in));
        }
      }
      else if (pars.topology == t_schroeder) {
        constexpr float q_scale = 0.99f / get_parameter (q_tag {}).max;

        std::array<float_x1, n_ap_channels>              acum {};
        std::array<float_x1, n_ap_channels * max_stages> to_push {};
        acum = vec1_array_wrap (vec_to_array (wet));

        float_x1 gain {pars.q * q_scale};
        for (uint g = 0; g < pars.n_stages; ++g) {
          for (uint c = 0; c < n_ap_channels; ++c) {
            auto n_spls = std::clamp (
              _stagesdl_spls.state.spls[g][c],
              (float) _stagesdl.min_delay_spls(),
              (float) _stagesdl.max_delay_spls());
            uint dl_channel = g * n_ap_channels + c;
            auto yn         = _stagesdl.get (n_spls, dl_channel);
            auto r          = allpass_fn::tick<float_x1> (acum[c], yn, gain);
            acum[c]         = r.out;
            to_push[dl_channel] = r.to_push;
          }
        }
        _stagesdl.push (to_push);
        wet = vec_from_array (vec1_array_unwrap (acum));
      }
      else if (pars.topology == t_comb) {
        constexpr float q_scale = 0.99f / get_parameter (q_tag {}).max;
        std::array<float, n_ap_channels> acum {};
        float                            spice_mix {pars.q * q_scale};
        spice_mix *= spice_mix;

        auto     spice = _allpass1p.tick_on_idx (0, wet);
        float_x4 spice_mul {0.f, 0.f, spice_mix, spice_mix};
        wet          = spice * spice_mul + wet * (1.f - spice_mul);
        auto to_push = vec1_array_wrap (vec_to_array (wet));
        _ff_comb.push (to_push);

        auto pan_start = make_array (
          0.999f, 0.95f, 0.98f, -0.942f, 0.999f, -0.921f, 0.98f, 0.952f);

        for (uint c = 0; c < n_ap_channels; ++c) {
          auto pans = pan_start;
          for (uint g = 0; g < pars.n_stages; ++g) {
            auto n_spls = std::clamp (
              _stagesdl_spls.state.spls[g][c],
              (float) _ff_comb.min_delay_spls(),
              (float) _ff_comb.max_delay_spls());
            auto  spl  = _ff_comb.get (n_spls, c);
            auto& panc = pans[(c + g) % pans.size()];
            acum[c] += spl[0] * panc;
            panc *= -panc;
          }
        }
        wet = vec_from_array (acum);
        wet /= (float) pars.n_stages;
      }
      else if (pars.topology == t_hybrid_1) {
        auto            stage   = pars.n_stages - 1;
        constexpr float q_scale = 0.19f / get_parameter (q_tag {}).max;

        std::array<float_x1, n_ap_channels>              acum {};
        std::array<float_x1, n_ap_channels * max_stages> to_push {};
        acum = vec1_array_wrap (vec_to_array (wet));

        float_x1 gain {0.8f + pars.q * q_scale};
        for (uint c = 0; c < n_ap_channels; ++c) {
          auto n_spls = std::clamp (
            _stagesdl_spls.state.spls[stage][c],
            (float) _stagesdl.min_delay_spls(),
            (float) _stagesdl.max_delay_spls());
          uint dl_channel     = stage * n_ap_channels + c;
          auto yn             = _stagesdl.get (n_spls, dl_channel);
          auto r              = allpass_fn::tick<float_x1> (acum[c], yn, gain);
          acum[c]             = r.out;
          to_push[dl_channel] = r.to_push;
        }
        _stagesdl.push (to_push);
        wet = vec_from_array (vec1_array_unwrap (acum));
      }
      else if (pars.topology == t_hybrid_2) {
        std::array<float, n_ap_channels> acum {};
        auto to_push = vec1_array_wrap (vec_to_array (wet));
        _ff_comb.push (to_push);

        for (uint c = 0; c < n_ap_channels; ++c) {
          auto n_spls = std::clamp (
            _stagesdl_spls.state.spls[0][c],
            (float) _ff_comb.min_delay_spls(),
            (float) _ff_comb.max_delay_spls());
          auto spl = _ff_comb.get (n_spls, c);
          acum[c]  = spl[0];
        }

        constexpr float q_scale = 0.3f / get_parameter (q_tag {}).max;
        float           gain {0.3f + pars.q * q_scale};
        wet += vec_from_array (acum) * gain;
        wet /= 1.f + gain;
      }
      wet         = _dc_blocker.tick (wet);
      auto del_fb = _feedback_delay_shelf.tick_on_idx (k_shelf_lo, wet);
      del_fb      = _feedback_delay_shelf.tick_on_idx (k_shelf_hi, wet);
      _delay.push (make_crange (del_fb));

      double_x2 outx2    = {(double) wet[0], (double) wet[1]};
      double_x2 parallel = {(double) wet[2], (double) wet[3]};

      outx2 *= 0.5f + (0.5f - abs (pars.parallel_mix));
      outx2 += pars.parallel_mix * parallel;

      auto main_fb = outx2;
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_lo, main_fb);
      main_fb      = _feedback_shelf.tick_on_idx (k_shelf_hi, main_fb);
      main_fb = main_fb / vec_sqrt (main_fb * main_fb * pars.feedback_drv + 1.);
      _feedback_samples = main_fb;

      // TODO: MIX
      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];
#endif
    }
  }
  //----------------------------------------------------------------------------
  // As feedbacks are smoothed, this might need to be done sample-wise.
  void scale_feedbacks (
    float& g1,
    float& g2,
    float  unscaled_g1,
    float  unscaled_g2)
  {
    // keep the sum of both feebacks below unity
    g1            = unscaled_g1;
    g2            = unscaled_g2;
    float fb_gain = abs (g1) + abs (g2);
    float att     = (fb_gain <= max_feedback) ? 1.f : (max_feedback / fb_gain);
    g1 *= att;
    g2 *= att;
  }
  //----------------------------------------------------------------------------
  void reset_memory_parts()
  {
    uint headroom  = _ff_comb.min_size_spls();
    auto comb_size = pow2_round_ceil (
      (uint) ceil (freq_to_delay_spls (min_ap_freq)) + headroom);
    comb_size = _ff_comb.n_required_elems (comb_size, n_ap_channels);
    comb_size = div_ceil<unsigned> (
      comb_size, sizeof _delay_mem[0] / sizeof _ff_comb.get (0.f, 0));

    auto ap_size
      = pow2_round_ceil ((uint) ceil (freq_to_delay_spls (min_ap_freq)));
    ap_size = _stagesdl.n_required_elems (ap_size, n_ap_channels * max_stages);
    ap_size = div_ceil<unsigned> (
      ap_size, sizeof _delay_mem[0] / sizeof _stagesdl.get (0.f, 0));

    auto delay_size
      = pow2_round_ceil ((uint) (_srate * max_delay_sec * max_delay_factor));
    delay_size = _delay.n_required_elems (delay_size, 1);
    delay_size = div_ceil<unsigned> (
      delay_size, sizeof _delay_mem[0] / sizeof _delay.get (0.f, 0));

    _delay_mem.clear();
    _delay_mem.resize (comb_size + delay_size + ap_size);
    auto mem = make_crange (_delay_mem);

    _ff_comb.reset (mem.cut_head (comb_size).cast<float_x1>(), n_ap_channels);
    _stagesdl.reset (
      mem.cut_head (ap_size).cast<float_x1>(), n_ap_channels * max_stages);
    _delay.reset (mem.cut_head (delay_size), 1);
  }
  //----------------------------------------------------------------------------
  void refresh_lfo_hz()
  {
    float hz = 0.f;
    if (_params.unsmoothed.lfo_eights != 0.f) {
      hz = (8. * _plugcontext->get_play_state().bpm)
        / (60. * _params.unsmoothed.lfo_eights);
    }
    hz += _params.unsmoothed.lfo_hz_user;
    _params_smooth.target().lfo_hz_final = hz;
  }
  //----------------------------------------------------------------------------
  template <class T>
  T freq_to_delay_spls (T freq)
  {
    if constexpr (is_vec_v<T>) {
      using VT = vec_value_type_t<T>;
      return (VT) _srate / ((VT) 2 * freq);
    }
    else {
      return (T) _srate / ((T) 2 * freq);
    }
  }
  //----------------------------------------------------------------------------
  void run_control_rate (all_parameters const& pars)
  {
    // this is very naive, note that we could e.g. be interpolating SVF
    // coefficients internally.
    constexpr uint vec_size = vec_traits<float_x4>().size;

    _lfos.set_freq (vec_set<2> (pars.lfo_hz_final), _lfo_t_spl);
    auto stereo_ph
      = phase<1> {vec_set<1> (pars.lfo_stereo), phase<1>::degrees {}}.get_raw (
        0);
    if (pars.lfo_hz_final != 0.f) {
      // updating stereo phase diff.
      auto ph = _lfos.get_phase();
      ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (ph);
    }
    else {
      auto start_ph = phase<n_channels> {
        vec_set<n_channels> (pars.lfo_start_phase),
        phase<n_channels>::degrees {}};
      start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
      _lfos.set_phase (start_ph);
    }
    vec<float, 2> lfov;
    switch (pars.lfo_wave) {
    case 0:
      lfov = _lfos.tick_sine();
      break;
    case 1:
      lfov = _lfos.tick_triangle();
      break;
    case 2:
      lfov = _lfos.tick_sample_hold();
      break;
    case 3:
      lfov = _lfos.tick_filt_sample_and_hold();
      break;
    case 4:
      lfov = _lfos.tick_trapezoid (vec_set<2> (0.3f));
      break;
    case 5:
      lfov = _lfos.tick_square();
      break;
    case 6:
      lfov = _lfos.tick_saw();
      break;
    }
    _params_smooth.target().lfo_last = lfov[0];
    auto                  interp_stages  = (double) (pars.n_stages * 2) - 1;
    std::array<double, 2> fconstant;

    auto freq_hi = std::max (pars.freq_hi, pars.freq_lo);
    auto freq_lo = std::min (pars.freq_hi, pars.freq_lo);

    switch (pars.mode) {
    case m_exp:
    case m_exp_spread:
    case m_exp_alternate:
      fconstant[0] = 1.;
      fconstant[1] = pow (freq_hi / freq_lo, 1. / interp_stages);
      break;
    case m_exp_stereo:
      interp_stages *= 2;
      fconstant[0] = pow (freq_hi / freq_lo, 1. / interp_stages);
      fconstant[1] = fconstant[0];
      break;
    case m_lin:
    case m_lin_spread:
    case m_lin_alternate:
      fconstant[0] = 0.f;
      fconstant[1] = (freq_hi - freq_lo) / interp_stages;
      break;
    case m_lin_stereo:
      interp_stages *= 2;
      fconstant[0] = (freq_hi - freq_lo) / interp_stages;
      fconstant[1] = fconstant[0];
      break;
    default:
      fconstant[0] = 0.;
      fconstant[1] = 0.;
      break;
    }
    auto     f         = freq_lo;
    auto     lfo_depth = pars.lfo_depth;
    float_x4 q_fact {lfov[0], lfov[1], lfov[1], lfov[0]};
    q_fact = vec_exp (q_fact * -lfo_depth * 0.23f);
    if (pars.topology == t_schroeder || pars.topology == t_comb) {
      // favor the lower range
      lfo_depth *= lfo_depth;
    }
    for (uint s = 0; s < pars.n_stages; ++s) {
      float_x4 freqs {};
      float_x4 qs = vec_set<float_x4> (pars.q);
      if (pars.topology >= t_2_pole_legacy) {
        qs *= q_fact;
      }
      // fill freqs
      for (uint ap = 0; ap < (vec_size / 2); ++ap) {
        float max_depth = pars.lin_lfo_mod ? 0.6 : 1.7;
        float k;
        uint  stage = (s * 4) + ap;

        switch (pars.mode) {
        case m_exp:
        case m_lin:
        case m_exp_stereo:
        case m_lin_stereo:
          k = max_depth;
          break;
        case m_exp_spread:
        case m_lin_spread:
          k = (stage < (pars.n_stages * 2)) ? max_depth : -max_depth;
          break;
        case m_exp_alternate:
        case m_lin_alternate:
          k = ((stage / 2) & 1) ? max_depth : -max_depth;
          break;
        default:
          break;
        }

        std::array<decltype (f), 2> fs;
        for (uint iv = 0; iv < fs.size(); ++iv) {
          fs[iv] = f;
          if (mode_is_exponential (pars.mode)) {
            f *= fconstant[iv];
          }
          else {
            f += fconstant[iv];
          }
        }

        if (pars.lin_lfo_mod) {
          freqs[ap * 2]     = fs[0] + (lfov[0] * lfo_depth * fs[0] * k);
          freqs[ap * 2 + 1] = fs[1] + (lfov[1] * lfo_depth * fs[1] * k);
        }
        else {
          freqs[ap * 2]     = fs[0] * exp (lfov[0] * lfo_depth * k);
          freqs[ap * 2 + 1] = fs[1] * exp (lfov[1] * lfo_depth * k);
        }
      }
      freqs
        = vec_clamp (freqs, min_ap_freq, std::min (_srate * 0.4999f, 22000.f));
      switch (pars.topology) {
      case t_1_pole_legacy:
        _allpass1p.reset_coeffs_on_idx (s, freqs, _t_spl);
        break;
      case t_2_pole_legacy:
        _allpass2p.reset_coeffs_on_idx (s, freqs, qs, _t_spl);
        break;
      case t_3_pole:
        _allpass1p_zdf.reset_coeffs_on_idx (s, freqs, _t_spl);
        _zdf_ap.reset_coeffs_on_idx (s, freqs, qs, _t_spl);
        break;
      case t_hybrid_1:
        _stagesdl_spls.target.spls[s]
          = vec_to_array (freq_to_delay_spls (freqs * 2));
        _zdf_ap.reset_coeffs_on_idx (
          s, freqs, 0.05f + qs * 0.3f, _t_spl);
        break;
      case t_hybrid_2:
        _stagesdl_spls.target.spls[s]
          = vec_to_array (freq_to_delay_spls (freqs));
        _zdf_ap.reset_coeffs_on_idx (
          s, freqs, 0.05f + qs * 0.5f, _t_spl);
        break;
      case t_2_pole:
        _zdf_ap.reset_coeffs_on_idx (s, freqs, qs, _t_spl);
        break;
      case t_1_pole:
        _allpass1p_zdf.reset_coeffs_on_idx (s, freqs, _t_spl);
        break;
      case t_schroeder:
        _stagesdl_spls.target.spls[s]
          = vec_to_array (freq_to_delay_spls (freqs));
        break;
      case t_comb:
        _stagesdl_spls.target.spls[s]
          = vec_to_array (freq_to_delay_spls (freqs));
        if (s == 0) {
          _allpass1p_zdf.reset_coeffs_on_idx (0, freqs, _t_spl);
        }
        break;
      }
      if (s == 0) {
        constexpr float cut_ratio = 0.3f;
        constexpr float gain_db   = -20.f;

        auto cutfreq = double_x2 {freqs[0], freqs[1]};
        cutfreq -= pars.feedback_hp * cutfreq * cut_ratio;
        _feedback_shelf.reset_coeffs_on_idx (
          k_shelf_lo,
          vec_max (210., cutfreq),
          vec_set<double_x2> (0.35),
          vec_set<double_x2> (pars.feedback_hp * gain_db),
          _t_spl,
          lowshelf_tag {});
        auto del_cutfreq = freqs;
        del_cutfreq -= pars.feedback_hp * freqs * cut_ratio;
        _feedback_delay_shelf.reset_coeffs_on_idx (
          k_shelf_lo,
          vec_max (210.f, del_cutfreq),
          vec_set<float_x4> (0.35),
          vec_set<float_x4> (pars.feedback_hp * gain_db),
          _t_spl,
          lowshelf_tag {});
        _zdf_lshelf.reset_coeffs (
          vec_max (210., del_cutfreq),
          vec_set<float_x4> (pars.feedback_hp * gain_db),
          _t_spl);
      }
      if (s == (pars.n_stages - 1)) {
        constexpr float cut_ratio = 0.3f;
        constexpr float gain_db   = -6.f;

        auto cutfreq = double_x2 {freqs[0], freqs[1]};
        cutfreq += pars.feedback_lp * cutfreq * cut_ratio;
        _feedback_shelf.reset_coeffs_on_idx (
          k_shelf_hi,
          vec_max (4000., cutfreq),
          vec_set<double_x2> (0.35),
          vec_set<double_x2> (pars.feedback_lp * gain_db),
          _t_spl,
          highshelf_tag {});

        auto del_cutfreq = freqs;
        del_cutfreq -= pars.feedback_hp * freqs * cut_ratio;
        _feedback_delay_shelf.reset_coeffs_on_idx (
          k_shelf_hi,
          vec_max (4000.f, del_cutfreq),
          vec_set<float_x4> (0.35),
          vec_set<float_x4> (pars.feedback_lp * gain_db),
          _t_spl,
          highshelf_tag {});
        _zdf_hshelf.reset_coeffs (
          vec_max (4000.f, del_cutfreq),
          vec_set<float_x4> (pars.feedback_lp * gain_db),
          _t_spl);
      }
    }
  }
  struct unsmoothed_parameters {
    float lfo_eights;
    float lfo_hz_user;
    float delay_feedback;
    uint  lfo_wave;
    uint  n_stages;
    uint  topology;
    uint  mode;
    bool  lin_lfo_mod;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float freq_lo;
    float freq_hi;
    float q;
    float lfo_hz_final;

    float lfo_stereo;
    float lfo_depth;
    float lfo_start_phase;
    float feedback;

    float delay_spls;
    float scaled_delay_feedback;
    float delay_lfo;
    float lfo_last;

    float parallel_mix;
    float feedback_hp;
    float feedback_lp;
    float feedback_drv;

    float scaled_feedback;
    float mix;
  };
  //----------------------------------------------------------------------------
  struct all_parameters : public unsmoothed_parameters,
                          public smoothed_parameters {};
  //----------------------------------------------------------------------------
  union smoothed_parameter_union {
    smoothed_parameters value;
    simd_array<float, sizeof (smoothed_parameters) / sizeof (float), sse_bytes>
      arr;
  };
  //----------------------------------------------------------------------------
  struct parameter_values {
    alignas (sse_bytes) smoothed_parameter_union smooth_target;
    alignas (sse_bytes) smoothed_parameter_union smooth_state;
    unsmoothed_parameters unsmoothed;
  };
  //----------------------------------------------------------------------------
  //
  union smoothed_stages_delay_samples {
    using array_type = array2d<float, n_ap_channels, max_stages>;
    array_type                                                         spls;
    simd_array<float, sizeof (array_type) / sizeof (float), sse_bytes> raw;
  };
  //----------------------------------------------------------------------------
  struct stages_delay_samples {
    alignas (sse_bytes) smoothed_stages_delay_samples target;
    alignas (sse_bytes) smoothed_stages_delay_samples state;
  };
  //----------------------------------------------------------------------------
  double_x2        _feedback_samples;
  parameter_values _params;

  std::vector<float_x4, overaligned_allocator<float_x4, 128>> _delay_mem;

  // TODO: use a variant?
  part_class_array<andy::svf_allpass, float_x4, max_stages> _allpass2p;
  part_class_array<onepole_allpass, float_x4, max_stages>   _allpass1p;

  zdf::feedback<
    float_x4,
    zdf::lin_pre_fb_node_nonlin_after_tag,
    zdf::lin_mystran_2_tag,
    sigmoid::rsqrt<true>>
    _zero_feedback;
  part_class_array<andy::svf_multimode_zdf<allpass_tag>, float_x4, max_stages>
    _zdf_ap;
  part_class_array<onepole_zdf<allpass_tag>, float_x4, max_stages>
                                                                 _allpass1p_zdf;
  part_class_array<onepole_zdf<lowshelf_naive_tag>, float_x4, 1> _zdf_lshelf;
  part_class_array<onepole_zdf<highshelf_naive_tag>, float_x4, 1> _zdf_hshelf;

  interpolated_delay_line<float_x1, linear_interp, true, false>       _stagesdl;
  interpolated_delay_line<float_x1, sinc_interp<8, 64>, false, false> _ff_comb;
  enum { k_shelf_lo, k_shelf_hi };
  part_class_array<andy::svf, double_x2, 2> _feedback_shelf;
  part_class_array<andy::svf, float_x4, 2>  _feedback_delay_shelf;

  stages_delay_samples                           _stagesdl_spls {};
  modulable_thiran1_delay_line<float_x4, 1>      _delay {};
  part_class_array<mystran_dc_blocker, float_x4> _dc_blocker;

  lfo<n_channels> _lfos; // 0 = L, 1 = R
  uint            _n_processed_samples;
  uint            _control_rate_mask;
  float           _lfo_t_spl;
  float           _t_spl;
  float           _srate;
  float           _lp_smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
#endif
} // namespace artv
