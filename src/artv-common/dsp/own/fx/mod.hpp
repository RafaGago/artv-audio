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
  struct stages_tag {};
  void set (stages_tag, float v) { _param_smooth.target().stages = v * 0.01f; }

  static constexpr auto get_parameter (stages_tag)
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
    _param_smooth.target().feedback_drive = v * 1400.f;
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
    if (v == _param.mode) {
      return;
    }
    _param.mode = v;
    switch (v) {
    case mode::zdf:
      break;
    case mode::schroeder:
      using TC = decltype (_scho)::value_type;
      _scho.reset (make_crange (_mem).cast<TC>(), max_scho_stages);
      //_scho.set_resync_delta (10.0);
      break;
    case mode::comb:
      using TF = decltype (_comb)::value_type;
      _comb.reset (make_crange (_mem).cast<TF>(), max_comb_stages);
      // initialize windowed sinc table
      _comb.reset_interpolator (0, false, 0.45f, 140.f);
      break;
    default:
      assert (false);
      break;
    }
  }

  struct mode {
    enum { zdf, schroeder, comb };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array ("Phaser (ZDF)", "Phaser (Schroeder)", "Flanger/Chorus"),
      16);
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
    _del_spls.reset (_t_spl, 10.f);

    // ensure that all parameter's new values will be different
    memset (&_param, -1, sizeof _param);
    memset (&_param_smooth.target(), -1, sizeof _param_smooth.target());

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _param_smooth.set_all_from_target();

    reset_mem();
    using phase = decltype (_rnd_lfo)::phase_type;
    _rnd_lfo.set_phase (phase {phase::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    _rnd_lfo.set_freq (vec_set<4> (0.3f), _t_spl);
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
    case mode::comb:
      tick_comb<T> (outs, ins, samples);
      break;
    default:
      assert (false);
      break;
    }
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    stages_tag,
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
  static constexpr uint max_comb_delay_ms = 30;
  static constexpr uint max_comb_stages   = 8;
  static constexpr uint n_channels        = 2;
  static constexpr uint n_ap_channels     = 4;
  //----------------------------------------------------------------------------
  struct all_parameters;
  using sweep_lfo_type  = lfo<n_channels>;
  using sweep_lfo_value = sweep_lfo_type::value_type;
  //----------------------------------------------------------------------------
  sweep_lfo_value run_lfo (all_parameters const& pars)
  {
    float hz;
    switch (pars.lfo_time_base) {
    case lfo_time_base::free: {
      float lfo_rate = pars.lfo_rate * 0.01f;
      lfo_rate       = 1.f - lfo_rate;
      lfo_rate *= lfo_rate * lfo_rate;
      hz = lfo_rate * 16.f;
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

    _sweep_lfo.set_freq (vec_set<2> (hz), _lfo_t_spl);
    auto stereo_ph
      = phase<1> {phase_tag::degrees {}, pars.lfo_stereo}.get_raw (0);
    if (!pars.lfo_off) {
      // updating stereo phase diff.
      auto ph = _sweep_lfo.get_phase();
      ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
      _sweep_lfo.set_phase (ph);
    }
    else {
      auto start_ph = phase<n_channels> {phase_tag::degrees {}, 0.f, 0.f};
      start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
      _sweep_lfo.set_phase (start_ph);
    }

    sweep_lfo_value lfov;
    switch (pars.lfo_wave) {
    case lfo_wave::sine:
      lfov = _sweep_lfo.tick_sine();
      break;
    case lfo_wave::triangle:
      lfov = _sweep_lfo.tick_triangle();
      break;
    case lfo_wave::sample_hold:
      lfov = _sweep_lfo.tick_sample_hold();
      break;
    case lfo_wave::noise:
      lfov = _sweep_lfo.tick_filt_sample_and_hold();
      break;
    case lfo_wave::trapezoid:
      lfov = _sweep_lfo.tick_trapezoid (vec_set<2> (0.3f));
      break;
    case lfo_wave::square:
      lfov = _sweep_lfo.tick_square();
      break;
    case lfo_wave::saw:
      lfov = _sweep_lfo.tick_saw();
      break;
    }
    return lfov;
  }
  //----------------------------------------------------------------------------
  void run_phaser_mod (
    all_parameters const& pars,
    uint                  n_stages,
    float                 detune_ratio,
    sweep_lfo_value       lfo)
  {
    float modlim = std::min (pars.center, 0.5f);
    modlim       = std::min (1.f - pars.center, modlim);

    sweep_lfo_value mod    = lfo * pars.lfo_depth * modlim;
    sweep_lfo_value center = pars.center + mod;
    // detune by a ramp
    sweep_lfo_value detuned = (1.f - abs (pars.detune * detune_ratio)) * center;
    sweep_lfo_value diff    = (center - detuned) / (float) (n_stages * 2);
    sweep_lfo_value start   = (pars.detune > 0.f) ? detuned : center;
    diff                    = (pars.detune > 0.f) ? diff : -diff;
    float_x4 curr           = vec_cat (start, start + diff);
    float_x4 add            = vec_cat (diff, diff) * 2.f;

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
        _n_stages = 1 + (uint) (pars.stages * (max_phaser_stages - 1));
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

      auto out   = mix (pars, wet, pars.b, 1.f, ins[0][i], ins[1][i]);
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
        _n_stages = 1 + (uint) (pars.stages * (max_scho_stages - 1));

        pars.center = 1.f - pars.center; // reverse range lf to hf
        pars.center /= _n_stages;
        run_phaser_mod (pars, _n_stages, 0.25f, run_lfo (pars));

        for (uint s = 0; s < _n_stages; ++s) {
          constexpr float sec_factor = (max_scho_delay_ms * 0.001f);
          float_x4        v          = vec_from_array (_mod[s]);
          float           t          = v[0];
          t                          = _srate * (t * t * sec_factor);
          _del_spls.target()[s]      = t;
          constexpr float min_g      = 0.0001f;
          constexpr float g_factor   = 0.15 - min_g;
          _g.scho[s]
            = min_g + pars.a * pars.a * 0.8f + v * v * g_factor * pars.a;
        }
        for (uint s = _n_stages; s < _del_spls.target().size(); ++s) {
          _del_spls.target()[s] = _del_spls.target()[s - 1];
        }
        float_x4 fv {1460.f, 1460.f, 1670.f, 1780.f};
        auto     f = vec_from_array (_mod[_n_stages / 2]);
        f *= f;
        f *= (1.f + pars.a) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      _del_spls.tick();

      float_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      std::array<float_x4, max_scho_stages> to_push {};
      auto&                                 del_spls = _del_spls.get();
      for (uint s = 0; s < _n_stages; ++s) {
        auto n_spls = std::clamp (
          del_spls[s],
          (float) _scho.min_delay_spls(),
          (float) _scho.max_delay_spls());
        auto yn    = _scho.get (n_spls, s);
        auto r     = allpass_fn::tick<float_x4> (wet, yn, _g.scho[s]);
        wet        = r.out;
        to_push[s] = r.to_push;
      }
      _scho.push (to_push);
      auto fb_val = zdf_type::nonlin::tick (
        wet,
        vec_set<4> (pars.feedback_drive),
        vec_set<4> (pars.feedback_curve));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);

      auto out   = mix (pars, wet, pars.b, 1.f, ins[0][i], ins[1][i]);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_comb (crange<T*> outs, crange<T const*> ins, uint samples)
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
        _n_stages   = 1 + (uint) (pars.stages * (max_comb_stages - 1));
        pars.center = 1.f - pars.center; // reverse range lf to hf
        run_phaser_mod (pars, _n_stages, 0.35f, run_lfo (pars));

        for (uint s = 0; s < _n_stages; ++s) {
          constexpr float sec_factor = (max_comb_delay_ms * 0.001f);

          float t                  = _mod[s][0];
          t                        = _srate * (t * t * sec_factor);
          _del_spls.target()[s]    = t;
          constexpr float g_factor = 0.15f;
          auto            m        = _mod[s][1];
          _g.comb[s][0] = pars.a * pars.a * 0.45f + m * m * g_factor * pars.a;
          _rnd_lfo.set_freq (vec_set<4> (0.3f + pars.a * 1.f), _t_spl);
        }
        for (uint s = _n_stages; s < _del_spls.target().size(); ++s) {
          _del_spls.target()[s] = _del_spls.target()[s - 1];
        }
        float_x4 fv {1460.f, 1460.f, 1670.f, 1780.f};
        auto     f = vec_from_array (_mod[_n_stages / 2]);
        f *= f;
        f *= (1.f + pars.a) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      _del_spls.tick();

      float_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      auto& del_spls = _del_spls.get();
      auto  comb_in  = vec1_array_wrap (vec_to_array (wet));
      wet            = vec_set<4> (0.f);

      float panwidth   = pars.a;
      float pan_d      = panwidth / (_n_stages - 1);
      float pan        = (1.f - panwidth) * 0.5f;
      auto  total_gain = vec_set<4> (0.f);

      std::array<float_x1, max_comb_stages> to_push {};
      for (uint s = 0; s < _n_stages; ++s) {
        uint idx    = s % vec_traits_t<float_x4>::size;
        auto n_spls = std::clamp (
          del_spls[s],
          (float) _comb.min_delay_spls(),
          (float) _comb.max_delay_spls());
        auto yn    = _comb.get (n_spls, s);
        auto r     = allpass_fn::tick<float_x1> (comb_in[idx], yn, _g.comb[s]);
        to_push[s] = r.to_push;
        // using -x^2+2x as a bad approximation of sin(x*pi/2) pan law. It
        // matches at the 1 extreme perfectly, which is the most important for a
        // feedback loop.
        float_x4 panv {pan, 1.f - pan, pan, 1.f - pan};
        panv = -panv * panv + 2.f * panv;
        total_gain += panv;
        pan += pan_d;

        auto parallel = vec_set<4> (r.out[0]);
        parallel *= panv;
        wet += parallel;
      }
      _comb.push (to_push);

      auto fb_val = zdf_type::nonlin::tick (
        wet / total_gain, // feedback with pan gains reverted.
        vec_set<4> (pars.feedback_drive), // modulate?
        vec_set<4> (pars.feedback_curve)); // modulate?
      _1spl_fb = _fb_filters.tick_cascade (fb_val);
      // 0.9 to 1 gain modulation
      auto lfov
        = (_rnd_lfo.tick_filt_sample_and_hold() * 0.05f * pars.b) + 0.95f;
      wet *= lfov;

      // probably this mixing parameter needs to be adjusted from 0, even though
      // the feedback lines should use the 4 channels...
      auto out
        = mix (pars, wet, pars.b, (float) _n_stages, ins[0][i], ins[1][i]);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  double_x2 mix (
    all_parameters const& p,
    float_x4              wet,
    float                 parallel_mix,
    float                 att,
    T                     dry_l,
    T                     dry_r)
  {
    float dry_gain = 1.f - p.depth;
    float wet_gain = p.depth;
    wet_gain /= ((get_fb_gain (p.feedback, p.feedback_drive)) * att);

    double_x2 wetdbl   = {(double) wet[0], (double) wet[1]};
    double_x2 parallel = {(double) wet[2], (double) wet[3]};
    double_x2 dry      = {(double) dry_l, (double) dry_r};
    wetdbl *= 1.f - abs (parallel_mix);
    wetdbl += parallel_mix * parallel;
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

    uint comb_size
      = std::ceil (max_comb_delay_ms * 0.001f * _srate) + _comb.min_size_spls();
    comb_size = _comb.n_required_elems (comb_size, max_comb_stages);

    schroeder_size *= vec_traits_t<decltype (_scho)::value_type>::size;
    comb_size *= vec_traits_t<decltype (_comb)::value_type>::size;

    _mem.clear();
    _mem.resize (std::max (schroeder_size, comb_size));
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
    float stages;
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
    sigmoid::mystran<4, true>>;
  struct fb_filter {
    enum { locut, hicut, count };
  };

  part_classes<mp_list<zdf_type>, float_x4>                     _feedback;
  part_class_array<allpass_type, float_x4, max_phaser_stages>   _phaser;
  part_classes<filters_list, float_x4>                          _fb_filters;
  part_classes<mp_list<onepole_allpass>, float_x4>              _onepole;
  interpolated_delay_line<float_x4, linear_interp, true, false> _scho;
  interpolated_delay_line<float_x1, sinc_interp<8, 64>, false, false> _comb;
  value_smoother<
    float,
    std::array<float, std::max (max_scho_stages, max_comb_stages)>>
    _del_spls;
  union {
    std::array<float_x4, max_scho_stages> scho;
    std::array<float_x1, max_comb_stages> comb;
  } _g;

  std::array<float_x1, max_comb_stages> _del_g;
  crange<float>                         _delay_mem;

  alignas (sse_bytes) array2d<float, n_ap_channels, max_phaser_stages> _mod;

  std::vector<float_x1, overaligned_allocator<float_x4, 128>> _mem;
  sweep_lfo_type _sweep_lfo; // 0 = L, 1 = R, control rate
  lfo<4>         _rnd_lfo; // 0 = L, 1 = R, audio rate
  float_x4       _1spl_fb;
  uint           _n_processed_samples;
  uint           _control_rate_mask;
  uint           _n_stages;
  float          _lfo_t_spl;
  float          _t_spl;
  float          _srate;
  float          _lp_smooth_coeff;
  float          _beat_hz;
};
//------------------------------------------------------------------------------
} // namespace artv
