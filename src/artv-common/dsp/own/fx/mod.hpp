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
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

#define ARTV_MOD_SCHO_TIRAN     1
#define ARTV_MOD_CHO_FLAN_TIRAN 1

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
    return float_param ("%", 0., 100., 75., 0.001);
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
  struct lfo_warp_tag {};
  void set (lfo_warp_tag, int v)
  {
    _param_smooth.target().lfo_warp = 0.5f - v * 0.01f * 0.5f;
  }

  static constexpr auto get_parameter (lfo_warp_tag)
  {
    return float_param ("%", -100., 100., 0., 0.01);
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
    return float_param ("%", 0., 100., 50., 0.001);
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
    return float_param ("%", -100., 100., 70., 0.1);
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
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    v *= 0.01f;
    v                            = v * v * v * v;
    _param_smooth.target().drive = v * 1400.f;
  }

  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("%", 0., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct drive_curve_tag {};
  void set (drive_curve_tag, float v)
  {
    _param_smooth.target().drive_curve = v *= 0.01f;
  }

  static constexpr auto get_parameter (drive_curve_tag)
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
  struct spread_tag {};
  void set (spread_tag, float v) { _param_smooth.target().spread = v * 0.01f; }

  static constexpr auto get_parameter (spread_tag)
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
    switch (v) {
    case mode::zdf:
      break;
    case mode::schroeder:
      [[fallthrough]];
    case mode::schroeder_nested:
      using TS = decltype (_scho)::value_type;
      _scho.reset (_mem_scho, max_scho_stages);
#if ARTV_MOD_SCHO_TIRAN
      _scho.set_resync_delta (10.0);
#endif
      break;
    case mode::flanger: {
      _flan.reset (_mem_flan, flan_stages);
      _dry.reset (_mem_dry, 1);
#if ARTV_MOD_CHO_FLAN_TIRAN
      _flan.set_resync_delta (10.0);
      _dry.set_resync_delta (10.0);
#else
      // pass the shared sinc interpolator coefficients
      _flan.reset_interpolator (0, false, _sinc_co.to_const());
      _dry.reset_interpolator (0, false, _sinc_co.to_const());
#endif
      _n_stages = flan_stages;
      break;
    }
    case mode::chorus: {
      _chor.reset (_mem_chor, max_chor_stages);
      _dry.reset (_mem_dry, 1);
#if ARTV_MOD_CHO_FLAN_TIRAN
      _chor.set_resync_delta (10.0);
      _dry.set_resync_delta (10.0);
#else
      // pass the shared sinc interpolator coefficients
      _chor.reset_interpolator (0, false, _sinc_co.to_const());
      _dry.reset_interpolator (0, false, _sinc_co.to_const());
#endif
      break;
    }
    default:
      assert (false);
      break;
    }
    if (_param.mode == mode::chorus) {
      // avoid high delay_spls values on change
      for (auto& v : _del_spls.target()) {
        v = 20.f;
      }
      _del_spls.set_all_from_target();
    }
    _1spl_fb             = vec_set<4> (0.f);
    _n_processed_samples = 0; // trigger the control block on first sample
    _param.mode          = v;
  }
  struct mode {
    enum { zdf, schroeder, schroeder_nested, flanger, chorus };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Phaser (ZDF)", "PhaserFlanger", "FlangerPhaser", "Flanger", "Chorus"),
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
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
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
      tick_schroeder<T, false> (outs, ins, samples);
      break;
    case mode::schroeder_nested:
      tick_schroeder<T, true> (outs, ins, samples);
      break;
    case mode::flanger:
      tick_flanger<T> (outs, ins, samples);
      break;
    case mode::chorus:
      tick_chorus<T> (outs, ins, samples);
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
    lfo_warp_tag,
    lfo_wave_tag,
    lfo_stereo_tag,
    center_tag,
    a_tag,
    feedback_tag,
    feedback_locut_tag,
    feedback_hicut_tag,
    drive_tag,
    drive_curve_tag,
    spread_tag,
    b_tag,
    depth_tag,
    mode_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr uint max_phaser_stages = 21;
  static constexpr uint max_scho_stages   = 22;
  static constexpr uint max_scho_delay_ms = 30;
  static constexpr uint max_chor_delay_ms = 48;
  static constexpr uint max_chor_stages   = 8;
  static constexpr uint max_flan_delay_ms = 21;
  static constexpr uint flan_stages       = 2;
  static constexpr uint max_dry_delay_ms  = max_chor_delay_ms / 2;
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
      hz = lfo_rate * 12.f;
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
    float           centerf,
    float           spread,
    float           lfo_depth,
    float           spread_ratio,
    uint            n_stages,
    sweep_lfo_value lfo)
  {
    float modlim = std::min (centerf, 0.5f);
    modlim       = std::min (1.f - centerf, modlim);

    sweep_lfo_value mod    = lfo * lfo_depth * modlim;
    sweep_lfo_value center = centerf + mod;
    // detune by a ramp
    sweep_lfo_value detuned = (1.f - abs (spread * spread_ratio)) * center;
    sweep_lfo_value diff    = (center - detuned) / (float) (n_stages * 2);
    sweep_lfo_value start   = (spread > 0.f) ? detuned : center;
    diff                    = (spread > 0.f) ? diff : -diff;
    f32_x4 curr             = vec_cat (start, start + diff);
    f32_x4 add              = vec_cat (diff, diff) * 2.f;

    for (uint i = 0; i < n_stages; ++i) {
      _mod[i] = vec_to_array (curr);
      curr += add;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_phaser (xspan<T*> outs, xspan<T const*> ins, uint samples)
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
        float stages = 1.f + pars.stages * (max_phaser_stages - 1);
        _n_stages    = (uint) stages;

        run_phaser_mod (
          pars.center,
          pars.spread,
          pars.lfo_depth,
          0.07f,
          _n_stages,
          run_lfo (pars));

        for (uint s = 0; s < _n_stages; ++s) {
          f32_x4 v = vec_from_array (_mod[s]);
          auto   f = lfo_freq_warp (v, 1.f - pars.lfo_warp);
          auto   q = f;
          f        = ((f * 0.9995f) + 0.0005f);
          f *= f32_x4 {21200.f, 21200.f, 20000.f, 20000.f};
          // A lot of freq-dependant weight on the Q when detuning
          q = (1.f + 3.f * q * abs (pars.spread));
          q *= pars.a * pars.a * f32_x4 {1.f, 1.f, 1.6f, 1.6f};
          q += 0.015f;
          q *= 1.f + (_n_stages) * (1.1f / max_phaser_stages);
          _phaser.reset_coeffs_on_idx (s, f, q, _t_spl);
        }
        for (uint s = _n_stages; s < max_phaser_stages; ++s) {
          _phaser.reset_states_on_idx (s);
        }
        f32_x4 f {450.f, 460.f, 670.f, 680.f};
        f *= vec_from_array (_mod[0]);
        _onepole.reset_coeffs<0> (f, _t_spl);
      }

      f32_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      auto   onep_out = _onepole.tick<0> (wet);
      wet             = vec_shuffle (wet, onep_out, 0, 1, 4, 5);

      // obtain response for the phaser
      zdf::response<f32_x4> aps_resp {};
      for (uint s = 0; s < _n_stages; ++s) {
        aps_resp *= _phaser.tick_on_idx (s, zdf::gs_coeffs_tag {});
      };
      // obtain response for the filters
      auto filt_resp
        = _fb_filters.tick<fb_filter::locut> (zdf::gs_coeffs_tag {});
      filt_resp *= _fb_filters.tick<fb_filter::hicut> (zdf::gs_coeffs_tag {});

      // Get ZDF feedback
      wet = _feedback.tick (
        wet,
        aps_resp,
        filt_resp,
        vec_set<4> (pars.feedback),
        vec_set<4> (pars.drive),
        vec_set<4> (pars.drive_curve));

      // Run regular filter processing
      assert (_n_stages > 0);
      for (uint s = 0; s < _n_stages; ++s) {
        wet = _phaser.tick_on_idx (s, wet);
      }
      // just running the shelves to update the states.
      auto fbv = wet * pars.feedback;
      fbv      = zdf_type::nonlin::tick (
        fbv, vec_set<4> (pars.drive), vec_set<4> (pars.drive_curve));
      _fb_filters.tick_cascade (fbv);

      auto out = mix (
        wet,
        f64_x2 {ins[0][i], ins[1][i]},
        pars.depth / get_fb_gain (pars.feedback, pars.drive),
        1.f - pars.depth,
        pars.b);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T, bool nested>
  void tick_schroeder (xspan<T*> outs, xspan<T const*> ins, uint samples)
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
        float stages = 2.f + pars.stages * (max_scho_stages - 2);
        _n_stages    = ((uint) stages) / 2;
        run_phaser_mod (
          (1.f - pars.center) / _n_stages,
          pars.spread,
          pars.lfo_depth,
          0.55f,
          _n_stages,
          run_lfo (pars));

        for (uint s = 0; s < _n_stages; ++s) {
          for (uint j = 0; j < 2; ++j) {
            uint            idx        = s * 2 + j;
            constexpr float sec_factor = (max_scho_delay_ms * 0.001f);
            float           t          = _mod[s][j];
            t                          = lfo_time_warp (t, pars.lfo_warp);
            t                          = _srate * t * sec_factor;
            t                          = std::clamp (
              t,
              (float) _scho.min_delay_spls(),
              (float) _scho.max_delay_spls());
            _del_spls.target()[idx] = t;
            f32_x2          v       = (j == 0) ? f32_x2 {_mod[s][1], _mod[s][3]}
                                               : f32_x2 {_mod[s][0], _mod[s][2]};
            constexpr float min_g   = 0.0001f;
            auto            g       = pars.a * 0.88f + min_g;
            // gain mod
            v            = (1.f - g) * (v - 0.5f) * 1.8f * pars.b * pars.b;
            _g.scho[idx] = g + v;
          }
        }
        for (uint s = (_n_stages * 2); s < _del_spls.target().size(); ++s) {
          _del_spls.target()[s] = _del_spls.target()[s - 1];
        }
        f32_x4 fv {1460.f, 1460.f, 1670.f, 1780.f};
        auto   f = vec_from_array (_mod[_n_stages / 2]);
        f *= f;
        f *= (1.f + pars.a) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      _del_spls.tick();

      f32_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      std::array<f32_x2, 2> sig {{{wet[0], wet[2]}, {wet[1], wet[3]}}};
      std::array<f32_x2, max_scho_stages> to_push {};
      auto&                               del_spls = _del_spls.get();
      if constexpr (!nested) {
        for (uint s = 0; s < _n_stages; ++s) {
          for (uint c = 0; c < n_channels; ++c) {
            uint idx     = s * 2 + c;
            auto yn      = _scho.get (del_spls[idx], idx);
            auto r       = allpass_fn::tick<f32_x2> (sig[c], yn, _g.scho[idx]);
            sig[c]       = r.out;
            to_push[idx] = r.to_push;
          }
        }
      }
      else {
        // Nested, as e.g. a lattice
        std::array<f32_x2, 2> fwd {};
        for (uint c = 0; c < n_channels; ++c) {
          uint idx = c;
          auto yn  = _scho.get (del_spls[idx], idx);

          fwd[c] = sig[c] + yn * _g.scho[idx];
          sig[c] = yn + fwd[c] * -_g.scho[idx];
        }
        for (uint s = 1; s < _n_stages; ++s) {
          for (uint c = 0; c < n_channels; ++c) {
            uint idx = s * 2 + c;
            auto yn  = _scho.get (del_spls[idx], idx);

            fwd[c] += yn * _g.scho[idx];
            to_push[idx - 2] = yn + fwd[c] * -_g.scho[idx];
          }
        }
        to_push[(_n_stages * 2) - 2] = fwd[0];
        to_push[(_n_stages * 2) - 1] = fwd[1];
        sig[0]                       = -sig[0];
        sig[1]                       = -sig[1];
      }
      _scho.push (to_push);
      wet = vec_shuffle (sig[0], sig[1], 0, 2, 1, 3);

      auto fb_val = zdf_type::nonlin::tick (
        wet, vec_set<4> (pars.drive), vec_set<4> (pars.drive_curve));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);

      auto out = mix (
        wet,
        f64_x2 {ins[0][i], ins[1][i]},
        pars.depth / get_fb_gain (pars.feedback, pars.drive),
        1.f - pars.depth,
        1.f);
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_chorus (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {
      _param_smooth.tick();
      // access parameters without caring if they are smoother or not.
      all_parameters pars;
      *((smoothed_parameters*) &pars)   = _param_smooth.get();
      *((unsmoothed_parameters*) &pars) = _param;

      auto rndlfo = _rnd_lfo.tick_filt_sample_and_hold();

      if ((_n_processed_samples & _control_rate_mask) == 0) {
        // stages are in pairs (TODO check)
        float stages
          = 1.f + pars.stages * ((float) (max_chor_stages / 2) - 1.0000000001f);
        _n_stages      = (uint) stages;
        _n_stages_frac = stages - (float) _n_stages;
        _n_stages += 1; // fractional stage

        run_phaser_mod (
          0.05f * rndlfo[0] + 0.85f - pars.center, // reverse range lf to hf
          pars.spread,
          pars.lfo_depth * (0.5f + 0.5f * pars.center * pars.center),
          0.8f + 0.1f * rndlfo[0],
          _n_stages,
          run_lfo (pars));

        for (uint s = 0; s < (_n_stages * 2); ++s) {
          constexpr float msec_offset = 9.f;
          constexpr float sec_factor
            = ((max_chor_delay_ms - msec_offset) * 0.001f);

          bool  neg = !!((s / 2) % 2);
          float t   = _mod[s][0];
          t         = lfo_time_warp (t, pars.lfo_warp);
          t         = neg ? 1.f - t : t;
          t         = _srate * (msec_offset * 0.001f + (t * sec_factor));
          _del_spls.target()[s]    = t;
          constexpr float g_factor = 0.1f;
          auto            m        = _mod[s][1];
          auto            gv       = pars.feedback;
          auto            negf     = neg ? 1.f : -1.f;
          _g.chor[s][0]            = gv * gv * 0.87f + m * g_factor * gv;
          _g.chor[s][0] *= neg;
          _rnd_lfo.set_freq (vec_set<4> (0.8f + abs (pars.b)), _t_spl);

          // set the phasers (no feeedback)
          auto v = vec_from_array (_mod[s]);
          auto f = v * v;
          auto q = v;
          f      = ((f * (0.45f + 0.3f * rndlfo)) + 0.15f);
          f *= f32_x4 {21200.f, 21200.f, 17000.f, 17000.f};
          q = (1.f + 2.f * q * abs (pars.spread));
          q += 0.09f;
          _phaser.reset_coeffs_on_idx (s, f, q, _t_spl);
        }
        for (uint s = (_n_stages * 2); s < _del_spls.target().size(); ++s) {
          _del_spls.target()[s] = _del_spls.target()[s - 1];
        }
        f32_x4 fv {460.f, 460.f, 970.f, 980.f};
        auto   f = vec_from_array (_mod[_n_stages / 2]);
        f *= f;
        f *= (1.f + pars.feedback) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      // create some difference on the inputs
      f32_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      _del_spls.tick();
      auto& del_spls = _del_spls.get();
      auto  chor_in  = vec1_array_wrap (vec_to_array (wet));
      wet            = vec_set<4> (0.f);

      // To allow stage crossfading constant panning positions are kept.
      // Processing is done in pairs.
      //
      // pan position   | L             R
      // stage ordering | 1 3 2 4 4 2 3 1

      constexpr auto kpanl = make_array (0.f, 2.f / 7.f, 1.f / 7.f, 3.f / 7.f);
      static_assert (kpanl.size() == (max_chor_stages / 2));

      auto rndmod = (rndlfo * 0.05f * abs (pars.b)) + 0.95f;

      auto panlv = vec_from_array (kpanl);
      // skew towards center based on pars.a
      auto st = pars.lfo_stereo * (1.f / get_parameter (lfo_stereo_tag {}).max);
      panlv *= st;
      panlv += (1.f - st) * 0.5f;
      // using -x^2+2x as a cheap approximation of the sin(x*pi/2) pan law.
      panlv = -panlv * panlv + 2.f * panlv;
      // some randomization
      panlv *= rndmod;
      auto panrv = 1.f - panlv;

      auto panl = vec_to_array (panlv);
      auto panr = vec_to_array (panrv);

      auto trnd = vec_to_array (rndlfo * _srate * 0.0005f);

      auto                                prev = wet;
      std::array<f32_x1, max_chor_stages> to_push {};
      for (uint s = 0; s < _n_stages; ++s) {

        uint idx  = s * 2;
        uint lane = idx % vec_traits_t<f32_x4>::size;

        std::array<float, n_channels> channel;

        for (uint j = 0; j < channel.size(); ++j, ++idx, ++lane) {
          auto n_spls = std::clamp (
            del_spls[idx] + trnd[lane],
            (float) _chor.min_delay_spls(),
            (float) _chor.max_delay_spls());
          auto g       = _g.chor[idx];
          auto yn      = _chor.get (n_spls, idx);
          auto r       = allpass_fn::tick<f32_x1> (chor_in[lane], yn, g);
          to_push[idx] = r.to_push;
          channel[j]   = r.out[0];
        }

        float l = channel[0] * panl[s] + channel[1] * panr[s];
        float r = channel[0] * panr[s] + channel[1] * panl[s];

        // TODO: probably run this on its own loop?
        f32_x4 v {l, r, l, r};
        v    = _phaser.tick_on_idx (s, v);
        prev = wet;
        wet += v;
      }
      //  crossfade last stage
      wet = prev + (wet - prev) * _n_stages_frac;
      _chor.push (to_push);

      // naive saturation
      wet = zdf_type::nonlin::tick (
        wet, vec_set<4> (pars.drive), vec_set<4> (pars.drive_curve) * rndmod);
      wet = _fb_filters.tick_cascade (wet);

      // delay dry signal based
      f32_x2 dry {ins[0][i], ins[1][i]};
      _dry.push (xspan {&dry, 1});
      constexpr auto ratio
        = (float) max_chor_delay_ms / (0.01f + (float) (max_dry_delay_ms));
      float n_spls = del_spls[0] * ratio * pars.a;
      n_spls       = std::clamp (
        n_spls, (float) _dry.min_delay_spls(), (float) _dry.max_delay_spls());
      auto dryz      = _dry.get (n_spls, 0);
      auto dryz_gain = abs (pars.b) * (float) _n_stages;
      wet += f32_x4 {dryz[0], dryz[1], dryz[0], dryz[1]} * dryz_gain;

      auto out = mix (
        wet,
        f64_x2 {dry[0], dry[1]},
        pars.depth / (_n_stages_frac + (float) (_n_stages - 1) + dryz_gain),
        (1.f - pars.depth),
        pars.b);

      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick_flanger (xspan<T*> outs, xspan<T const*> ins, uint samples)
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
        run_phaser_mod (
          1.f - pars.center, // reverse range lf to hf
          pars.spread,
          pars.lfo_depth,
          0.45f,
          1, // flanger is always one stage
          run_lfo (pars));

        for (uint s = 0; s < flan_stages; ++s) {
          static_assert (flan_stages == 2, "this loop assumes s < 2");
          constexpr float sec_factor = (max_flan_delay_ms * 0.001f);

          float t = _mod[0][s];
          t       = lfo_time_warp (t, pars.lfo_warp);
          t       = _srate * t * sec_factor;
          t       = std::clamp (
            t, (float) _flan.min_delay_spls(), (float) _flan.max_delay_spls());
          _del_spls.target()[s]    = t;
          constexpr float g_factor = 0.1f;
          auto            m        = _mod[0][2 + s];
          auto            gv       = pars.stages;
          _g.flan[s][0]            = -gv * gv * 0.75f + m * m * g_factor * gv;
        }
        for (uint s = flan_stages; s < _del_spls.target().size(); ++s) {
          _del_spls.target()[s] = _del_spls.target()[s - 1];
        }
        f32_x4 fv {460.f, 460.f, 970.f, 980.f};
        auto   f = vec_from_array (_mod[0]);
        f *= f;
        f *= (1.f + pars.feedback) * fv;
        _onepole.reset_coeffs<0> (f, _t_spl);
      }
      _del_spls.tick();

      f32_x4 wet {ins[0][i], ins[1][i], ins[0][i], ins[1][i]};
      wet = _onepole.tick<0> (wet);
      wet -= _1spl_fb * pars.feedback;

      auto& del_spls = _del_spls.get();
      auto  flan_in  = vec_split<flan_stages> (wet);

      std::array<f32_x2, flan_stages> to_push {};
      for (uint s = 0; s < flan_stages; ++s) {
        auto yn    = _flan.get (del_spls[s], s);
        auto r     = allpass_fn::tick<f32_x2> (flan_in[s], yn, _g.flan[s]);
        to_push[s] = r.to_push;
        wet[0 + s] = r.out[0];
        wet[2 + s] = r.out[1];
      }
      _flan.push (to_push);

      auto fb_val = zdf_type::nonlin::tick (
        wet, vec_set<4> (pars.drive), vec_set<4> (pars.drive_curve));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);

      f32_x2 dry {ins[0][i], ins[1][i]};
      _dry.push (xspan {&dry, 1});
      constexpr float kdry
        = 0.001f * std::max (max_dry_delay_ms, max_flan_delay_ms / 2);
      float n_spls = _srate * kdry * (1.f - pars.center) * pars.a;
      n_spls       = std::clamp (
        n_spls, (float) _dry.min_delay_spls(), (float) _dry.max_delay_spls());
      dry = _dry.get (n_spls, 0);

      auto out = mix (
        wet,
        f64_x2 {dry[0], dry[1]},
        pars.depth / get_fb_gain (pars.feedback, pars.drive),
        1.f - pars.depth,
        pars.b);

      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  f64_x2 mix (
    f32_x4 wet,
    f64_x2 dry,
    float  wet_gain,
    float  dry_gain,
    float  parallel_mix)
  {
    f64_x2 wetdbl   = {(double) wet[0], (double) wet[1]};
    f64_x2 parallel = {(double) wet[2], (double) wet[3]};
    wetdbl *= 1.f - abs (parallel_mix);
    wetdbl += parallel_mix * parallel;
    return wetdbl * wet_gain + dry * dry_gain;
  }
  //----------------------------------------------------------------------------
  // hardness makes a sqrt sigmmoid
  float get_fb_gain (float fb, float hardness)
  {
    float lim = zdf_type::nonlin::limit_inf<f32_x1> (
      make_vec (hardness), make_vec (0.f))[0];
    return 1.f + std::min (lim, abs (fb));
  }
  //----------------------------------------------------------------------------
  void reset_mem()
  {
    uint schroeder_size
      = std::ceil (max_scho_delay_ms * 0.001f * _srate) + _scho.min_size_spls();
    schroeder_size
      = _scho.n_required_elems (schroeder_size + 1, max_scho_stages);

    uint chor_size
      = std::ceil (max_chor_delay_ms * 0.001f * _srate) + _chor.min_size_spls();
    chor_size = _chor.n_required_elems (chor_size + 1, max_chor_stages);

    uint flan_size
      = std::ceil (max_flan_delay_ms * 0.001f * _srate) + _flan.min_size_spls();
    flan_size = _flan.n_required_elems (flan_size + 1, flan_stages);

    uint dry_size = std::ceil (max_dry_delay_ms * 0.001f * _srate);
    dry_size      = _dry.n_required_elems (dry_size + 1, n_channels);

    using T_mem = decltype (_mem)::value_type;
#if !ARTV_MOD_CHO_FLAN_TIRAN
    using T_sinc = decltype (_sinc_co)::value_type;
#endif
    using T_scho = decltype (_scho)::value_type;
    using T_chor = decltype (_chor)::value_type;
    using T_flan = decltype (_flan)::value_type;
    using T_dry  = decltype (_dry)::value_type;

    schroeder_size *= vec_traits_t<T_scho>::size;
    chor_size *= vec_traits_t<T_chor>::size;
    dry_size *= vec_traits_t<T_dry>::size;
    flan_size *= vec_traits_t<T_flan>::size;

    auto max_size = std::max (schroeder_size, chor_size + dry_size);
    max_size      = std::max (max_size, flan_size + dry_size);

    _mem.clear();
#if !ARTV_MOD_CHO_FLAN_TIRAN
    _mem.resize (sinc_t::n_coeffs + max_size);
#else
    _mem.resize (max_size);
#endif

    auto mem = xspan {_mem};
#if !ARTV_MOD_CHO_FLAN_TIRAN
    static_assert (std::is_same_v<T_sinc, T_mem>);
    // overaligned to 128, so the tables aren't on cache line boundaries
    _sinc_co = mem.cut_head (sinc_t::n_coeffs);
#endif
    _mem_scho = mem.get_head (schroeder_size).cast<T_scho>();
    _mem_dry  = mem.cut_head (dry_size).cast<T_dry>();
    _mem_chor = mem.get_head (chor_size).cast<T_chor>();
    _mem_flan = mem.get_head (flan_size).cast<T_flan>();

#if !ARTV_MOD_CHO_FLAN_TIRAN
    // initialize windowed sinc table
    sinc_t::reset_coeffs (_sinc_co, 0.45f, 140.f);
#endif
  }
//----------------------------------------------------------------------------
#if 0
// from linear to  cubic.
  template <class T>
  static T lfo_freq_warp (T v, float warp)
  {
    return (1.f - warp) * v + warp * v * v * v;
  }

  //----------------------------------------------------------------------------
  // from cubic to inverse cubic.  0.5 warp = crosses 0.5,0.5 but snakey.
  template <class T>
  static T lfo_freq_warp (T v, float warp)
  {
    auto v1   = v * v * v;
    auto vm12 = v - 1.f;
    auto v2   = vm12 * vm12 * vm12 + 1.f;

    return (1.f - warp) * v2 + warp * v1;
  }
    //----------------------------------------------------------------------------
  // from squared to inverse squared. 0.5 warp = linear
  template <class T>
  static T lfo_freq_warp (T v, float warp)
  {
    auto v1   = v * v;
    auto vm12 = v - 1.f;
    auto v2   = 1.f - vm12 * vm12;

    return (1.f - warp) * v2 + warp * v1;
  }
    //----------------------------------------------------------------------------
  // from quadractic to inverse squared. 0.5 warp = crosses 0.5,0.5 but snakey.
  template <class T>
  static T lfo_freq_warp (T v, float warp)
  {
    auto v1   = v * v * v * v;
    auto vm12 = v - 1.f;
    auto v2   = 1.f - vm12 * vm12;

    return (1.f - warp) * v2 + warp * v1;
  }
#endif
  //----------------------------------------------------------------------------
  // from linear to cubic.
  template <class T>
  static T lfo_freq_warp (T v, float warp)
  {
    return (1.f - warp) * v * v * v + warp * v;
  }
  //----------------------------------------------------------------------------
  // from linear to cubic. but weighted towards the cubic side
  template <class T>
  static T lfo_time_warp (T v, float warp)
  {
    warp *= warp;
    return (1.f - warp) * v * v * v + warp * v;
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
    float lfo_warp;
    float lfo_stereo;
    float center;
    float a;
    float b;
    float depth; // aka mix (make explicit?)
    float feedback;
    float drive;
    float drive_curve;
    float spread;
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

  using sinc_t = sinc_interp<8, 96>;

  part_classes<mp_list<zdf_type>, f32_x4>                   _feedback;
  part_class_array<allpass_type, f32_x4, max_phaser_stages> _phaser;
  part_classes<filters_list, f32_x4>                        _fb_filters;
  part_classes<mp_list<onepole_allpass>, f32_x4>            _onepole;
#if ARTV_MOD_SCHO_TIRAN
  modulable_thiran1_delay_line<f32_x2, 4, true, false> _scho {};
#else
  interpolated_delay_line<f32_x2, linear_interp, true, false> _scho;
#endif
#if ARTV_MOD_CHO_FLAN_TIRAN
  modulable_thiran1_delay_line<f32_x1, 4, true, false> _chor;
  modulable_thiran1_delay_line<f32_x2, 4, true, false> _dry;
  modulable_thiran1_delay_line<f32_x2, 4, true, false> _flan;
#else
  interpolated_delay_line<f32_x1, sinc_t, true, false, true> _chor;
  interpolated_delay_line<f32_x2, sinc_t, true, false, true> _dry;
  interpolated_delay_line<f32_x2, sinc_t, true, false, true> _flan;
#endif
  value_smoother<
    float,
    std::array<float, std::max (max_scho_stages, max_chor_stages)>>
    _del_spls;
  union {
    std::array<f32_x2, max_scho_stages> scho;
    std::array<f32_x1, max_chor_stages> chor;
    std::array<f32_x2, flan_stages>     flan;
  } _g;

  std::array<f32_x1, max_chor_stages> _del_g;
  xspan<float>                        _delay_mem;

  alignas (sse_bytes) array2d<float, n_ap_channels, max_phaser_stages> _mod;

  std::vector<float, overaligned_allocator<float, 128>> _mem;
  sweep_lfo_type _sweep_lfo; // 0 = L, 1 = R, control rate
  lfo<4>         _rnd_lfo; // 0 = L, 1 = R, audio rate
  f32_x4         _1spl_fb;
  uint           _n_processed_samples;
  uint           _control_rate_mask;
  uint           _n_stages;
  float          _n_stages_frac;
  float          _lfo_t_spl;
  float          _t_spl;
  float          _srate;
  float          _lp_smooth_coeff;
  float          _beat_hz;
#if !ARTV_MOD_CHO_FLAN_TIRAN
  xspan<float> _sinc_co;
#endif
  xspan<f32_x2> _mem_dry;
  xspan<f32_x1> _mem_chor;
  xspan<f32_x2> _mem_flan;
  xspan<f32_x2> _mem_scho;
};
//------------------------------------------------------------------------------
} // namespace artv
