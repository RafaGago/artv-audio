#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/filters/andy_smooth.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/biquad.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/zdf.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/oscillators/phasor.hpp"
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
#include "artv-common/misc/vec_math.hpp"
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
  struct mod_warp_tag {};
  void set (mod_warp_tag, int v)
  {
    _param_smooth.target().mod_warp = 0.5f - v * 0.01f * 0.5f;
  }

  static constexpr auto get_parameter (mod_warp_tag)
  {
    return float_param ("%", -100., 100., 50., 0.01);
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
  struct env_depth_tag {};
  void set (env_depth_tag, float v)
  {
    v *= 0.01;
    _param.env_depth = v;
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (env_depth_tag)
  {
    return float_param ("%", -100., 100., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct env_speed_tag {};
  void set (env_speed_tag, float v)
  {
    v *= 0.01;
    if (v == _param.env_speed) {
      return;
    }
    _param.env_speed = v;
    v += 0.005;
    _env.reset_coeffs<env::fast> (
      vec_set<2> (0.0003f + v * 0.05f), vec_set<2> (0.0009f + v), _ctrl_t_spl);
    _env.reset_coeffs<env::slow> (
      vec_set<2> (0.0006f + v * 0.4f),
      vec_set<2> (0.0018f + v * 8.f),
      _ctrl_t_spl);
    xspan_memcpy (
      _env.get_coeffs<env::smooth1>(), _env.get_coeffs<env::fast>());
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (env_speed_tag)
  {
    return float_param ("%", 0., 100., 20., 0.001);
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
  struct stereo_tag {};
  void set (stereo_tag, float v) { _param_smooth.target().stereo = v * 0.01f; }

  static constexpr auto get_parameter (stereo_tag)
  {
    return float_param ("%", -100., 100., 0., 0.1);
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
  void set (feedback_tag, float in)
  {
    auto v = abs (in);
    v *= 0.01;
    // parabola
    v -= 1.f;
    v = 1.f - v * v;
    v = std::copysign (v, in);

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
      break;
    case mode::flanger: {
      _flan.reset (_mem_flan, flan_stages);
      _dry.reset (_mem_dry, 1);
#if !ARTV_MOD_CHO_FLAN_TIRAN
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
#if !ARTV_MOD_CHO_FLAN_TIRAN
      // pass the shared sinc interpolator coefficients
      _chor.reset_interpolator (0, false, _sinc_co.to_const());
      _dry.reset_interpolator (0, false, _sinc_co.to_const());
#endif
      break;
    }
    default:
      assert (false);
      return;
    }
    if (_param.mode == mode::chorus) {
      // avoid high delay_spls values on change
      for (auto& v : _del_spls.target()) {
        v = 20.f;
      }
      _del_spls.set_all_from_target();
    }
    _delay.reset (_mem_delay, n_channels);
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
  static constexpr uint max_delay_4ths = 8;
  struct delay_time_tag {};
  void set (delay_time_tag, float v)
  {
    _param_smooth.target().del_4beats = std::max (v, 0.02f);
  }

  static constexpr auto get_parameter (delay_time_tag)
  {
    return float_param ("quarters", 0., max_delay_4ths, 3., 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_gain_tag {};
  void set (delay_gain_tag, float v)
  {
    _param_smooth.target().del_gain = v * 0.01f;
  }

  static constexpr auto get_parameter (delay_gain_tag)
  {
    return float_param ("%", -100., 100, 50., 0.1);
  }
  //----------------------------------------------------------------------------
  struct delay_mode_tag {};
  void set (delay_mode_tag, int v) { _param.delay_mode = v; }

  struct delay_mode {
    enum {
      off,
      stereo,
      stereo_inverted,
      pingpong_l,
      pingpong_r,
      pingpong,
      pongping
    };
  };

  static constexpr auto get_parameter (delay_mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Off",
        "Stereo",
        "Stereo-inv",
        "Ping-Pong L",
        "Ping-Pong R",
        "Ping-Pong",
        "Pong-Ping"),
      20);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _srate               = pc.get_sample_rate();
    _t_spl               = 1.f / _srate;
    _beat_hz             = pc.get_play_state().bpm * (1.f / 60.f);
    _1_4beat_spls        = (0.25f / _beat_hz) * _srate;
    _n_processed_samples = 0;
    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 5;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _ctrl_t_spl          = _t_spl * (_control_rate_mask + 1);
    _1spl_fb             = decltype (_1spl_fb) {};

    _phaser.reset_states_cascade();
    _fb_filters.reset_states_cascade();
    _feedback.reset_states_cascade();
    _onepole.reset_states_cascade();
    _env.reset_states_cascade();
    _param_smooth.reset (_t_spl, 10.f);
    _del_spls.reset (_t_spl, 10.f);

    reset_mem();
    using phase = decltype (_rnd_lfo)::phase_type;
    _rnd_lfo.set_phase (phase {phase::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    _rnd_lfo.set_freq (vec_set<4> (0.3f), _t_spl);
    _tremolo_lfo.reset();
    _env.reset_coeffs<env::smooth2> (
      vec_set<2> (1.f), vec_set<2> (0.05f), _ctrl_t_spl);

    mp11::mp_for_each<parameters> ([&] (auto param) {
      set (param, get_parameter (param).min);
      if constexpr (!is_choice<decltype (get_parameter (param))>) {
        set (param, get_parameter (param).max); // max might be not yet impl.
      }
      else {
        set (param, get_parameter (param).min + 1);
      }
      set (param, get_parameter (param).defaultv);
    });
    _param_smooth.set_all_from_target();
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
      tick<T> (
        outs,
        ins,
        samples,
        [this] (auto&&... args) { phaser_control_block<T> (args...); },
        [this] (auto&&... args) { phaser_tick<T> (args...); });
      break;
    case mode::schroeder:
      tick<T> (
        outs,
        ins,
        samples,
        [this] (auto&&... args) { schroeder_control_block<T> (args...); },
        [this] (auto&&... args) { schroeder_tick<T, false> (args...); });
      break;
    case mode::schroeder_nested:
      tick<T> (
        outs,
        ins,
        samples,
        [this] (auto&&... args) { schroeder_control_block<T> (args...); },
        [this] (auto&&... args) { schroeder_tick<T, true> (args...); });
      break;
    case mode::flanger:
      tick<T> (
        outs,
        ins,
        samples,
        [this] (auto&&... args) { flanger_control_block<T> (args...); },
        [this] (auto&&... args) { flanger_tick<T> (args...); });
      break;
    case mode::chorus:
      tick<T> (
        outs,
        ins,
        samples,
        [this] (auto&&... args) { chorus_control_block<T> (args...); },
        [this] (auto&&... args) { chorus_tick<T> (args...); });
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
    mod_warp_tag,
    lfo_wave_tag,
    stereo_tag,
    center_tag,
    a_tag,
    feedback_tag,
    feedback_locut_tag,
    feedback_hicut_tag,
    drive_tag,
    spread_tag,
    b_tag,
    depth_tag,
    mode_tag,
    env_depth_tag,
    env_speed_tag,
    delay_time_tag,
    delay_gain_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size    = 16;
  static constexpr uint max_phaser_stages = 21;
  static constexpr uint max_scho_stages   = 22;
  static constexpr uint max_scho_delay_ms = 30;
  static constexpr uint max_chor_delay_ms = 45;
  static constexpr uint max_chor_stages   = 8;
  static constexpr uint max_flan_delay_ms = 21;
  static constexpr uint flan_stages       = 2;
  static constexpr uint max_dry_delay_ms  = 30;
  static constexpr uint n_channels        = 2;
  static constexpr uint n_ap_channels     = 4;
  //----------------------------------------------------------------------------
  struct all_parameters;
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  using sweep_lfo_type  = lfo<n_channels>;
  using sweep_lfo_value = sweep_lfo_type::value_type;
  //----------------------------------------------------------------------------
  sweep_lfo_value run_lfo (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par)
  {
    float hz;
    switch (par.lfo_time_base) {
    case lfo_time_base::free: {
      float lfo_rate = s_par.lfo_rate * 0.01f;
      lfo_rate       = 1.f - lfo_rate;
      lfo_rate *= lfo_rate * lfo_rate;
      hz = lfo_rate * 12.f;
    } break;
    case lfo_time_base::quarter_beat:
      // a quarter beat for each 10%
      hz = (4.f * _beat_hz) / (1.f + s_par.lfo_rate * 0.1f);
      break;
    case lfo_time_base::two_beats:
      // two beats for each 10%
      hz = (1 / 2.f * _beat_hz) / (1.f + s_par.lfo_rate * 0.1f);
      break;
    default:
      assert (false);
      break;
    }

    _sweep_lfo.set_freq (vec_set<2> (hz), _ctrl_t_spl);
    _tremolo_lfo.set_freq (vec_set<4> (hz), _t_spl);
    auto stereo_ph
      = phase<1> {phase_tag::shifted_bipolar {}, s_par.stereo}.get_raw (0);
    if (!par.lfo_off) {
      // updating stereo phase diff.
      auto ph = _sweep_lfo.get_phase();
      ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
      _sweep_lfo.set_phase (ph);
      auto tph = _tremolo_lfo.get_phase();
      tph.set_raw (tph.get_raw (0) + stereo_ph, 1);
      tph.set_raw (tph.get_raw (0), 2);
      tph.set_raw (tph.get_raw (1), 3);
      _tremolo_lfo.set_phase (tph);
    }
    else {
      auto start_ph = phase<n_channels> {phase_tag::degrees {}, 0.f, 0.f};
      start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
      _sweep_lfo.set_phase (start_ph);
    }

    sweep_lfo_value lfov;
    switch (par.lfo_wave) {
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
    return lfov * s_par.lfo_depth;
  }
  //----------------------------------------------------------------------------
  sweep_lfo_value run_env (unsmoothed_parameters const& par, sweep_lfo_value in)
  {
    if (par.env_depth == 0.f) {
      return sweep_lfo_value {};
    }
    auto fast = _env.tick<env::fast> (in * in);
    auto slow = _env.tick<env::slow> (in * in);
    fast -= slow; // bipolar
    slow = zero_to_lowest (slow); // NaN "sweep_lfo_value" is based on builtins
    auto ret   = fast / slow;
    ret        = _env.tick<env::smooth1> (ret);
    ret        = _env.tick<env::smooth2> (ret);
    float fact = (par.env_depth < 0) ? -4.f : 4.f;
    ret *= par.env_depth * par.env_depth * fact;
    ret = vec_clamp (ret, -1.f, 1.f);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class T>
  sweep_lfo_value run_mod_srcs (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par,
    T                            l,
    T                            r)
  {
    auto lfov = run_lfo (par, s_par);
    auto envv = run_env (par, f32_x2 {(float) l, (float) r});
    return vec_clamp (lfov + envv, -1.f, 1.f);
  }
  //----------------------------------------------------------------------------
  void compute_mod (
    float           centerf,
    float           spread,
    float           lfo_depth,
    float           spread_ratio,
    float           spread_mod,
    uint            n_stages,
    sweep_lfo_value mod)
  {
    float modlim = std::min (centerf, 0.5f);
    modlim       = std::min (1.f - centerf, modlim);

    auto mod_unipolar = 0.5f + mod * 0.5f;
    mod_unipolar      = (spread_mod < 0) ? 1.f - mod_unipolar : mod_unipolar;
    spread_mod        = abs (spread_mod);
    auto spread_r     = (1.f - spread_mod) + (spread_mod * mod_unipolar);
    spread_r *= spread_ratio;

    mod *= modlim;
    sweep_lfo_value center = centerf + mod;
    // detune by a ramp
    sweep_lfo_value detuned = (1.f - vec_abs (spread * spread_r)) * center;
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
  void delay_push (
    xspan<smoothed_parameters const>            s_par,
    uint                                        mode,
    xspan<std::array<f32_x1, n_channels> const> taps_tail,
    xspan<f32_x4 const>                         wet,
    xspan<T const>                              inl,
    xspan<T const>                              inr)
  {
    auto const block_size = wet.size();
    auto const push
      = [this, s_par] (std::array<f32_x1, n_channels> spls, uint i) {
          float max_g    = 0.999999f - abs (s_par[i].feedback);
          float g        = s_par[i].del_gain;
          auto  del_gain = max_g * std::copysign (g * g, g);
          spls[0] *= del_gain;
          spls[1] *= del_gain;
          _delay.push (spls);
        };
    std::array<f32_x1, n_channels> spls;
    switch (mode) {
    case delay_mode::stereo:
      // tail already included in the signal
      for (uint i = 0; i < block_size; ++i) {
        spls[0][0] = wet[i][0];
        spls[1][0] = wet[i][1];
        push (spls, i);
      }
      break;
    case delay_mode::stereo_inverted:
      // tail already included in the signal
      for (uint i = 0; i < block_size; ++i) {
        spls[0][0] = wet[i][1];
        spls[1][0] = wet[i][0];
        push (spls, i);
      }
      break;
      // self note, this is not a regular delay but a modulation effect. Both
      // channels (L-R) are unconditionally pased through the fx chain, while
      // a real ping pong delay would only get (mid) signal on one of the
      // lanes.
      //
      // The requirements are:
      // -The channels at the feedback point (this call) have to be a stereo
      //  phaser/flanger/chorus
      // -The effects have to be on the feedback loop of the delay.
      //
      // As the whole FX chain has to be run as stereo, the possibility to do
      // feed the mid signal is lost, as the other channel is also mixed with
      // the input, so the only way to get ping pong with the full wet signal
      // is to take one of the sides. To get the mid channel a compromise has
      // to be made and to take opposite side dry.
    case delay_mode::pingpong_l:
      for (uint i = 0; i < block_size; ++i) {
        spls[0][0] = wet[i][1];
        spls[1]    = taps_tail[i][0];
        push (spls, i);
      }
      break;
    case delay_mode::pingpong_r:
      for (uint i = 0; i < block_size; ++i) {
        spls[0]    = taps_tail[i][1];
        spls[1][0] = wet[i][0];
        push (spls, i);
      }
      break;
    case delay_mode::pingpong:
      for (uint i = 0; i < block_size; ++i) {
        spls[0][0] = wet[i][1] + inl[i];
        spls[1]    = taps_tail[i][0];
        push (spls, i);
      }
      break;
    case delay_mode::pongping:
      for (uint i = 0; i < block_size; ++i) {
        spls[0]    = taps_tail[i][1];
        spls[1][0] = wet[i][0] + inr[i];
        push (spls, i);
      }
      break;
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class FC, class FP>
  void tick (
    xspan<T*>       outs,
    xspan<T const*> ins,
    uint            samples,
    FC              control_func,
    FP              process_func)
  {
    std::array<f32_x4, max_block_size>              wet;
    array2d<f32_x1, n_channels, max_block_size>     delay_taps_tail;
    std::array<smoothed_parameters, max_block_size> s_par;
    // signal no changes/reloads on the whole block
    unsmoothed_parameters par = _param;

    auto inl  = xspan {ins[0], samples};
    auto inr  = xspan {ins[1], samples};
    auto outl = xspan {outs[0], samples};
    auto outr = xspan {outs[1], samples};

    while (inl.size()) {
      uint n_next_ctrl = ~_n_processed_samples & _control_rate_mask;
      bool run_ctrl    = (n_next_ctrl == _control_rate_mask); // was zero
      n_next_ctrl += 1;
      uint block_size = std::min<uint> (max_block_size, inl.size());
      block_size      = std::min<uint> (block_size, n_next_ctrl);

      if (run_ctrl) {
        control_func (par, _param_smooth.get(), inl[0], inr[0]);
      }
      if (par.delay_mode != delay_mode::off) {
        for (uint i = 0; i < block_size; ++i) {
          _param_smooth.tick();
          s_par[i]   = _param_smooth.get();
          auto& tail = delay_taps_tail[i];
          float spls = _1_4beat_spls * s_par[i].del_4beats;
          tail[0]    = _delay.get (spls, 0, i);
          tail[1]    = _delay.get (spls, 1, i);
          wet[i]     = vec_cat (tail[0], tail[1], tail[0], tail[1]);
        }
      }
      else {
        memset (wet.data(), 0, sizeof wet);
        for (uint i = 0; i < block_size; ++i) {
          _param_smooth.tick();
          s_par[i] = _param_smooth.get();
        }
      }
      process_func (
        outl.cut_head (block_size),
        outr.cut_head (block_size),
        xspan {wet.data(), block_size},
        inl.cut_head (block_size),
        inr.cut_head (block_size),
        xspan {delay_taps_tail.data(), block_size},
        par,
        xspan {s_par.data(), block_size});

      _n_processed_samples += block_size;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void phaser_control_block (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par,
    T                            inl,
    T                            inr)
  {
    float stages = 1.f + s_par.stages * (max_phaser_stages - 1);
    _n_stages    = (uint) stages;

    compute_mod (
      s_par.center,
      s_par.spread,
      s_par.lfo_depth,
      0.07f,
      0.f,
      _n_stages,
      run_mod_srcs (par, s_par, inl, inr));

    for (uint s = 0; s < _n_stages; ++s) {
      f32_x4 v = vec_from_array (_mod[s]);
      auto   f = mod_freq_warp (v, 1.f - s_par.mod_warp);
      auto   q = f;
      f        = ((f * 0.9995f) + 0.0005f);
      f *= f32_x4 {21200.f, 21200.f, 20000.f, 20000.f};
      // A lot of freq-dependant weight on the Q when detuning
      q = (1.f + 3.f * q * abs (s_par.spread));
      q *= s_par.a * s_par.a * f32_x4 {1.f, 1.f, 1.6f, 1.6f};
      q += 0.015f;
      q *= 1.f + (_n_stages) * (1.1f / max_phaser_stages);
      _phaser.reset_coeffs_on_idx (s, f, q, _t_spl);
    }
    for (uint s = _n_stages; s < max_phaser_stages; ++s) {
      _phaser.reset_states_on_idx (s);
    }
    f32_x4 f {250.f, 250.f, 340.f, 341.f};
    f *= 1.5f + vec_from_array (_mod[0]);
    _onepole.reset_coeffs<0> (f, _t_spl);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void phaser_tick (
    xspan<T>       outl,
    xspan<T>       outr,
    xspan<f32_x4>  wet, // might contain the delay feedback
    xspan<T const> inl,
    xspan<T const> inr,
    xspan<std::array<f32_x1, n_channels> const> taps_tail,
    unsmoothed_parameters const&                par,
    xspan<smoothed_parameters const>            s_par)
  {
    uint const block_size = wet.size();

    // gathering into vectors, summing with delay + hp block
    for (uint i = 0; i < block_size; ++i) {
      auto mid = (inl[i] + inr[i]) * 0.5f;
      wet[i] += f32_x4 {inl[i], inr[i], mid, mid};
      auto hp_out = _onepole.tick<0> (wet[i])[1];
      wet[i]      = vec_shuffle (wet[i], hp_out, 0, 1, 6, 7);
    }
    // main filtering block
    for (uint i = 0; i < block_size; ++i) {
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
      wet[i] = _feedback.tick (
        wet[i],
        aps_resp,
        filt_resp,
        vec_set<4> (s_par[i].feedback),
        vec_set<4> (s_par[i].drive),
        vec_set<4> (tanh_like_hardness));
      // Run regular filter processing
      assert (_n_stages > 0);
      for (uint s = 0; s < _n_stages; ++s) {
        wet[i] = _phaser.tick_on_idx (s, wet[i]);
      }
      // just running the shelves to update the states.
      auto fbv = wet[i] * s_par[i].feedback;
      fbv      = zdf_type::nonlin::tick (
        fbv, vec_set<4> (s_par[i].drive), vec_set<4> (tanh_like_hardness));
      _fb_filters.tick_cascade (fbv);
    }
    // delay_block
    if (par.delay_mode != delay_mode::off) {
      delay_push (s_par, par.delay_mode, taps_tail, wet, inl, inr);
    }
    // final mixing
    for (uint i = 0; i < block_size; ++i) {
      auto out = mix (
        wet[i],
        f64_x2 {inl[i], inr[i]},
        s_par[i].depth / get_fb_gain (s_par[i].feedback, s_par[i].drive),
        1.f - s_par[i].depth,
        s_par[i].b);
      outl[i] = out[0];
      outr[i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void schroeder_control_block (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par,
    T                            inl,
    T                            inr)
  {
    float stages = 2.f + s_par.stages * (max_scho_stages - 2);
    _n_stages    = ((uint) stages) / 2;

    compute_mod (
      (1.f - s_par.center) / _n_stages,
      s_par.spread,
      s_par.lfo_depth,
      0.5f,
      s_par.b,
      _n_stages,
      run_mod_srcs (par, s_par, inl, inr));

    constexpr float min_g         = 0.0001f;
    constexpr float max_g         = 0.88f;
    auto const      g             = s_par.a * max_g + min_g;
    auto const      mod_right_max = max_g - g;
    auto const      mod_left_max  = -g;
    auto            modmax = (s_par.b < 0.f) ? mod_left_max : mod_right_max;
    modmax *= s_par.b * s_par.b;

    for (uint s = 0; s < _n_stages; ++s) {
      for (uint j = 0; j < 2; ++j) {
        uint            idx        = s * 2 + j;
        constexpr float sec_factor = (max_scho_delay_ms * 0.001f);
        float           t          = _mod[s][j];
        t                          = mod_time_warp (t, s_par.mod_warp);
        t                          = _srate * t * sec_factor;
        t                          = std::clamp (
          t, (float) _scho.min_delay_spls(), (float) _scho.max_delay_spls());
        _del_spls.target()[idx] = t;
        f32_x2 v                = (j == 0) ? f32_x2 {_mod[s][1], _mod[s][3]}
                                           : f32_x2 {_mod[s][0], _mod[s][2]};
        auto   mod              = v * v * 1.5f * modmax;
        mod                     = vec_min (modmax, mod);
        _g.scho[idx]            = g + mod;
      }
    }
    for (uint s = (_n_stages * 2); s < _del_spls.target().size(); ++s) {
      _del_spls.target()[s] = _del_spls.target()[s - 1];
    }
    f32_x4 fv {1460.f, 1460.f, 1670.f, 1780.f};
    auto   f = vec_from_array (_mod[_n_stages / 2]);
    f *= f;
    f *= (1.f + s_par.a) * fv;
    _onepole.reset_coeffs<0> (f, _t_spl);
  }
  //----------------------------------------------------------------------------
  template <class T, bool Nested>
  void schroeder_tick (
    xspan<T>       outl,
    xspan<T>       outr,
    xspan<f32_x4>  wet, // might contain the delay feedback
    xspan<T const> inl,
    xspan<T const> inr,
    xspan<std::array<f32_x1, n_channels> const> taps_tail,
    unsmoothed_parameters const&                par,
    xspan<smoothed_parameters const>            s_par)
  {
    uint const block_size = wet.size();

    std::array<decltype (_del_spls)::value_type, max_block_size> del_spls_arr;

    for (uint i = 0; i < block_size; ++i) {
      _del_spls.tick();
      del_spls_arr[i] = _del_spls.get();
      wet[i] += f32_x4 {inl[i], inr[i], inl[i], inr[i]};
      wet[i] = _onepole.tick<0> (wet[i])[0];
    }

    for (uint i = 0; i < block_size; ++i) {
      std::array<f32_x2, max_scho_stages> to_push {};
      wet[i] -= _1spl_fb * s_par[i].feedback;
      std::array<f32_x2, 2> sig {
        {{wet[i][0], wet[i][2]}, {wet[i][1], wet[i][3]}}};
      auto& del_spls = del_spls_arr[i];
      if constexpr (!Nested) {
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
#if 1
        // Nested, as e.g. a lattice
        std::array<f32_x2, 2> fwd {};
        for (uint c = 0; c < n_channels; ++c) {
          uint idx = c;
          auto yn  = _scho.get (del_spls[idx], idx);

          fwd[c] = sig[c] + yn * _g.scho[idx];
          sig[c] = yn - fwd[c] * _g.scho[idx];
        }
        for (uint s = 1; s < _n_stages; ++s) {
          for (uint c = 0; c < n_channels; ++c) {
            uint idx = s * 2 + c;
            auto yn  = _scho.get (del_spls[idx], idx);

            fwd[c] += yn * _g.scho[idx];
            to_push[idx - 2] = yn - fwd[c] * _g.scho[idx];
          }
        }
        to_push[(_n_stages * 2) - 2] = fwd[0];
        to_push[(_n_stages * 2) - 1] = fwd[1];
        sig[0]                       = -sig[0];
        sig[1]                       = -sig[1];
#else
        // Cascaded Nested N=2 allpasses, not that different from the
        // regular allpass...
        std::array<f32_x2, 2> fwd {};
        for (uint s = 0; s < _n_stages; ++s) {
          for (uint c = 0; c < n_channels; ++c) {
            uint idx = s * 2 + c;
            if ((s % 2) == 0) {
              auto yn = _scho.get (del_spls[idx], idx);
              auto g  = _g.scho[idx];
              fwd[c]  = sig[c] + yn * g;
              sig[c]  = yn - fwd[c] * g;
            }
            else {
              auto yn = _scho.get (del_spls[idx], idx);
              auto g  = _g.scho[idx];
              fwd[c] += yn * g;
              to_push[idx - 2] = yn - fwd[c] * g;
              to_push[idx]     = fwd[c];
            }
          }
        }
        if ((_n_stages % 2) != 0) {
          // last stage is not nested
          to_push[(_n_stages * 2) - 2] = fwd[0];
          to_push[(_n_stages * 2) - 1] = fwd[1];
        }
#endif
      }
      _scho.push (to_push);
      wet[i] = vec_shuffle (sig[0], sig[1], 0, 2, 1, 3);

      auto fb_val = zdf_type::nonlin::tick (
        wet[i], vec_set<4> (s_par[i].drive), vec_set<4> (tanh_like_hardness));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);
    }
    // delay_block
    if (par.delay_mode != delay_mode::off) {
      delay_push (s_par, par.delay_mode, taps_tail, wet, inl, inr);
    }
    // final mixing
    for (uint i = 0; i < block_size; ++i) {
      auto out = mix (
        wet[i],
        f64_x2 {inl[i], inr[i]},
        s_par[i].depth / get_fb_gain (s_par[i].feedback, s_par[i].drive),
        1.f - s_par[i].depth,
        1.f);
      outl[i] = out[0];
      outr[i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void chorus_control_block (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par,
    T                            inl,
    T                            inr)
  {
    // stages are in pairs (TODO check)
    float stages
      = 1.f + s_par.stages * ((float) (max_chor_stages / 2) - 1.0000000001f);
    _n_stages      = (uint) stages;
    _n_stages_frac = stages - (float) _n_stages;
    _n_stages += 1; // fractional stage

    float ratefact = s_par.lfo_rate * 0.01f;

    compute_mod (
      0.05f * _rnd_lfo_last[0] + 0.85f - s_par.center, // reverse range lf to hf
      s_par.spread,
      s_par.lfo_depth * (0.45f + 0.2f * s_par.center + 0.25f * ratefact),
      0.3f + 0.1f * _rnd_lfo_last[0] + 0.5f * s_par.center,
      s_par.stereo * 0.5f,
      _n_stages,
      run_mod_srcs (par, s_par, inl, inr));

    for (uint s = 0; s < (_n_stages * 2); ++s) {
      constexpr float msec_offset = 15.f;
      constexpr float sec_factor = ((max_chor_delay_ms - msec_offset) * 0.001f);

      bool  neg = !!((s / 2) % 2);
      float t   = _mod[s][0];
      t         = mod_time_warp (t, s_par.mod_warp);
      t         = neg ? 1.f - t : t;
      t         = _srate * (msec_offset * 0.001f + (t * sec_factor));
      _del_spls.target()[s]    = t;
      constexpr float g_factor = 0.1f;
      auto            m        = _mod[s][1];
      auto            gv       = s_par.b;
      auto            negf     = neg ? 1.f : -1.f;
      _g.chor[s][0]            = gv * gv * 0.407f + m * g_factor * gv;
      _g.chor[s][0] *= neg;
      _rnd_lfo.set_freq (vec_set<4> (0.8f + abs (s_par.b)), _t_spl);

      // set the phasers (no feeedback)
      auto v = vec_from_array (_mod[s]);
      auto f = v * v;
      auto q = v;
      f      = ((f * (0.45f + 0.3f * _rnd_lfo_last)) + 0.15f);
      f *= f32_x4 {21200.f, 21200.f, 17000.f, 17000.f};
      q = (1.f + 2.f * q * abs (s_par.spread));
      q += 0.09f;
      _phaser.reset_coeffs_on_idx (s, f, q, _t_spl);
    }
    for (uint s = (_n_stages * 2); s < _del_spls.target().size(); ++s) {
      _del_spls.target()[s] = _del_spls.target()[s - 1];
    }
    f32_x4 fv {460.f, 460.f, 970.f, 980.f};
    auto   f = vec_from_array (_mod[_n_stages / 2]);
    f *= f;
    f *= (1.f + s_par.b) * fv;
    _onepole.reset_coeffs<0> (f, _t_spl);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void chorus_tick (
    xspan<T>       outl,
    xspan<T>       outr,
    xspan<f32_x4>  wet, // might contain the delay feedback
    xspan<T const> inl,
    xspan<T const> inr,
    xspan<std::array<f32_x1, n_channels> const> taps_tail,
    unsmoothed_parameters const&                par,
    xspan<smoothed_parameters>                  s_par)
  {
    uint const block_size = wet.size();

    std::array<decltype (_del_spls)::value_type, max_block_size> del_spls_arr;
    std::array<decltype (_rnd_lfo_last), max_block_size>         tremolo;
    array2d<
      std::array<float, vec_traits_t<f32_x4>::size>,
      n_channels,
      max_block_size>
                                                               pan;
    array2d<float, vec_traits_t<f32_x4>::size, max_block_size> trnd;

    for (uint i = 0; i < block_size; ++i) {
      wet[i] += f32_x4 {inl[i], inr[i], inl[i], inr[i]};
      wet[i] = _onepole.tick<0> (wet[i])[0];
      _del_spls.tick();
      del_spls_arr[i] = _del_spls.get();
    }

    for (uint i = 0; i < block_size; ++i) {
      // To allow stage crossfading constant panning positions are kept.
      // Processing is done in pairs.
      //
      // pan position   | L             R
      // stage ordering | 1 3 2 4 4 2 3 1
      auto rndlfo   = _rnd_lfo.tick_filt_sample_and_hold();
      _rnd_lfo_last = rndlfo;

      tremolo[i]        = _tremolo_lfo.tick_sine() * s_par[i].b;
      auto tr_crossfade = (tremolo[i] + 1.f) * 0.5f;

      constexpr auto kpanl1 = f32_x4 {0.f, 0.07f, 0.3f, 0.15f};
      constexpr auto kpanl2 = f32_x4 {0.3f, 0.15f, 0.07f, 0.f};

      auto panlv = (1.f - tr_crossfade) * kpanl1;
      panlv += tr_crossfade * kpanl2;

      auto rndmod = (rndlfo * 0.05f * abs (s_par[i].b)) + 0.95f;

      // crossfade towards 0.5 (center panning)
      auto stereo_norm = abs (s_par[i].stereo);
      panlv *= stereo_norm;
      panlv += (1.f - stereo_norm) * 0.5f;
      // using -x^2+2x as a cheap approximation of the sin(x*pi/2) pan
      // law.
      panlv = -panlv * panlv + 2.f * panlv;
      // some randomization
      panlv *= rndmod;
      auto panrv = 1.f - panlv;

      pan[i][0] = vec_to_array (panlv);
      pan[i][1] = vec_to_array (panrv);
      trnd[i]   = vec_to_array (rndlfo * _srate * 0.0005f);
    }

    for (uint i = 0; i < block_size; ++i) {
      auto& del_spls = del_spls_arr[i];
      auto  chor_in  = vec1_array_wrap (vec_to_array (wet[i]));
      wet[i]         = vec_set<4> (0.f);
      auto                                prev = wet[i];
      std::array<f32_x1, max_chor_stages> to_push {};

      for (uint s = 0; s < _n_stages; ++s) {
        uint idx  = s * 2;
        uint lane = idx % vec_traits_t<f32_x4>::size;

        std::array<float, n_channels> channel;

        // TODO: no feedback, delay lines can be read contiguosly

        for (uint j = 0; j < channel.size(); ++j, ++idx, ++lane) {
          auto n_spls = std::clamp (
            del_spls[idx] + trnd[i][lane],
            (float) _chor.min_delay_spls(),
            (float) _chor.max_delay_spls());
          auto g       = _g.chor[idx];
          auto yn      = _chor.get (n_spls, idx);
          auto r       = allpass_fn::tick<f32_x1> (chor_in[lane], yn, g);
          to_push[idx] = r.to_push;
          channel[j]   = r.out[0];
        }
        auto& panl = pan[i][0];
        auto& panr = pan[i][1];
        float l    = channel[0] * panl[s] + channel[1] * panr[s];
        float r    = channel[0] * panr[s] + channel[1] * panl[s];

        // TODO: probably run this on its own loop?
        f32_x4 v {l, r, l, r};
        v    = _phaser.tick_on_idx (s, v);
        prev = wet[i];
        wet[i] += v;
      }
      auto last_stage = wet[i] - prev;
      wet[i]          = prev + last_stage * _n_stages_frac;
      _chor.push (to_push);
    }

    // dry (from flanger) reused as a different mod line
    for (uint i = 0; i < block_size; ++i) {
      f32_x2 dry {inl[i], inr[i]};
      _dry.push (xspan {&dry, 1});
      constexpr auto ratio
        = (float) max_chor_delay_ms / (0.015f + (float) (max_dry_delay_ms));
      float n_spls = del_spls_arr[i][0] * ratio * s_par[i].a;
      n_spls       = std::clamp (
        n_spls, (float) _dry.min_delay_spls(), (float) _dry.max_delay_spls());
      auto ffwd      = _dry.get (n_spls, 0);
      auto fbg       = abs (s_par[i].feedback);
      auto ffwd_gain = fbg * fbg * fbg * fbg * (float) _n_stages;
      ffwd_gain      = std::copysign (ffwd_gain, s_par[i].feedback);
      auto ffwd4     = f32_x4 {ffwd[0], ffwd[1], ffwd[0], ffwd[1]} * ffwd_gain;
      ffwd4 *= 1.f - ((tremolo[i] + 1.f) * 0.5f);
      wet[i] += ffwd4;
      // gain compensation for all gains in the process before entering the
      // saturation stage
      wet[i] /= _n_stages_frac + (float) (_n_stages - 1) + abs (ffwd_gain);
    }

    // saturation + filtering block
    for (uint i = 0; i < block_size; ++i) {
      wet[i] = zdf_type::nonlin::tick (
        wet[i], vec_set<4> (s_par[i].drive), vec_set<4> (tanh_like_hardness));
      wet[i] = _fb_filters.tick_cascade (wet[i]);
    }
    // delay_block
    if (par.delay_mode != delay_mode::off) {
      // this has no feedback, clear the parameter used for other purpuses for
      // "delay_push" (hackish)
      for (uint i = 0; i < block_size; ++i) {
        s_par[i].feedback = 0.f;
      }
      delay_push (s_par, par.delay_mode, taps_tail, wet, inl, inr);
    }
    // dry signal delaying + mixing
    for (uint i = 0; i < block_size; ++i) {
      float depth = (s_par[i].depth - 1.f);
      depth       = 1.f - depth * depth;
      auto out
        = mix (wet[i], f64_x2 {inl[i], inr[i]}, depth, 1.f - depth, s_par[i].b);
      outl[i] = out[0];
      outr[i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void flanger_control_block (
    unsmoothed_parameters const& par,
    smoothed_parameters const&   s_par,
    T                            inl,
    T                            inr)
  {
    compute_mod (
      1.f - s_par.center, // reverse range lf to hf
      s_par.spread,
      s_par.lfo_depth,
      0.45f,
      0.f,
      1, // flanger is always one stage
      run_mod_srcs (par, s_par, inl, inr));

    for (uint s = 0; s < flan_stages; ++s) {
      static_assert (flan_stages == 2, "this loop assumes s < 2");
      constexpr float sec_factor = (max_flan_delay_ms * 0.001f);

      float t = _mod[0][s];
      t       = mod_time_warp (t, s_par.mod_warp);
      t       = _srate * t * sec_factor;
      t       = std::clamp (
        t, (float) _flan.min_delay_spls(), (float) _flan.max_delay_spls());
      _del_spls.target()[s]    = t;
      constexpr float g_factor = 0.1f;
      auto            m        = _mod[0][2 + s];
      auto            gv       = s_par.stages;
      _g.flan[s][0]            = -gv * gv * 0.75f + m * m * g_factor * gv;
    }
    for (uint s = flan_stages; s < _del_spls.target().size(); ++s) {
      _del_spls.target()[s] = _del_spls.target()[s - 1];
    }
    f32_x4 fv {460.f, 460.f, 970.f, 980.f};
    auto   f = vec_from_array (_mod[0]);
    f *= f;
    f *= (1.f + s_par.feedback) * fv;
    _onepole.reset_coeffs<0> (f, _t_spl);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void flanger_tick (
    xspan<T>       outl,
    xspan<T>       outr,
    xspan<f32_x4>  wet, // might contain the delay feedback
    xspan<T const> inl,
    xspan<T const> inr,
    xspan<std::array<f32_x1, n_channels> const> taps_tail,
    unsmoothed_parameters const&                par,
    xspan<smoothed_parameters const>            s_par)
  {
    uint const block_size = wet.size();

    std::array<decltype (_del_spls)::value_type, max_block_size> del_spls_arr;

    for (uint i = 0; i < block_size; ++i) {
      _del_spls.tick();
      del_spls_arr[i] = _del_spls.get();
      auto mid        = (inl[i] + inr[i]) * 0.5f;
      wet[i] += f32_x4 {inl[i], inr[i], mid, mid};
      auto hp_out = _onepole.tick<0> (wet[i])[1];
      wet[i]      = vec_shuffle (wet[i], hp_out, 0, 1, 6, 7);
    }

    for (uint i = 0; i < block_size; ++i) {
      wet[i] -= _1spl_fb * s_par[i].feedback;

      auto& del_spls = del_spls_arr[i];
      auto  flan_in  = vec_split<flan_stages> (wet[i]);

      std::array<f32_x2, flan_stages> to_push {};
      for (uint s = 0; s < flan_stages; ++s) {
        auto yn       = _flan.get (del_spls[s], s);
        auto r        = allpass_fn::tick<f32_x2> (flan_in[s], yn, _g.flan[s]);
        to_push[s]    = r.to_push;
        wet[i][0 + s] = r.out[0];
        wet[i][2 + s] = r.out[1];
      }
      _flan.push (to_push);

      auto fb_val = zdf_type::nonlin::tick (
        wet[i], vec_set<4> (s_par[i].drive), vec_set<4> (tanh_like_hardness));
      _1spl_fb = _fb_filters.tick_cascade (fb_val);
    }
    // delay_block
    if (par.delay_mode != delay_mode::off) {
      delay_push (s_par, par.delay_mode, taps_tail, wet, inl, inr);
    }
    // dry signal delaying + mixing
    for (uint i = 0; i < block_size; ++i) {
      f32_x2 dry {inl[i], inr[i]};
      _dry.push (xspan {&dry, 1});
      constexpr float kdry
        = 0.001f * std::max (max_dry_delay_ms, max_flan_delay_ms / 2);
      float n_spls = _srate * kdry * (1.f - s_par[i].center) * s_par[i].a;
      n_spls       = std::clamp (
        n_spls, (float) _dry.min_delay_spls(), (float) _dry.max_delay_spls());
      dry = _dry.get (n_spls, 0);

      auto out = mix (
        wet[i],
        f64_x2 {dry[0], dry[1]},
        s_par[i].depth / get_fb_gain (s_par[i].feedback, s_par[i].drive),
        1.f - s_par[i].depth,
        s_par[i].b);
      outl[i] = out[0];
      outr[i] = out[1];
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

    float spls_beat = (1.f / _beat_hz) * _srate;
    uint  dly_size  = std::ceil (max_delay_4ths * 4.f * spls_beat);
    dly_size        = _delay.n_required_elems (dly_size + 1, n_channels);

    using T_mem = decltype (_mem)::value_type;
#if !ARTV_MOD_CHO_FLAN_TIRAN
    using T_sinc = decltype (_sinc_co)::value_type;
#endif
    using T_scho  = decltype (_scho)::value_type;
    using T_chor  = decltype (_chor)::value_type;
    using T_flan  = decltype (_flan)::value_type;
    using T_dry   = decltype (_dry)::value_type;
    using T_delay = decltype (_delay)::value_type;

    schroeder_size *= vec_traits_t<T_scho>::size;
    chor_size *= vec_traits_t<T_chor>::size;
    dry_size *= vec_traits_t<T_dry>::size;
    flan_size *= vec_traits_t<T_flan>::size;
    dly_size *= vec_traits_t<T_delay>::size;

    auto max_size
      = std::max (schroeder_size + dly_size, chor_size + dry_size + dly_size);
    max_size = std::max (max_size, flan_size + dry_size + dly_size);

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
    _mem_delay = mem.cut_head (dly_size).cast<T_delay>();
    _mem_scho  = mem.get_head (schroeder_size).cast<T_scho>();
    _mem_dry   = mem.cut_head (dry_size).cast<T_dry>();
    _mem_chor  = mem.get_head (chor_size).cast<T_chor>();
    _mem_flan  = mem.get_head (flan_size).cast<T_flan>();

#if !ARTV_MOD_CHO_FLAN_TIRAN
    // initialize windowed sinc table
    sinc_t::reset_coeffs (_sinc_co, 0.45f, 140.f);
#endif
  }
//----------------------------------------------------------------------------
#if 0
// from linear to  cubic.
  template <class T>
  static T mod_freq_warp (T v, float warp)
  {
    return (1.f - warp) * v + warp * v * v * v;
  }

  //----------------------------------------------------------------------------
  // from cubic to inverse cubic.  0.5 warp = crosses 0.5,0.5 but snakey.
  template <class T>
  static T mod_freq_warp (T v, float warp)
  {
    auto v1   = v * v * v;
    auto vm12 = v - 1.f;
    auto v2   = vm12 * vm12 * vm12 + 1.f;

    return (1.f - warp) * v2 + warp * v1;
  }
    //----------------------------------------------------------------------------
  // from squared to inverse squared. 0.5 warp = linear
  template <class T>
  static T mod_freq_warp (T v, float warp)
  {
    auto v1   = v * v;
    auto vm12 = v - 1.f;
    auto v2   = 1.f - vm12 * vm12;

    return (1.f - warp) * v2 + warp * v1;
  }
    //----------------------------------------------------------------------------
  // from quadractic to inverse squared. 0.5 warp = crosses 0.5,0.5 but snakey.
  template <class T>
  static T mod_freq_warp (T v, float warp)
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
  static T mod_freq_warp (T v, float warp)
  {
    return (1.f - warp) * v * v * v + warp * v;
  }
  //----------------------------------------------------------------------------
  // from linear to cubic. but weighted towards the cubic side
  template <class T>
  static T mod_time_warp (T v, float warp)
  {
    warp *= warp;
    return (1.f - warp) * v * v * v + warp * v;
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   lfo_wave;
    u32   lfo_time_base;
    u32   mode;
    u32   delay_mode;
    float feedback_locut;
    float feedback_hicut;
    float env_speed;
    float env_depth;
    bool  lfo_off; // The lfo is smoothed, this signals if the lfo is off
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float lfo_rate;
    float lfo_depth;
    float mod_warp;
    float stereo;

    float center;
    float a;
    float b;
    float depth; // aka mix (make explicit?)

    float feedback;
    float drive;
    float spread;
    float stages;

    float del_4beats;
    float del_gain;
  };
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  using allpass_type = andy::svf_zdf_allpass;
  using filters_list = mp_list<andy::svf_zdf_highpass, andy::svf_zdf_lowpass>;
  using zdf_type     = zdf::feedback<
    zdf::lin_pre_fb_node_nonlin_after_tag,
    zdf::lin_mystran_tag<2>,
    sigmoid::mystran<1, true>>;

  static constexpr float tanh_like_hardness = 0.19f;
  struct fb_filter {
    enum { locut, hicut, count };
  };

  using sinc_t = sinc_interp<8, 96>;

  struct env {
    enum { fast, slow, smooth1, smooth2, count };
  };
  part_classes<
    mp_list<slew_limiter, slew_limiter, slew_limiter, andy::smoother>,
    sweep_lfo_value>
                                                                    _env;
  part_classes<mp_list<zdf_type>, f32_x4>                           _feedback;
  part_class_array<allpass_type, f32_x4, max_phaser_stages>         _phaser;
  part_classes<filters_list, f32_x4>                                _fb_filters;
  part_classes<mp_list<onepole<allpass_tag, highpass_tag>>, f32_x4> _onepole;
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
  modulable_thiran1_delay_line<f32_x1, 4, false, false> _delay;
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
  sweep_lfo_type                  _sweep_lfo; // 0 = L, 1 = R, control rate
  lfo<4>                          _rnd_lfo; // 0 = L, 1 = R, audio rate
  lfo<4>                          _tremolo_lfo; // 0 = L, 1 = R, audio rate
  uint                            _n_processed_samples;
  uint                            _control_rate_mask;
  uint                            _n_stages;
  f32_x4                          _1spl_fb;
  decltype (_rnd_lfo)::value_type _rnd_lfo_last;
  float                           _lp_smooth_coeff;
  float                           _n_stages_frac;
  float                           _beat_hz;
  float                           _1_4beat_spls;
  float                           _ctrl_t_spl;
  float                           _t_spl;
  float                           _srate;
#if !ARTV_MOD_CHO_FLAN_TIRAN
  xspan<float> _sinc_co;
#endif
  xspan<f32_x1> _mem_delay;
  xspan<f32_x2> _mem_dry;
  xspan<f32_x1> _mem_chor;
  xspan<f32_x2> _mem_flan;
  xspan<f32_x2> _mem_scho;
};
//------------------------------------------------------------------------------
} // namespace artv
