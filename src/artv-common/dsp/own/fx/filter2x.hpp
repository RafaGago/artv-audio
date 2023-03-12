#pragma once

#include <algorithm>

#include <bl/base/static_integer_math.h>

#include "artv-common/dsp/own/classes/control_rate.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/presence.hpp"
#include "artv-common/dsp/own/parts/filters/saike.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
class filter2x {
private:
  enum class paramtype {
    band_type,
    topology,
    frequency,
    reso,
    drive,
    gain,
    tolerance,
    feedback,
    ef_to_frequency,
    ef_to_resonance
  };

public:
  template <uint Band, paramtype Type>
  struct param {
    static constexpr uint band  = Band;
    static constexpr uint ptype = (uint) Type;
  };
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::eq;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::band_type>, int v)
  {
    static_assert (band < n_bands, "");

    if (v >= (int) bandtype::size || v < 0) {
      v = 0;
    }
    _cfg[band].has_changes |= (int) _cfg[band].type != (int) v;
    _cfg[band].reset_band_state |= _cfg[band].has_changes;
    _cfg[band].type = (u8) v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::band_type>)
  {
    return choice_param (
      band == 0 ? 1 : 0,
      make_cstr_array (
        "Off",
        "K35",
        "K35 Asym",
        "Steiner",
        "Steiner Asym",
        "Ladder 4pole",
        "Ladder 2pole"),
      20);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::topology>, int v)
  {
    static_assert (band < n_bands, "");
    if (v >= n_topologies || v < 0) {
      v = (int) bandtype::off;
    }
    _cfg[band].has_changes |= (int) _cfg[band].topology != (int) v;
    _cfg[band].reset_band_state |= _cfg[band].has_changes;
    _cfg[band].topology = (u8) v;
  }

  static constexpr uint n_topologies = 4;

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::topology>)
  {
    return choice_param (
      0, make_cstr_array ("Lowpass", "Highpass", "Bandpass", "Bandreject"), 8);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::frequency>, float v)
  {
    static_assert (band < n_bands, "");
    v = midi_note_to_hz (v);
    _cfg[band].has_changes |= _cfg[band].freq != v;
    _cfg[band].freq = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::frequency>)
  {
    return frequency_parameter (20.0, 20000.0, 440.0);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::reso>, float v)
  {
    static_assert (band < n_bands, "");
    auto& b = _cfg[band];
    b.has_changes |= b.reso != v;
    _cfg[band].reso = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::reso>)
  {
    return float_param ("", 0.f, 1.f, 0.5f, 0.0001f);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::drive>, float v)
  {
    static_assert (band < n_bands, "");
    _cfg[band].has_changes |= _cfg[band].drive_db != v;
    _cfg[band].drive_db = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::drive>)
  {
    return float_param ("dB", -40.f, 30.f, 0.f, 0.3f, 0.6f, true);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::gain>, float v)
  {
    static_assert (band < n_bands, "");
    _cfg[band].has_changes |= _cfg[band].gain_db != v;
    _cfg[band].gain_db = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::gain>)
  {
    return float_param ("dB", -30.f, 20.f, 0.f, 0.3f, 0.6f, true);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::tolerance>, float v)
  {
    static_assert (band < n_bands, "");
    v *= 0.01;
    _cfg[band].has_changes |= _cfg[band].tolerance != v;
    _cfg[band].tolerance = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::tolerance>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::feedback>, float v)
  {
    static_assert (band < n_bands, "");
    v *= (0.01f * 0.93f);
    _cfg[band].feedback = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::feedback>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  struct envfollow_attack_tag {};

  void set (envfollow_attack_tag, float v)
  {
    v *= 0.001f; // to seconds
    if (v != _ef.attack) {
      _ef.attack = v;
      update_envelope_follower();
    }
  }

  static constexpr auto get_parameter (envfollow_attack_tag)
  {
    return float_param ("ms", 1, 180., 20., 0.1);
  }
  //----------------------------------------------------------------------------
  struct envfollow_sensitivity_tag {};

  void set (envfollow_sensitivity_tag, float v) { _ef.gain = db_to_gain (v); }

  static constexpr auto get_parameter (envfollow_sensitivity_tag)
  {
    return float_param ("dB", -10., 30., 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct envfollow_release_tag {};

  void set (envfollow_release_tag, float v)
  {
    v *= 0.001f; // to seconds
    if (v != _ef.release) {
      _ef.release = v;
      update_envelope_follower();
    }
  }

  static constexpr auto get_parameter (envfollow_release_tag)
  {
    return float_param ("ms", 40., 1000., 150., 1.);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::ef_to_frequency>, float v)
  {
    static_assert (band < n_bands, "");
    v *= (0.01f * 4.f);
    _cfg[band].ef_to_freq = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::ef_to_frequency>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::ef_to_resonance>, float v)
  {
    static_assert (band < n_bands, "");
    v *= (0.01f * 0.5f);
    _cfg[band].ef_to_reso = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::ef_to_resonance>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  using band1_type_tag       = param<0, paramtype::band_type>;
  using band2_type_tag       = param<1, paramtype::band_type>;
  using band1_topology_tag   = param<0, paramtype::topology>;
  using band2_topology_tag   = param<1, paramtype::topology>;
  using band1_freq_tag       = param<0, paramtype::frequency>;
  using band2_freq_tag       = param<1, paramtype::frequency>;
  using band1_reso_tag       = param<0, paramtype::reso>;
  using band2_reso_tag       = param<1, paramtype::reso>;
  using band1_drive_tag      = param<0, paramtype::drive>;
  using band2_drive_tag      = param<1, paramtype::drive>;
  using band1_gain_tag       = param<0, paramtype::gain>;
  using band2_gain_tag       = param<1, paramtype::gain>;
  using band1_tolerance_tag  = param<0, paramtype::tolerance>;
  using band2_tolerance_tag  = param<1, paramtype::tolerance>;
  using band1_feedback_tag   = param<0, paramtype::feedback>;
  using band2_feedback_tag   = param<1, paramtype::feedback>;
  using band1_ef_to_freq_tag = param<0, paramtype::ef_to_frequency>;
  using band2_ef_to_freq_tag = param<1, paramtype::ef_to_frequency>;
  using band1_ef_to_reso_tag = param<0, paramtype::ef_to_resonance>;
  using band2_ef_to_reso_tag = param<1, paramtype::ef_to_resonance>;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    band1_type_tag,
    band1_topology_tag,
    band1_freq_tag,
    band1_drive_tag,
    band1_reso_tag,
    band1_gain_tag,
    band1_tolerance_tag,
    band1_feedback_tag,
    band2_type_tag,
    band2_topology_tag,
    band2_freq_tag,
    band2_drive_tag,
    band2_reso_tag,
    band2_gain_tag,
    band2_tolerance_tag,
    band2_feedback_tag,
    envfollow_attack_tag,
    envfollow_release_tag,
    envfollow_sensitivity_tag,
    band1_ef_to_freq_tag,
    band1_ef_to_reso_tag,
    band2_ef_to_freq_tag,
    band2_ef_to_reso_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    using x1_t = vec<double, 1>;

    _t_spl = 1.f / pc.get_sample_rate();
    _cfg   = decltype (_cfg) {};

    smoother::reset_coeffs (
      xspan {&_smooth_coeff, 1}.cast (x1_t {}),
      vec_set<x1_t> (1. / 0.08),
      _t_spl);

    memset (&_smooth_pars, 0, sizeof _smooth_pars);
    memset (&_last_sample, 0, sizeof _last_sample);

    _filter.zero_all_states();
    _proc.reset_states<proc_envfollow>();
    _proc.reset_states<proc_dc_block>();
    _proc.reset_coeffs<proc_dc_block> (vec_set<f64_x2> (1.), _t_spl);

    _ef_value = vec_set<f64_x2> (0.);

    constexpr uint base_order = bl_static_log2_floor_u64 ((u64) blocksize);
    uint sr_order = get_samplerate_order (pc.get_sample_rate()) + base_order;
    _control.reset (1 << (sr_order + 2));
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint block_smp = 0; block_smp < samples; block_smp += blocksize) {
      uint n_samples = std::min<uint> (blocksize, samples - block_smp);

      // both filters are parallel/operate on the same input
      std::array<f64_x2, blocksize> in_cp;
      for (uint i = 0; i < blocksize; ++i) {
        in_cp[i] = f64_x2 {ins[0][block_smp + i], ins[1][block_smp + i]};
      }
      memset (&outs[0][block_smp], 0, n_samples * sizeof (T));
      memset (&outs[1][block_smp], 0, n_samples * sizeof (T));

      int  control_offset   = _control.tick (n_samples);
      bool is_control_block = control_offset >= 0;

      // control block
      if (is_control_block) {
        // Envelope follower/modulation
        auto prev = _ef_value;
        _ef_value = _proc.tick<proc_envfollow> (in_cp[control_offset]);
        // Make unipolar and clip at 1.
        _ef_value = vec_min (vec_abs (_ef_value * _ef.gain), 1.);
        _ef_value = (_ef_value + prev) * 0.5; // Box LP
        _ef_value *= _ef_value;
      }

      for (uint b = 0; b < n_bands; ++b) {
        if (_cfg[b].has_changes || is_control_block) {
          reset_band (b);
        }

        auto btype = get_band_type (b);
        if (btype == bandtype::off) {
          continue;
        }

        alignas (sse_bytes) smoothed_params target_pars;

        for (uint i = 0; i < n_samples; ++i) {
          // parameter update block
          static_assert (sizeof _smooth_pars[b] / sizeof (float) == 4, "");

          set_target_params (target_pars, b);
          auto smcoeff = (float) _smooth_coeff;

          smoother::tick<f32_x4> (
            xspan {&smcoeff, 1},
            xspan {_smooth_pars[b].arr}.cast (f32_x4 {}),
            vec_load<f32_x4> (target_pars.arr.data()));
        }

        xspan<f64_x2> internal = _filter.get_coeffs (b);
        for (uint j = 0; j < _target_coeffs[b].size(); ++j) {
          for (uint i = 0; i < n_samples; ++i) {
            internal[j] = smoother::tick (
              xspan {&_smooth_coeff, 1},
              xspan {&internal[j], 1},
              _target_coeffs[b][j]);
          }
        }

        auto& smoothed_p = _smooth_pars[b];

        f64_x2 fb_gain = vec_set<f64_x2> (smoothed_p.vars.feedback);
        if (is_moog_1 (btype)) {
          fb_gain *= 0.03; // moogs are extremely unstable reduce feedback
        }
        else if (is_moog_2 (btype)) {
          fb_gain *= 0.03; // moogs are extremely unstable reduce feedback
        }

        f64_x2 pre_drive
          = f64_x2 {smoothed_p.vars.pre_drive_l, smoothed_p.vars.pre_drive_r};

        std::array<f64_x2, blocksize> out;
        f64_x2                        last = _last_sample[b];

        switch (btype) {
        case bandtype::off: // not reachable
          jassert (false);
          break;
        case bandtype::ms20_lp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_lowpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_hp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_highpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_bp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_bandpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_notch> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_asym_lp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_asym_lowpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_asym_hp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_asym_highpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_asym_bp:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_asym_bandpass> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::ms20_asym_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::ms20_asym_notch> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::steiner_1_lp:
        case bandtype::steiner_1_hp:
        case bandtype::steiner_1_bp:
        case bandtype::steiner_1_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::steiner_1> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::steiner_2_lp:
        case bandtype::steiner_2_hp:
        case bandtype::steiner_2_bp:
        case bandtype::steiner_2_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::steiner_2> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::moog_1_lp:
        case bandtype::moog_1_hp:
        case bandtype::moog_1_bp:
        case bandtype::moog_1_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::moog_1> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        case bandtype::moog_2_lp:
        case bandtype::moog_2_hp:
        case bandtype::moog_2_bp:
        case bandtype::moog_2_br:
          for (uint i = 0; i < n_samples; ++i) {
            auto in = (in_cp[i] * pre_drive) + (last * fb_gain);
            out[i]  = _filter.tick_on_idx<saike::moog_2> (b, in);
            last    = _proc.tick<proc_dc_block> (out[i]);
          }
          break;
        default:
          jassert (false);
          break;
        }
        _last_sample[b] = last;

        for (uint i = 0; i < n_samples; ++i) {
          static constexpr double limit = constexpr_db_to_gain (10.);

          out[i] *= smoothed_p.vars.post_drive;
          vec_clamp (out[i], -limit, limit);

          outs[0][block_smp + i] += out[i][0];
          outs[1][block_smp + i] += out[i][1];
        }
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  enum class bandtype {
    off,
    ms20_beg,
    ms20_lp = ms20_beg,
    ms20_hp,
    ms20_bp,
    ms20_br,
    ms20_end,
    ms20_asym_beg = ms20_end,
    ms20_asym_lp  = ms20_asym_beg,
    ms20_asym_hp,
    ms20_asym_bp,
    ms20_asym_br,
    ms20_asym_end,
    steiner_1_beg = ms20_asym_end,
    steiner_1_lp  = steiner_1_beg,
    steiner_1_hp,
    steiner_1_bp,
    steiner_1_br,
    steiner_1_end,
    steiner_2_beg = steiner_1_end,
    steiner_2_lp  = steiner_2_beg,
    steiner_2_hp,
    steiner_2_bp,
    steiner_2_br,
    steiner_2_end,
    moog_1_beg = steiner_2_end,
    moog_1_lp  = moog_1_beg,
    moog_1_hp,
    moog_1_bp,
    moog_1_br,
    moog_1_end,
    moog_2_beg = moog_1_end,
    moog_2_lp  = moog_2_beg,
    moog_2_hp,
    moog_2_bp,
    moog_2_br,
    moog_2_end,
    size = moog_2_end
  };
  //----------------------------------------------------------------------------
  static constexpr bool is_ms20 (bandtype type)
  {
    return type >= bandtype::ms20_beg && type < bandtype::ms20_end;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_ms20_asym (bandtype type)
  {
    return type >= bandtype::ms20_asym_beg && type < bandtype::ms20_asym_end;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_steiner_1 (bandtype type)
  {
    return type >= bandtype::steiner_1_beg && type < bandtype::steiner_1_end;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_steiner_2 (bandtype type)
  {
    return type >= bandtype::steiner_2_beg && type < bandtype::steiner_2_end;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_moog_1 (bandtype type)
  {
    return type >= bandtype::moog_1_beg && type < bandtype::moog_1_end;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_moog_2 (bandtype type)
  {
    return type >= bandtype::moog_2_beg && type < bandtype::moog_2_end;
  }
  //----------------------------------------------------------------------------
  bandtype get_band_type (uint band)
  {
    uint offset = (_cfg[band].type * n_topologies) - (n_topologies - 1);
    uint t      = _cfg[band].type == 0 ? 0 : offset + _cfg[band].topology;
    return (bandtype) t;
  }
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    constexpr double min_hz
      = constexpr_midi_note_to_hz (get_parameter (band1_freq_tag {}).min);
    constexpr double max_hz
      = constexpr_midi_note_to_hz (get_parameter (band1_freq_tag {}).max);
    constexpr double min_reso = get_parameter (band1_reso_tag {}).min;
    constexpr double max_reso = get_parameter (band1_reso_tag {}).max;

    auto& b = _cfg[band];

    auto ef_mod = f64_x2 {};
    if (b.ef_to_freq != 0.f || b.ef_to_reso != 0.f) {
      ef_mod = ((_ef_value * _ef_value) + vec_sqrt (_ef_value)) * 0.5;
    }

    f64_x2 freq = {b.freq, b.freq * (1.f + (b.tolerance * 0.2f))};
    freq *= vec_exp2 (ef_mod * b.ef_to_freq);
    freq = vec_clamp (freq, min_hz, max_hz);

    f64_x2 reso = {b.reso, b.reso * (1.f - (b.tolerance * 0.11f))};
    reso += reso * ef_mod * b.ef_to_reso;
    reso = vec_clamp (reso, min_reso, max_reso);

    auto btype = get_band_type (band);

    auto co = xspan {_target_coeffs[band]};

    switch (btype) {
    case bandtype::off:
      _filter.reset_states_on_idx (band);
      break;
    case bandtype::ms20_lp:
      _filter.reset_coeffs_ext<saike::ms20_lowpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_hp:
      _filter.reset_coeffs_ext<saike::ms20_highpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_bp:
      _filter.reset_coeffs_ext<saike::ms20_bandpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_br:
      _filter.reset_coeffs_ext<saike::ms20_notch> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_asym_lp:
      _filter.reset_coeffs_ext<saike::ms20_asym_lowpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_asym_hp:
      _filter.reset_coeffs_ext<saike::ms20_asym_highpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_asym_bp:
      _filter.reset_coeffs_ext<saike::ms20_asym_bandpass> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::ms20_asym_br:
      _filter.reset_coeffs_ext<saike::ms20_asym_notch> (
        band, co, freq, reso, _t_spl);
      break;
    case bandtype::steiner_1_lp:
      _filter.reset_coeffs_ext<saike::steiner_1> (
        band, co, freq, reso, _t_spl, lowpass_tag {});
      break;
    case bandtype::steiner_1_hp:
      _filter.reset_coeffs_ext<saike::steiner_1> (
        band, co, freq, reso, _t_spl, highpass_tag {});
      break;
    case bandtype::steiner_1_bp:
      _filter.reset_coeffs_ext<saike::steiner_1> (
        band, co, freq, reso, _t_spl, bandpass_tag {});
      break;
    case bandtype::steiner_1_br:
      _filter.reset_coeffs_ext<saike::steiner_1> (
        band, co, freq, reso, _t_spl, notch_tag {});
      break;
    case bandtype::steiner_2_lp:
      _filter.reset_coeffs_ext<saike::steiner_2> (
        band, co, freq, reso, _t_spl, lowpass_tag {});
      break;
    case bandtype::steiner_2_hp:
      _filter.reset_coeffs_ext<saike::steiner_2> (
        band, co, freq, reso, _t_spl, highpass_tag {});
      break;
    case bandtype::steiner_2_bp:
      _filter.reset_coeffs_ext<saike::steiner_2> (
        band, co, freq, reso, _t_spl, bandpass_tag {});
      break;
    case bandtype::steiner_2_br:
      _filter.reset_coeffs_ext<saike::steiner_2> (
        band, co, freq, reso, _t_spl, notch_tag {});
      break;
    case bandtype::moog_1_lp:
      _filter.reset_coeffs_ext<saike::moog_1> (
        band, co, freq, reso, _t_spl, lowpass_tag {});
      break;
    case bandtype::moog_1_hp:
      _filter.reset_coeffs_ext<saike::moog_1> (
        band, co, freq, reso, _t_spl, highpass_tag {});
      break;
    case bandtype::moog_1_bp:
      _filter.reset_coeffs_ext<saike::moog_1> (
        band, co, freq, reso, _t_spl, bandpass_tag {});
      break;
    case bandtype::moog_1_br:
      _filter.reset_coeffs_ext<saike::moog_1> (
        band, co, freq, reso, _t_spl, notch_tag {});
      break;
    case bandtype::moog_2_lp:
      _filter.reset_coeffs_ext<saike::moog_2> (
        band, co, freq, reso, _t_spl, lowpass_tag {});
      break;
    case bandtype::moog_2_hp:
      _filter.reset_coeffs_ext<saike::moog_2> (
        band, co, freq, reso, _t_spl, highpass_tag {});
      break;
    case bandtype::moog_2_bp:
      _filter.reset_coeffs_ext<saike::moog_2> (
        band, co, freq, reso, _t_spl, bandpass_tag {});
      break;
    case bandtype::moog_2_br:
      _filter.reset_coeffs_ext<saike::moog_2> (
        band, co, freq, reso, _t_spl, notch_tag {});
      break;
    default:
      jassert (false);
    }
    // drive handling
    // The idea is to hit the filter channel at 0dBVU and that the drive has
    // an useful range at that level. I adjust with a Q before self-oscillation
    // and try to slightly correct very wild gain differences. I used a loop
    // hoovering on trakmeter's 0dB RMS and Reaper's white noise get at -18dB

    static constexpr float moog_gain    = constexpr_db_to_gain (-28.);
    static constexpr float steiner_gain = constexpr_db_to_gain (15.);

    float correction = is_moog_1 (btype) || is_moog_2 (btype) ? moog_gain : 1.;
    float pre_correction
      = is_steiner_1 (btype) || is_steiner_1 (btype) ? steiner_gain : 1.;

    float drive            = db_to_gain (b.drive_db);
    _cfg[band].pre_drive_l = drive;
    _cfg[band].pre_drive_l *= correction;
    _cfg[band].pre_drive_l *= pre_correction;

    // wanting Drive to act as a weak gain on the positive side too so
    // applying on the numerator:
    // (((x^2) / 64 + (63/64), which crosses at (x=1, y=1).
    //
    // Instead of just applying (1. / pre_drive).
    _cfg[band].post_drive
      = (((drive * drive) * 0.015625) + 0.984375) / _cfg[band].pre_drive_l;
    _cfg[band].post_drive *= correction;

    _cfg[band].post_drive *= db_to_gain (_cfg[band].gain_db);
    _cfg[band].pre_drive_r
      = _cfg[band].pre_drive_l * (1.f - (b.tolerance * 0.13f));

    // reset smoothing
    if (_cfg[band].reset_band_state) {
      _filter.reset_states_on_idx (band);
      xspan_copy<f64_x2> (_filter.get_coeffs (band), co);
      set_target_params (_smooth_pars[band], band);
    }
    _cfg[band].has_changes      = false;
    _cfg[band].reset_band_state = false;
  }
  //----------------------------------------------------------------------------
  void update_envelope_follower()
  {
    _proc.reset_coeffs<proc_envfollow> (
      vec_set<f64_x2> (_ef.attack),
      vec_set<f64_x2> (_ef.release),
      _t_spl * (float) _control.get_period());
  }
  //----------------------------------------------------------------------------
  static constexpr uint blocksize = 32;
  //----------------------------------------------------------------------------
  enum ef_type {
    ef_sqrt,
    ef_linear,
    ef_pow2,
    ef_pow3,
  };
  //----------------------------------------------------------------------------
  static constexpr uint n_bands    = 2;
  static constexpr uint n_channels = 2;
  //----------------------------------------------------------------------------
  using filters = parts_union_array<
    mp_list<
      saike::ms20_lowpass,
      saike::ms20_highpass,
      saike::ms20_bandpass,
      saike::ms20_notch,
      saike::ms20_asym_lowpass,
      saike::ms20_asym_highpass,
      saike::ms20_asym_bandpass,
      saike::ms20_asym_notch,
      saike::steiner_1,
      saike::steiner_2,
      saike::moog_1,
      saike::moog_2>,
    f64_x2,
    n_bands>;
  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  struct bandconfig {
    float freq             = 440.f;
    float reso             = 0.7f;
    float drive_db         = 0.f;
    float gain_db          = 0.f;
    float tolerance        = 0.f;
    float dry_wet          = 1.f;
    float pre_drive_l      = 1.f;
    float pre_drive_r      = 1.f;
    float post_drive       = 1.f;
    float feedback         = 0.f;
    float ef_to_freq       = 0.f;
    float ef_to_reso       = 0.f;
    u8    type             = 0;
    u8    topology         = 0;
    bool  has_changes      = false;
    bool  reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  union smoothed_params {
    struct {
      float pre_drive_l, pre_drive_r, post_drive, feedback;
    } vars;
    std::array<float, 4> arr;
  };
  //----------------------------------------------------------------------------
  void set_target_params (smoothed_params& dst, uint band)
  {
    dst.vars.pre_drive_l = _cfg[band].pre_drive_l;
    dst.vars.pre_drive_r = _cfg[band].pre_drive_r;
    dst.vars.post_drive  = _cfg[band].post_drive;
    dst.vars.feedback    = _cfg[band].feedback;
  }
  //----------------------------------------------------------------------------
  struct ef_cfg {
    float attack  = 0;
    float release = 0.;
    // uint  mode    = 0;
    float gain = 0.;
  };
  //----------------------------------------------------------------------------
  std::array<bandconfig, n_bands> _cfg;

  alignas (sse_bytes) std::array<smoothed_params, n_bands> _smooth_pars;

  using coeff_array = std::array<filters::value_type, filters::n_coeffs>;

  alignas (sse_bytes) std::array<coeff_array, n_bands> _target_coeffs;
  filters _filter;
  enum { proc_envfollow, proc_dc_block };
  part_classes<mp_list<slew_limiter, mystran_dc_blocker>, f64_x2> _proc;

  alignas (sse_bytes) std::array<f64_x2, n_bands> _last_sample;
  f64_x2       _ef_value;
  ef_cfg       _ef;
  control_rate _control;

  double _smooth_coeff;

  float _t_spl;
};
//------------------------------------------------------------------------------
} // namespace artv
