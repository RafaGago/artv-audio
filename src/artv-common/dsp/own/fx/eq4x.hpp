#pragma once

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/presence.hpp"
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
class eq4x {
private:
  enum class paramtype { band_type, frequency, q, gain, diff };

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
  void reset (plugin_context& pc)
  {
    using x1_t = vec<decltype (_smooth_coeff), 1>;
    _t_spl     = 1.f / pc.get_sample_rate();
    _cfg       = decltype (_cfg) {};
    smoother::reset_coeffs (
      xspan {&_smooth_coeff, 1}.cast (x1_t {}),
      vec_set<x1_t> (1. / 0.02),
      _t_spl);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    uint n_enabled_bands = 0;

    std::array<bool, 2>    is_copied {false, false};
    std::array<f64_x2, 32> io;

    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].has_changes) {
        reset_band (b);
      }
      if (_cfg[b].type == bandtype::off) {
        continue;
      }

      ++n_enabled_bands;

      for (uint offset = 0; offset < samples; offset += io.size()) {
        uint blocksize = std::min<uint> (io.size(), samples - offset);
        // interleaving

        std::array<T const*, 2> src;
        src[0] = is_copied[0] ? &outs[0][offset] : &ins[0][offset];
        src[1] = is_copied[1] ? &outs[1][offset] : &ins[1][offset];

        for (uint i = 0; i < blocksize; ++i) {
          io[i][0] = src[0][i];
          io[i][1] = src[1][i];
        }

        // smoothing (at blocksize rate)
        xspan<f64_x2> internal = _eq.get_coeffs (b);
        for (uint j = 0; j < _target_coeffs[b].size(); ++j) {
          for (uint i = 0; i < blocksize; ++i) {
            internal[j] = smoother::tick (
              xspan {&_smooth_coeff, 1},
              xspan {&internal[j], 1},
              _target_coeffs[b][j]);
          }
        }
        // processing
        switch (_cfg[b].type) {
        case bandtype::off: // not reachable
        case bandtype::svf_bell:
        case bandtype::svf_lshelf:
        case bandtype::svf_hshelf:
        case bandtype::svf_allpass:
        case bandtype::svf_bell_bandpass:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<andy::svf> (b, io[i]);
          }
          break;
        case bandtype::butterworth_lp:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<btw_lp> (
              b, io[i], q_to_butterworth_order (_cfg[b].q));
          }
          break;
        case bandtype::butterworth_hp:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<btw_hp> (
              b, io[i], q_to_butterworth_order (_cfg[b].q));
          }
          break;
        case bandtype::svf_tilt:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<tilt_eq> (b, io[i]);
          }
          break;
        case bandtype::presence:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<liteon::presence_high_shelf> (b, io[i]);
          }
          break;
        case bandtype::onepole_allpass:
          for (uint i = 0; i < blocksize; ++i) {
            io[i] = _eq.tick_on_idx<onepole_allpass> (b, io[i]);
          }
          break;
        default:
          jassert (false);
          break;
        }
        // deinterleaving
        switch (_topology) {
        case topology::stereo:
          for (uint i = 0; i < blocksize; ++i) {
            outs[0][offset + i] = io[i][0];
            outs[1][offset + i] = io[i][1];
          }
          break;
        case topology::l:
          for (uint i = 0; i < blocksize; ++i) {
            outs[0][offset + i] = io[i][0];
          }
          break;
        case topology::r:
          for (uint i = 0; i < blocksize; ++i) {
            outs[1][offset + i] = io[i][1];
          }
          break;
        case topology::lr_half: {
          uint chnl = b & 1;
          for (uint i = 0; i < blocksize; ++i) {
            outs[chnl][offset + i] = io[i][chnl];
          }
        } break;
        default:
          assert (false);
          break;
        }
      }

      switch (_topology) {
      case topology::stereo:
        is_copied[0] = is_copied[1] = true;
        break;
      case topology::l:
        is_copied[0] = true;
        break;
      case topology::r:
        is_copied[1] = true;
        break;
      case topology::lr_half:
        is_copied[b & 1] = true;
        break;
      default:
        break;
      }
    }
    for (uint c = 0; c < 2; ++c) {
      if (unlikely (!is_copied[c])) {
        memcpy (outs[c], ins[c], samples * sizeof outs[c][0]);
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::band_type>, int v)
  {
    static_assert (band < n_bands, "");
    if (v >= (int) bandtype::size || v < 0) {
      v = (int) bandtype::off;
    }
    if ((int) _cfg[band].type == (int) v) {
      return;
    }
    _cfg[band].has_changes      = true;
    _cfg[band].reset_band_state = true;
    _cfg[band].type             = (bandtype) v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::band_type>)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Off",
        "Peak",
        "LowShelf",
        "HighShelf",
        "Allpass",
        "Butterworth LP",
        "Butterworth HP",
        "Tilt",
        "Presence HiShelf",
        "Allpass 1-pole",
        "Peak BP"),
      30);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::frequency>, float v)
  {
    static_assert (band < n_bands, "");
    _cfg[band].has_changes |= _cfg[band].freq_note != v;
    _cfg[band].freq_note = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::frequency>)
  {
    return frequency_parameter (20.0, 20000.0, 440.0);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::q>, float v)
  {
    static_assert (band < n_bands, "");
    auto& b = _cfg[band];
    // changing butterworth slope requires a state reset, as they are
    // different filters under the same menu.
    bool change = false;
    bool reset  = (b.type == bandtype::butterworth_hp);
    reset |= (b.type == bandtype::butterworth_lp);
    if (reset) {
      reset &= (q_to_butterworth_order (b.q) != q_to_butterworth_order (v));
    }
    else {
      change = _cfg[band].q != v;
    }
    b.reset_band_state |= reset;
    b.has_changes |= reset | change;
    _cfg[band].q = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::q>)
  {
    return float_param ("", 0.1f, 20.f, 0.5f, 0.0001f, 0.3f);
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
    return float_param ("dB", -20.f, 20.f, 0.f, 0.25f, 0.6f, true);
  }
  //----------------------------------------------------------------------------
  template <uint band>
  void set (param<band, paramtype::diff>, float v)
  {
    static_assert (band < n_bands, "");
    v *= 0.01 * 2.; // 2 octaves up and down
    _cfg[band].has_changes |= _cfg[band].diff != v;
    _cfg[band].diff = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::diff>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.05f);
  }
  //----------------------------------------------------------------------------
  struct topology_tag {};
  void set (topology_tag, int v)
  {
    auto t = (topology) v;
    if (t < topology::count) {
      _topology = t;
    }
  }

  static constexpr auto get_parameter (topology_tag)
  {
    return choice_param (
      0, make_cstr_array ("Stereo", "L", "R", "L:1-3, R:2-4"), 16);
  }
  //----------------------------------------------------------------------------
  using band1_type_tag = param<0, paramtype::band_type>;
  using band2_type_tag = param<1, paramtype::band_type>;
  using band3_type_tag = param<2, paramtype::band_type>;
  using band4_type_tag = param<3, paramtype::band_type>;
  using band1_freq_tag = param<0, paramtype::frequency>;
  using band2_freq_tag = param<1, paramtype::frequency>;
  using band3_freq_tag = param<2, paramtype::frequency>;
  using band4_freq_tag = param<3, paramtype::frequency>;
  using band1_q_tag    = param<0, paramtype::q>;
  using band2_q_tag    = param<1, paramtype::q>;
  using band3_q_tag    = param<2, paramtype::q>;
  using band4_q_tag    = param<3, paramtype::q>;
  using band1_gain_tag = param<0, paramtype::gain>;
  using band2_gain_tag = param<1, paramtype::gain>;
  using band3_gain_tag = param<2, paramtype::gain>;
  using band4_gain_tag = param<3, paramtype::gain>;
  using band1_diff_tag = param<0, paramtype::diff>;
  using band2_diff_tag = param<1, paramtype::diff>;
  using band3_diff_tag = param<2, paramtype::diff>;
  using band4_diff_tag = param<3, paramtype::diff>;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    band1_type_tag,
    band1_freq_tag,
    band1_gain_tag,
    band1_q_tag,
    band2_type_tag,
    band2_freq_tag,
    band2_gain_tag,
    band2_q_tag,
    band3_type_tag,
    band3_freq_tag,
    band3_gain_tag,
    band3_q_tag,
    band4_type_tag,
    band4_freq_tag,
    band4_gain_tag,
    band4_q_tag,
    band1_diff_tag,
    band2_diff_tag,
    band3_diff_tag,
    band4_diff_tag,
    topology_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    using x1_t = vec<double, 1>;

    auto& b             = _cfg[band];
    auto  bandtype_prev = b.type;

    auto freq = vec_set<f64_x2> (midi_note_to_hz (b.freq_note));
    freq[1] *= exp2 (b.diff);
    auto q = vec_set<f64_x2> (b.q);
    q[1] += q[1] * b.diff * 0.05;
    auto gain = vec_set<f64_x2> (b.gain_db);

    switch (b.type) {
    case bandtype::off:
      _eq.reset_states_on_idx (band);
      break;
    case bandtype::svf_bell:
      _eq.reset_coeffs_ext<andy::svf> (
        band, _target_coeffs[band], freq, q, gain, _t_spl, bell_tag {});
      break;
    case bandtype::svf_lshelf:
      _eq.reset_coeffs_ext<andy::svf> (
        band, _target_coeffs[band], freq, q, gain, _t_spl, lowshelf_tag {});
      break;
    case bandtype::svf_hshelf:
      _eq.reset_coeffs_ext<andy::svf> (
        band, _target_coeffs[band], freq, q, gain, _t_spl, highshelf_tag {});
      break;
    case bandtype::svf_allpass:
      _eq.reset_coeffs_ext<andy::svf> (
        band, _target_coeffs[band], freq, q, _t_spl, allpass_tag {});
      break;
    case bandtype::butterworth_lp:
      _eq.reset_coeffs_ext<btw_lp> (
        band, _target_coeffs[band], freq, _t_spl, q_to_butterworth_order (b.q));
      break;
    case bandtype::butterworth_hp:
      _eq.reset_coeffs_ext<btw_hp> (
        band, _target_coeffs[band], freq, _t_spl, q_to_butterworth_order (b.q));
      break;
    case bandtype::svf_tilt:
      _eq.reset_coeffs_ext<tilt_eq> (
        band, _target_coeffs[band], freq, q, gain, _t_spl);
      break;
    case bandtype::presence: {
      q[1]              = q[0]; // undo diff
      auto qnorm_0_to_1 = q / get_parameter (band1_q_tag {}).max;
      _eq.reset_coeffs_ext<liteon::presence_high_shelf> (
        band, _target_coeffs[band], freq, qnorm_0_to_1, gain, _t_spl);
    } break;
    case bandtype::onepole_allpass: {
      _eq.reset_coeffs_ext<onepole_allpass> (
        band, _target_coeffs[band], freq, _t_spl);
    } break;
    case bandtype::svf_bell_bandpass: {
      _eq.reset_coeffs_ext<andy::svf> (
        band,
        _target_coeffs[band],
        freq,
        q,
        gain,
        _t_spl,
        bell_bandpass_tag {});
      break;
    } break;
    default:
      jassert (false);
    }
    // reset smoothing
    if (_cfg[band].reset_band_state) {
      _eq.reset_states_on_idx (band);
      xspan_copy<f64_x2> (_eq.get_coeffs (band), _target_coeffs[band]);
    }
    _cfg[band].has_changes      = false;
    _cfg[band].reset_band_state = false;
  }
  //----------------------------------------------------------------------------
  uint q_to_butterworth_order (double q)
  {
    auto order = (uint) q;
    order      = std::max (1u, order);
    order      = std::min (max_butterworth_order, order);
    return order;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_bands    = 4;
  static constexpr uint n_channels = 2;
  //----------------------------------------------------------------------------
  static constexpr uint max_butterworth_order = 6;
  //----------------------------------------------------------------------------
  enum class bandtype {
    off,
    svf_bell,
    svf_lshelf,
    svf_hshelf,
    svf_allpass,
    butterworth_lp,
    butterworth_hp,
    svf_tilt,
    presence,
    onepole_allpass,
    svf_bell_bandpass,
    size
  };
  //----------------------------------------------------------------------------
  struct bandconfig {
    bandtype type             = bandtype::off;
    float    freq_note        = constexpr_midi_note_to_hz (440.f);
    float    q                = 0.7f;
    float    gain_db          = 0.f;
    float    diff             = 0.f;
    bool     has_changes      = false;
    bool     reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  using btw_lp = butterworth_any_order<lowpass_tag, max_butterworth_order>;
  using btw_hp = butterworth_any_order<highpass_tag, max_butterworth_order>;

  using eqs = parts_union_array<
    mp_list<
      btw_lp,
      btw_hp,
      andy::svf,
      onepole_allpass,
      tilt_eq,
      liteon::presence_high_shelf>,
    f64_x2,
    n_bands>;
  //----------------------------------------------------------------------------
  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  enum class topology { stereo, l, r, lr_half, count };
  //----------------------------------------------------------------------------
  using coeff_array = std::array<eqs::value_type, eqs::n_coeffs>;
  std::array<bandconfig, n_bands> _cfg;

  alignas (eqs::value_type) std::array<coeff_array, n_bands> _target_coeffs;
  eqs      _eq;
  double   _smooth_coeff;
  topology _topology;

  float _t_spl;
};
//------------------------------------------------------------------------------
} // namespace artv
