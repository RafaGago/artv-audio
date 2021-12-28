#pragma once

#include "artv-common/dsp/own/classes/jsfx.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
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
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

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

    _plugcontext = &pc;
    _cfg         = decltype (_cfg) {};
    smoother::reset_coeffs (
      make_crange (_smooth_coeff).cast (x1_t {}),
      vec_set<x1_t> (1. / 0.02),
      pc.get_sample_rate());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    uint n_enabled_bands = 0;

    for (uint b = 0; b < n_bands; ++b) {
      if (_cfg[b].has_changes) {
        reset_band (b);
      }
      if (_cfg[b].type == bandtype::off) {
        continue;
      }

      ++n_enabled_bands;

      for (uint i = 0; i < samples; ++i) {
        // coefficient change smoothing.TODO: if filters start differing
        // broadly on coefficient sizes this might have to be optimized
        // (maybe to a) small/medium/high iteration or to the exact number
        // of elements.
        crange<double_x2> internal = _eq.get_coeffs (b);
        for (uint j = 0; j < _target_coeffs[b].size(); ++j) {
          internal[j] = smoother::tick (
            make_crange (_smooth_coeff),
            make_crange (internal[j]),
            _target_coeffs[b][j]);
        }
        double_x2 in, out;
        in[0] = ins[0][i];
        in[1] = ins[1][i];

        switch (_cfg[b].type) {
        case bandtype::off: // not reachable
        case bandtype::svf_bell:
        case bandtype::svf_lshelf:
        case bandtype::svf_hshelf:
        case bandtype::svf_allpass:
          out = _eq.tick_on_idx<andy::svf> (b, in);
          break;
        case bandtype::butterworth_lp:
          out = _eq.tick_on_idx<btw_lp> (
            b, in, q_to_butterworth_order (_cfg[b].q));
          break;
        case bandtype::butterworth_hp:
          out = _eq.tick_on_idx<btw_hp> (
            b, in, q_to_butterworth_order (_cfg[b].q));
          break;
        case bandtype::svf_tilt:
          out = _eq.tick_on_idx<tilt_eq> (b, in);
          break;
        case bandtype::presence:
          out = _eq.tick_on_idx<liteon::presence_high_shelf> (b, in);
          break;
        case bandtype::onepole_allpass:
          out = _eq.tick_on_idx<onepole_allpass> (b, in);
          break;
        default:
          jassert (false);
          break;
        }
        outs[0][i] = out[0];
        outs[1][i] = out[1];
      }
    }
    if (unlikely (n_enabled_bands == 0)) {
      for (uint i = 0; i < (uint) bus_type; ++i) {
        memcpy (outs[i], ins[i], samples * sizeof outs[0][0]);
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
    _cfg[band].has_changes |= (int) _cfg[band].type != (int) v;
    _cfg[band].reset_band_state |= _cfg[band].has_changes;
    _cfg[band].type = (bandtype) v;
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
        "Allpass OnePole"),
      30);
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
  void set (param<band, paramtype::q>, float v)
  {
    static_assert (band < n_bands, "");
    auto& b = _cfg[band];
    // changing butterworth slope requires a state reset, as they are different
    // filters under the same menu.
    bool reset = (b.type == bandtype::butterworth_hp);
    reset |= (b.type == bandtype::butterworth_lp);
    reset &= (q_to_butterworth_order (b.q) != q_to_butterworth_order (v));

    b.reset_band_state |= reset;
    b.has_changes |= reset;
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
    band4_diff_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  enum class bandtype {
    off,
    svf_bell,
    svf_lshelf,
    svf_hshelf,
    svf_nodrive_last = svf_hshelf,
    svf_allpass,
    butterworth_lp,
    butterworth_hp,
    svf_tilt,
    presence,
    onepole_allpass,
    size
  };
  //----------------------------------------------------------------------------
  struct bandconfig {
    bandtype type             = bandtype::off;
    float    freq             = 440.f;
    float    q                = 0.7f;
    float    gain_db          = 0.f;
    float    diff             = 0.f;
    bool     has_changes      = false;
    bool     reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    using x1_t = vec<double, 1>;

    auto& b             = _cfg[band];
    auto  bandtype_prev = b.type;
    auto  sr            = (float) _plugcontext->get_sample_rate();

    auto freq = vec_set<double_x2> (b.freq);
    freq[1] *= exp2 (b.diff);
    auto q = vec_set<double_x2> (b.q);
    q[1] += q[1] * b.diff * 0.05;
    auto gain = vec_set<double_x2> (b.gain_db);

    switch (b.type) {
    case bandtype::off:
      _eq.reset_states_on_idx (band);
      break;
    case bandtype::svf_bell:
      _eq.reset_target_coeffs<andy::svf> (
        _target_coeffs[band], freq, q, gain, sr, bell_tag {});
      break;
    case bandtype::svf_lshelf:
      _eq.reset_target_coeffs<andy::svf> (
        _target_coeffs[band], freq, q, gain, sr, lowshelf_tag {});
      break;
    case bandtype::svf_hshelf:
      _eq.reset_target_coeffs<andy::svf> (
        _target_coeffs[band], freq, q, gain, sr, highshelf_tag {});
      break;
    case bandtype::svf_allpass:
      _eq.reset_target_coeffs<andy::svf> (
        _target_coeffs[band], freq, q, sr, allpass_tag {});
      break;
    case bandtype::butterworth_lp:
      _eq.reset_target_coeffs<btw_lp> (
        _target_coeffs[band], freq, sr, q_to_butterworth_order (b.q));
      break;
    case bandtype::butterworth_hp:
      _eq.reset_target_coeffs<btw_hp> (
        _target_coeffs[band], freq, sr, q_to_butterworth_order (b.q));
      break;
    case bandtype::svf_tilt:
      _eq.reset_target_coeffs<tilt_eq> (
        _target_coeffs[band], freq, q, gain, sr);
      break;
    case bandtype::presence: {
      q[1]              = q[0]; // undo diff
      auto qnorm_0_to_1 = q / get_parameter (band1_q_tag {}).max;
      _eq.reset_target_coeffs<liteon::presence_high_shelf> (
        _target_coeffs[band], freq, qnorm_0_to_1, gain, sr);
    } break;
    case bandtype::onepole_allpass: {
      _eq.reset_target_coeffs<onepole_allpass> (_target_coeffs[band], freq, sr);
    } break;
    default:
      jassert (false);
    }
    // reset smoothing
    if (_cfg[band].reset_band_state) {
      _eq.reset_states_on_idx (band);
      crange_copy<double_x2> (_eq.get_coeffs (band), _target_coeffs[band]);
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
  using btw_lp = butterworth_any_order<lowpass_tag, max_butterworth_order>;
  using btw_hp = butterworth_any_order<highpass_tag, max_butterworth_order>;

  using eqs = parts_to_class_one_of<
    double_x2,
    n_bands,
    btw_lp,
    btw_hp,
    andy::svf,
    onepole_allpass,
    tilt_eq,
    liteon::presence_high_shelf>;
  //----------------------------------------------------------------------------
  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  std::array<bandconfig, n_bands>       _cfg;
  std::array<eqs::coeff_array, n_bands> _target_coeffs;
  eqs                                   _eq;
  double                                _smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
