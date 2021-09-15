#pragma once

#include "artv-common/dsp/own/classes/jsfx.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/filters/presence.hpp"
#include "artv-common/dsp/own/parts/filters/saike.hpp"
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
  enum class paramtype { band_type, frequency, q, gain };

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
    _plugcontext = &pc;
    _cfg         = decltype (_cfg) {};
    smoother::reset_coeffs (
      make_crange (_smooth_coeff),
      make_vec_x1 (1. / 0.02),
      pc.get_sample_rate());
    memset (&_smooth_state, 0, sizeof _smooth_state);
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
      alignas (sse_bytes) coeff_array smoothed_band_coefs;

      for (uint i = 0; i < samples; ++i) {
        constexpr uint sse_step = vec_traits<double_x2>().size;
        static_assert (smoothed_band_coefs.size() % sse_step == 0, "");
        for (uint j = 0; j < smoothed_band_coefs.size(); j += sse_step) {
          // coefficient change smoothing.TODO: if filters start differing
          // broadly on coefficient sizes this might have to be optimized
          // (maybe to a) small/medium/high iteration or to the exact number
          // of elements.
          double_x2 out = smoother::tick (
            make_crange (_smooth_coeff),
            make_crange (&_smooth_state[b].filter[j], sse_step),
            vec_load<double_x2> (&_coeffs[b][j]),
            single_coeff_set_tag {});
          vec_store (&smoothed_band_coefs[j], out);
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
          out = andy::svf::tick (
            smoothed_band_coefs, _state[b], in, single_coeff_set_tag {});
          break;
        case bandtype::butterworth_lp:
          out = butterworth_lowpass_any_order::tick (
            smoothed_band_coefs,
            _state[b],
            in,
            q_to_butterworth_order (_cfg[b].q),
            single_coeff_set_tag {});
          break;
        case bandtype::butterworth_hp:
          out = butterworth_highpass_any_order::tick (
            smoothed_band_coefs,
            _state[b],
            in,
            q_to_butterworth_order (_cfg[b].q),
            single_coeff_set_tag {});
          break;
        case bandtype::svf_tilt:
          out = tilt_eq::tick (
            smoothed_band_coefs, _state[b], in, single_coeff_set_tag {});
          break;
        case bandtype::presence:
          out = liteon::presence_high_shelf::tick (
            smoothed_band_coefs, _state[b], in, single_coeff_set_tag {});
          break;
        case bandtype::onepole_allpass:
          out = onepole::tick (
            smoothed_band_coefs, _state[b], in, single_coeff_set_tag {});
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
    band4_q_tag>;
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
    bool     has_changes      = false;
    bool     reset_band_state = false;
  };
  //----------------------------------------------------------------------------
  void reset_band (uint band)
  {
    auto& b             = _cfg[band];
    auto  bandtype_prev = b.type;
    auto  sr            = (float) _plugcontext->get_sample_rate();
    switch (b.type) {
    case bandtype::off:
      memset (&_coeffs[band], 0, sizeof _coeffs[band]);
      break;
    case bandtype::svf_bell:
      andy::svf::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (b.q),
        make_vec_x1<double> (b.gain_db),
        sr,
        bell_tag {});
      break;
    case bandtype::svf_lshelf:
      andy::svf::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (b.q),
        make_vec_x1<double> (b.gain_db),
        sr,
        lowshelf_tag {});
      break;
    case bandtype::svf_hshelf:
      andy::svf::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (b.q),
        make_vec_x1<double> (b.gain_db),
        sr,
        highshelf_tag {});
      break;
    case bandtype::svf_allpass:
      andy::svf::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (b.q),
        sr,
        allpass_tag {});
      break;
    case bandtype::butterworth_lp:
      butterworth_lowpass_any_order::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        sr,
        q_to_butterworth_order (b.q));
      break;
    case bandtype::butterworth_hp:
      butterworth_highpass_any_order::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        sr,
        q_to_butterworth_order (b.q));
      break;
    case bandtype::svf_tilt:
      tilt_eq::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (b.q),
        make_vec_x1<double> (b.gain_db),
        sr);
      break;
    case bandtype::presence: {
      float qnorm_0_to_1 = b.q / get_parameter (band1_q_tag {}).max;
      liteon::presence_high_shelf::reset_coeffs (
        _coeffs[band],
        make_vec_x1<double> (b.freq),
        make_vec_x1<double> (qnorm_0_to_1),
        make_vec_x1<double> (b.gain_db),
        sr);
    } break;
    case bandtype::onepole_allpass: {
      onepole::reset_coeffs (
        _coeffs[band], make_vec_x1<double> (b.freq), sr, allpass_tag {});
    } break;
    default:
      jassert (false);
    }
    // reset smoothing
    if (_cfg[band].reset_band_state) {
      memset (&_state[0][band], 0, sizeof _state[0][band]);
      memset (&_state[1][band], 0, sizeof _state[1][band]);
      memset (&_smooth_state[band], 0, sizeof _smooth_state[band]);
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
  using filters = mp_list<
    andy::svf,
    onepole,
    butterworth_lowpass<max_butterworth_order>,
    butterworth_highpass<max_butterworth_order>,
    liteon::presence_high_shelf>;

  template <class T>
  using to_n_coeffs = k_int<T::n_coeffs>;
  template <class T>
  using to_n_states = k_int<T::n_states>;

  using n_coeffs_types = mp11::mp_transform<to_n_coeffs, filters>;
  using n_states_types = mp11::mp_transform<to_n_states, filters>;

  static constexpr uint max_coeffs
    = mp11::mp_max_element<n_coeffs_types, mp11::mp_less>::value;

  static constexpr uint max_states
    = mp11::mp_max_element<n_states_types, mp11::mp_less>::value;

  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  using state_array = simd_array<double, max_states * n_channels, sse_bytes>;
  using coeff_array = simd_array<double, max_coeffs, sse_bytes>;
  //----------------------------------------------------------------------------
  struct smooth_states {
    alignas (sse_bytes) coeff_array filter;
  };
  //----------------------------------------------------------------------------
  std::array<bandconfig, n_bands>    _cfg;
  std::array<smooth_states, n_bands> _smooth_state;
  alignas (sse_bytes) std::array<coeff_array, n_bands> _coeffs;
  alignas (sse_bytes) std::array<state_array, n_bands> _state;
  double _smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
