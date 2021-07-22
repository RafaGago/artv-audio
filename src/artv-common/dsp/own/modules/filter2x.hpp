#pragma once

#include <algorithm>

#include "artv-common/dsp/own/blocks/filters/andy_svf.hpp"
#include "artv-common/dsp/own/blocks/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/blocks/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/blocks/filters/onepole.hpp"
#include "artv-common/dsp/own/blocks/filters/presence.hpp"
#include "artv-common/dsp/own/blocks/filters/saike.hpp"
#include "artv-common/dsp/own/jsfx.hpp"
#include "artv-common/dsp/own/misc.hpp"
#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

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
  };

public:
  template <uint Band, paramtype Type>
  struct param {
    static constexpr uint band  = Band;
    static constexpr uint ptype = (uint) Type;
  };
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::eq;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _cfg         = decltype (_cfg) {};
    smoother::lowpass<double> (
      make_crange (_smooth_coeff), 1. / 0.08, pc.get_sample_rate());
    memset (&_smooth_state, 0, sizeof _smooth_state);
    memset (&_state, 0, sizeof _state);
    memset (&_last_sample, 0, sizeof _last_sample);
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    constexpr uint                          blocksize = 128;
    std::array<std::array<T, blocksize>, 2> in_cp;

    for (uint block_smp = 0; block_smp < samples; block_smp += blocksize) {
      uint n_samples = std::min<uint> (in_cp[0].size(), samples - block_smp);

      // both filters are parallel/operate on the same input
      memcpy (in_cp[0].data(), &chnls[0][block_smp], n_samples * sizeof (T));
      memcpy (in_cp[1].data(), &chnls[1][block_smp], n_samples * sizeof (T));
      memset (&chnls[0][block_smp], 0, n_samples * sizeof (T));
      memset (&chnls[1][block_smp], 0, n_samples * sizeof (T));

      for (uint b = 0; b < n_bands; ++b) {
        if (_cfg[b].has_changes) {
          reset_band (b);
        }

        auto btype = get_band_type (b);
        if (btype == bandtype::off) {
          continue;
        }

        alignas (sse_bytes) std::array<coeff_array, 2> smoothed_band_coefs;
        alignas (sse_bytes) smoothed_params            smoothed_p;

        for (uint i = 0; i < n_samples; ++i) {

          static_assert (
            sizeof _smooth_state[b].params / sizeof (float) == 4, "");

          smoothed_p.vars.pre_drive_l = _cfg[b].pre_drive_l;
          smoothed_p.vars.pre_drive_r = _cfg[b].pre_drive_r;
          smoothed_p.vars.post_drive  = _cfg[b].post_drive;
          smoothed_p.vars.feedback    = _cfg[b].feedback;

          auto     smcoeff = (float) _smooth_coeff;
          simd_flt params  = smoother::tick_aligned<sse_bytes, float> (
            make_crange (smcoeff),
            _smooth_state[b].params.arr,
            simd_flt {smoothed_p.arr.data(), xsimd::aligned_mode {}});

          memcpy (smoothed_p.arr.data(), &params, sizeof params);

          constexpr uint sse_step = simd_dbl::size;
          static_assert (smoothed_band_coefs.size() % sse_step == 0, "");
          for (uint chnl = 0; chnl < 2; ++chnl) {
            for (uint j = 0; j < smoothed_band_coefs[0].size(); j += sse_step) {
              simd_dbl out = smoother::tick_aligned<sse_bytes, double> (
                make_crange (_smooth_coeff),
                make_crange (&_smooth_state[b].filter[chnl][j], sse_step),
                simd_dbl {&_coeffs[b][chnl][j], xsimd::aligned_mode {}});
              out.store_aligned (&smoothed_band_coefs[chnl][j]);
            }
          }

          simd_dbl in, out;
          in[0] = in_cp[0][i];
          in[1] = in_cp[1][i];
          in *= simd_dbl {
            smoothed_p.vars.pre_drive_l, smoothed_p.vars.pre_drive_r};

          simd_dbl prev {_last_sample[b][0], _last_sample[b][1]};
          prev *= simd_dbl {smoothed_p.vars.feedback};
          bool is_moog
            = btype >= bandtype::moog_1_lp && btype <= bandtype::moog_2_br;
          // moogs are extremely unstable reduce feedback
          prev *= simd_dbl {is_moog ? 0.03 : 1.};
          in += prev;

          switch (btype) {
          case bandtype::off: // not reachable
            jassert (false);
            break;
          case bandtype::ms20_lp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_lowpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_hp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_highpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_bp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_bandpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_br:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_notch::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_asym_lp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_asym_lowpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_asym_hp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_asym_highpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_asym_bp:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_asym_bandpass::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::ms20_asym_br:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::ms20_asym_notch::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::steiner_1_lp:
          case bandtype::steiner_1_hp:
          case bandtype::steiner_1_bp:
          case bandtype::steiner_1_br:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::steiner_1::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::steiner_2_lp:
          case bandtype::steiner_2_hp:
          case bandtype::steiner_2_bp:
          case bandtype::steiner_2_br:
            for (uint ch = 0; ch < 2; ++ch) {
              out[ch] = saike::steiner_2::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::moog_1_lp:
          case bandtype::moog_1_hp:
          case bandtype::moog_1_bp:
          case bandtype::moog_1_br:
            for (uint ch = 0; ch < 2; ++ch) {
              saike::moog_1::repair_unsmoothable_coeffs (
                smoothed_band_coefs[ch], _coeffs[b][ch]);
              out[ch] = saike::moog_1::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          case bandtype::moog_2_lp:
          case bandtype::moog_2_hp:
          case bandtype::moog_2_bp:
          case bandtype::moog_2_br:
            for (uint ch = 0; ch < 2; ++ch) {
              saike::moog_2::repair_unsmoothable_coeffs (
                smoothed_band_coefs[ch], _coeffs[b][ch]);
              out[ch] = saike::moog_2::tick (
                smoothed_band_coefs[ch], _state[b][ch], in[ch]);
            }
            break;
          default:
            jassert (false);
            break;
          }

          _last_sample[b][0] = out[0];
          _last_sample[b][1] = out[1];

          out *= smoothed_p.vars.post_drive;
          chnls[0][block_smp + i] += out[0];
          chnls[1][block_smp + i] += out[1];
        }
      }
    }
  }
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
      0,
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
    v *= 0.0093;
    _cfg[band].feedback = v;
  }

  template <uint band>
  static constexpr auto get_parameter (param<band, paramtype::feedback>)
  {
    return float_param ("%", -100.f, 100.f, 0.f, 0.1f);
  }
  //----------------------------------------------------------------------------
  using band1_type_tag      = param<0, paramtype::band_type>;
  using band2_type_tag      = param<1, paramtype::band_type>;
  using band1_topology_tag  = param<0, paramtype::topology>;
  using band2_topology_tag  = param<1, paramtype::topology>;
  using band1_freq_tag      = param<0, paramtype::frequency>;
  using band2_freq_tag      = param<1, paramtype::frequency>;
  using band1_reso_tag      = param<0, paramtype::reso>;
  using band2_reso_tag      = param<1, paramtype::reso>;
  using band1_drive_tag     = param<0, paramtype::drive>;
  using band2_drive_tag     = param<1, paramtype::drive>;
  using band1_gain_tag      = param<0, paramtype::gain>;
  using band2_gain_tag      = param<1, paramtype::gain>;
  using band1_tolerance_tag = param<0, paramtype::tolerance>;
  using band2_tolerance_tag = param<1, paramtype::tolerance>;
  using band1_feedback_tag  = param<0, paramtype::feedback>;
  using band2_feedback_tag  = param<1, paramtype::feedback>;
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
    band2_feedback_tag>;
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
  struct bandconfig {
    float freq             = 440.f;
    float reso             = 0.7f;
    float drive_db         = 0.f;
    float gain_db          = 0.f;
    float tolerance        = 0.f;
    float feedback         = 0.f;
    float dry_wet          = 1.f;
    float pre_drive_l      = 1.f;
    float pre_drive_r      = 1.f;
    float post_drive       = 1.f;
    u8    type             = 0;
    u8    topology         = 0;
    bool  has_changes      = false;
    bool  reset_band_state = false;
  };
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
    auto& b  = _cfg[band];
    auto  sr = (float) _plugcontext->get_sample_rate();

    std::array<float, 2> freq, reso;
    freq[0] = b.freq;
    freq[1] = b.freq * (1.f + (b.tolerance * 0.2f));

    freq[1] = std::clamp<float> (
      freq[1],
      midi_note_to_hz (get_parameter (band1_freq_tag {}).min),
      midi_note_to_hz (get_parameter (band1_freq_tag {}).max));

    reso[0] = b.reso;
    reso[1] = b.reso * (1.f - (b.tolerance * 0.11f));

    reso[1] = std::clamp (
      reso[1],
      get_parameter (band1_reso_tag {}).min,
      get_parameter (band1_reso_tag {}).max);

    auto btype = get_band_type (band);
    switch (btype) {
    case bandtype::off:
      memset (&_coeffs[band], 0, sizeof _coeffs[band]);
      break;
    case bandtype::ms20_lp:
      saike::ms20_lowpass::lowpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_lowpass::lowpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_hp:
      saike::ms20_highpass::highpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_highpass::highpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_bp:
      saike::ms20_bandpass::bandpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_bandpass::bandpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_br:
      saike::ms20_notch::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_notch::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_asym_lp:
      saike::ms20_asym_lowpass::lowpass (
        _coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_asym_lowpass::lowpass (
        _coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_asym_hp:
      saike::ms20_asym_highpass::highpass (
        _coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_asym_highpass::highpass (
        _coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_asym_bp:
      saike::ms20_asym_bandpass::bandpass (
        _coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_asym_bandpass::bandpass (
        _coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::ms20_asym_br:
      saike::ms20_asym_notch::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::ms20_asym_notch::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_1_lp:
      saike::steiner_1::lowpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_1::lowpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_1_hp:
      saike::steiner_1::highpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_1::highpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_1_bp:
      saike::steiner_1::bandpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_1::bandpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_1_br:
      saike::steiner_1::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_1::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_2_lp:
      saike::steiner_2::lowpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_2::lowpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_2_hp:
      saike::steiner_2::highpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_2::highpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_2_bp:
      saike::steiner_2::bandpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_2::bandpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::steiner_2_br:
      saike::steiner_2::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::steiner_2::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_1_lp:
      saike::moog_1::lowpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_1::lowpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_1_hp:
      saike::moog_1::highpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_1::highpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_1_bp:
      saike::moog_1::bandpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_1::bandpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_1_br:
      saike::moog_1::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_1::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_2_lp:
      saike::moog_2::lowpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_2::lowpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_2_hp:
      saike::moog_2::highpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_2::highpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_2_bp:
      saike::moog_2::bandpass (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_2::bandpass (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    case bandtype::moog_2_br:
      saike::moog_2::notch (_coeffs[band][0], freq[0], reso[0], sr);
      saike::moog_2::notch (_coeffs[band][1], freq[1], reso[1], sr);
      break;
    default:
      jassert (false);
    }
    // drive handling
    // The idea is to hit the filter channel at 0dBVU and that the drive has
    // an useful range at that level. I adjust with a Q before self-oscillation
    // and try to slightly correct very wild gain differences. I used a loop
    // hoovering on trakmeter's 0dB RMS and Reaper's white noise get at -18dB
    if (
      is_ms20 (btype) || is_ms20_asym (btype) || is_steiner_1 (btype)
      || is_steiner_2 (btype)) {
      float drive            = db_to_gain (b.drive_db);
      _cfg[band].pre_drive_l = (drive * 16.);
      // wanting Drive to act as a weak gain on the positive side too so
      // applying on the numerator:
      // (((x^2) / 64 + (63/64), which crosses at (x=1, y=1).
      //
      // Instead of just applying (1. / pre_drive).
      _cfg[band].post_drive
        = (((drive * drive) * 0.015625) + 0.984375) / _cfg[band].pre_drive_l;
    }
    else if (is_moog_1 (btype) || is_moog_2 (btype)) {
      float drive            = db_to_gain (b.drive_db);
      _cfg[band].pre_drive_l = drive;
      // wanting Drive to act as a weak gain on the positive side too so
      // applying on the numerator:
      // (((x^2) / 64 + (63/64), which crosses at (x=1, y=1).
      //
      // Instead of just applying (1. / drive).
      _cfg[band].post_drive
        = (((drive * drive) * 0.015625) + 0.984375) / _cfg[band].pre_drive_l;
      // These have the phase polarity inverted, so using them in parallel
      // caused cancellations
      _cfg[band].post_drive = -_cfg[band].post_drive;
      // these come very hot too, so they need taming.
      double tame = 1.;
      switch (btype) {
        // TODO: substitute these corections by constants when they have
        // settled.
      case bandtype::moog_1_lp: {
        constexpr double tame_k = constexpr_db_to_gain (-23.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_1_hp: {
        _cfg[band].pre_drive_l *= 0.25;
        constexpr double tame_k = constexpr_db_to_gain (-37.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_1_bp: {
        _cfg[band].pre_drive_l *= 0.25;
        constexpr double tame_k = constexpr_db_to_gain (-44.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_1_br: {
        _cfg[band].pre_drive_l *= 0.25;
        constexpr double tame_k = constexpr_db_to_gain (-39.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_2_lp: {
        constexpr double tame_k = constexpr_db_to_gain (-25.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_2_hp: {
        constexpr double tame_k = constexpr_db_to_gain (-43.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_2_bp: {
        constexpr double tame_k = constexpr_db_to_gain (-62.);
        tame                    = tame_k;
      } break;
      case bandtype::moog_2_br: {
        constexpr double tame_k = constexpr_db_to_gain (-44.);
        tame                    = tame_k;
      } break;
      }
      _cfg[band].post_drive *= tame;
    }
    else {
      _cfg[band].pre_drive_l = _cfg[band].post_drive = 1.f;
    }
    _cfg[band].post_drive *= db_to_gain (_cfg[band].gain_db);
    _cfg[band].pre_drive_r
      = _cfg[band].pre_drive_l * (1.f - (b.tolerance * 0.13f));

    // reset smoothing
    if (_cfg[band].reset_band_state) {
      memset (&_state[band][0], 0, sizeof _state[band][0]);
      memset (&_state[band][1], 0, sizeof _state[band][1]);
      memset (&_smooth_state[band], 0, sizeof _smooth_state[band]);
    }
    _cfg[band].has_changes      = false;
    _cfg[band].reset_band_state = false;
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_bands = 4;
  //----------------------------------------------------------------------------
  using filters = mp_list<
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
    saike::moog_2>;

  template <class T>
  using to_n_coeffs = std::integral_constant<int, T::n_coeffs>;
  template <class T>
  using to_n_states = std::integral_constant<int, T::n_states>;

  using n_coeffs_types = mp11::mp_transform<to_n_coeffs, filters>;
  using n_states_types = mp11::mp_transform<to_n_states, filters>;

  static constexpr uint max_coeffs
    = mp11::mp_max_element<n_coeffs_types, mp11::mp_less>::value;

  static constexpr uint max_states
    = mp11::mp_max_element<n_states_types, mp11::mp_less>::value;

  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  using state_array = simd_array<double, max_states, sse_bytes>;
  using coeff_array = simd_array<double, max_coeffs, sse_bytes>;
  //----------------------------------------------------------------------------
  union smoothed_params {
    struct {
      float pre_drive_l, pre_drive_r, post_drive, feedback;
    } vars;
    std::array<float, 4> arr;
  };
  //----------------------------------------------------------------------------
  struct smooth_states {
    alignas (sse_bytes) std::array<coeff_array, 2> filter;
    smoothed_params params;
  };
  //----------------------------------------------------------------------------
  std::array<bandconfig, n_bands>                 _cfg;
  std::array<std::array<state_array, 2>, n_bands> _state;
  std::array<smooth_states, n_bands>              _smooth_state;
  alignas (sse_bytes) std::array<std::array<coeff_array, 2>, n_bands> _coeffs;
  std::array<std::array<double, 2>, n_bands> _last_sample;

  double _smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
