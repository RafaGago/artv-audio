#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>

#include "artv-common/dsp/own/blocks/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/blocks/filters/onepole.hpp"
#include "artv-common/dsp/own/blocks/misc/interpolators.hpp"
#include "artv-common/dsp/own/blocks/waveshapers/hardclip.hpp"
#include "artv-common/dsp/own/blocks/waveshapers/sqrt_sigmoid.hpp"
#include "artv-common/dsp/own/blocks/waveshapers/sqrt_sin.hpp"
#include "artv-common/dsp/own/blocks/waveshapers/tanh.hpp"
#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

// TODO: -tanh
//       -parameter smoothing.
//       -As there is a downsampling FIR already, the high frequency correction
//        could be done simpultaneously on the downsampling filter by adding the
//        coefficients of a matched filter computed using Parks–McClellan/Remez
//        (scipy.signal.remez) with the same amount of taps. As both filters are
//        linear phase, the coefficients can just be summed (and the resulting
//        magnitude halved).

namespace artv {

//------------------------------------------------------------------------------
class saturation {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::waveshaper;
  //----------------------------------------------------------------------------
  struct type_tag {};

  void set (type_tag, int v) { _p.type = (decltype (_p.type)) v; }

  static constexpr auto get_parameter (type_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Sqrt",
        "Sqrt ADAA",
        "Tanh",
        "Tanh ADAA",
        "Hardclip",
        "Hardclip ADAA",
        "SqrtSin",
        "SqrtSin ADAA"),
      40);
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};

  void set (drive_tag, float v) { _p.drive = db_to_gain (v); }

  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -15.0, 15., 0.0, 0.25, 0.6, true);
  }
  //----------------------------------------------------------------------------
  struct compensated_drive_tag {};

  void set (compensated_drive_tag, float v)
  {
    _p.compensated_drive = db_to_gain (v);
  }

  static constexpr auto get_parameter (compensated_drive_tag)
  {
    return float_param ("dB", -20.0, 20, 0.0, 0.25, 0.6, true);
  }
  //----------------------------------------------------------------------------
  static constexpr float lo_cut_min_hz = 10.;

  struct lo_cut_tag {};

  void set (lo_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    if (v != _p.lo_cut_hz) {
      _p.lo_cut_hz = v;
      butterworth_type::lowpass (
        get_filt_coeffs (lo_lp), v, _plugcontext->get_sample_rate());
    }
  }

  static constexpr auto get_parameter (lo_cut_tag)
  {
    return frequency_parameter (lo_cut_min_hz, 10000., lo_cut_min_hz);
  }
  //----------------------------------------------------------------------------
  static constexpr float hi_cut_max_hz = 20000.;

  struct hi_cut_tag {};

  void set (hi_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    if (v != _p.hi_cut_hz) {
      _p.hi_cut_hz = v;
      butterworth_type::highpass (
        get_filt_coeffs (hi_hp), v, _plugcontext->get_sample_rate());
    }
  }

  static constexpr auto get_parameter (hi_cut_tag)
  {
    return frequency_parameter (60., hi_cut_max_hz, hi_cut_max_hz);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    drive_tag,
    compensated_drive_tag,
    type_tag,
    lo_cut_tag,
    hi_cut_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    pc.set_delay_compensation (1);
    memset (&_wvsh_states, 0, sizeof _wvsh_states);
    memset (&_filt_states, 0, sizeof _filt_states);
    memset (&_filt_coeffs, 0, sizeof _filt_coeffs);
    memset (&_allpass_states, 0, sizeof _allpass_states);

    allpass_interpolator::init<float> (make_crange (_allpass_coeff), 0.5);

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _p.type_prev = sat_type_count; // force initial click removal routine.
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint block_samples)
  {
    params p = _p;

    if (unlikely (p.type_prev != p.type)) {
      // some waveshapers will create peaks, as the integral on 0 might not be
      // 0. Running them for some samples of silence to initialize.
      static constexpr uint n_samples = 8;
      _p.type_prev                    = _p.type;
      std::array<T, 2> value_increment {
        chnls[0][0] * 1 / n_samples, chnls[1][0] * 1 / n_samples};
      std::array<T, n_samples * 2> in {};
      for (uint i = 1; i < n_samples; ++i) {
        in[i] = in[i - 1] + value_increment[0];
      }
      for (uint i = n_samples + 1; i < (n_samples * 2); ++i) {
        in[i] = in[i - 1] + value_increment[1];
      }
      process_block_replacing<T> ({&in[0], &in[n_samples]}, n_samples);
    }

    float inv_cdrive = 1 / p.compensated_drive;

    for (uint i = 0; i < block_samples; ++i) {
      // TODO: drive and filter change smoothing
      simd_dbl sat {chnls[0][i], chnls[1][i]};
      sat *= simd_dbl {p.drive * p.compensated_drive};

      simd_dbl lo {0.};
      simd_dbl hi {0.};

      if (p.lo_cut_hz > lo_cut_min_hz) {
        lo = butterworth_type::tick (
          get_filt_coeffs (lo_lp),
          {get_filt_states (lo_lp, 0), get_filt_states (lo_lp, 1)},
          {sat[0], sat[1]});
        sat -= lo;
        // 1 sample delay mix. Hack by retrieving a previous state...
        static_assert (std::is_same_v<butterworth_type, butterworth<1>>, "");
        lo = simd_dbl {
          get_filt_states (lo_lp, 0)[onepole::z1],
          get_filt_states (lo_lp, 1)[onepole::z1]};
      }

      if (p.hi_cut_hz < hi_cut_max_hz) {
        hi = butterworth_type::tick (
          get_filt_coeffs (hi_hp),
          {get_filt_states (hi_hp, 0), get_filt_states (hi_hp, 1)},
          {sat[0], sat[1]});
        sat -= hi;
        // 1 sample delay mix. Hack by retrieving a previous state...
        static_assert (std::is_same_v<butterworth_type, butterworth<1>>, "");
        hi = simd_dbl {
          get_filt_states (hi_hp, 0)[onepole::z1],
          get_filt_states (hi_hp, 1)[onepole::z1]};
      }

      // ARTV_SATURATION_USE_SSE has detrimental effects when ADAA enabled, the
      // current implementation always takes both branches.
      // The code is kept as a reminder for future me that this was already
      // tried.
      switch (p.type) {
      case sat_sqrt: {
        sat = sqrt_waveshaper_adaa<0>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
      } break;
      case sat_sqrt_adaa: {
#if ARTV_SATURATION_USE_SSE
        sat = sqrt_waveshaper_adaa<adaa_order>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
#else
        sat[0] = sqrt_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (0), sat[0]);
        sat[1] = sqrt_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (1), sat[1]);
#endif
        // Found empirically. TODO: improve
        static constexpr double gain = constexpr_db_to_gain (0.2);
        sat *= simd_dbl {gain};
      } break;
      case sat_tanh: {
        // At the point of writing my code or xsimd seems broken on tanh? works
        // well with others.
        sat[0] = tanh_waveshaper_adaa<0>::tick (
          {}, get_waveshaper_states (0), sat[0]);
        sat[1] = tanh_waveshaper_adaa<0>::tick (
          {}, get_waveshaper_states (1), sat[1]);
        // Found empirically. TODO: improve
        // static constexpr double gain = constexpr_db_to_gain (-2.6);
        // sat *= simd_dbl {gain};
      } break;
      case sat_tanh_adaa: {
#if ARTV_SATURATION_USE_SSE
        // At the point of writing my code or xsimd seems broken on tanh? works
        // well with others.
        sat = tanh_waveshaper_adaa<adaa_order>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
#else
        sat[0] = tanh_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (0), sat[0]);
        sat[1] = tanh_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (1), sat[1]);
#endif
        // Found empirically. TODO: improve
        // static constexpr double gain = constexpr_db_to_gain (-2.6);
        // sat *= simd_dbl {gain};
      } break;
      case sat_hardclip:
        sat = hardclip_waveshaper_adaa<0>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
        break;
      case sat_hardclip_adaa:
#if ARTV_SATURATION_USE_SSE
        sat = hardclip_waveshaper_adaa<adaa_order>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
#else
        sat[0] = hardclip_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (0), sat[0]);
        sat[1] = hardclip_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (1), sat[1]);
#endif
        break;
      case sat_sqrt_sin:
        sat = sqrt_sin_waveshaper_adaa<0>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
        break;
      case sat_sqrt_sin_adaa:
#if ARTV_SATURATION_USE_SSE
        sat = sqrt_sin_waveshaper_adaa<adaa_order>::tick_aligned<sse_bytes> (
          {},
          make_crange (_wvsh_states),
          make_crange ((const double*) &sat[0], decltype (sat)::size));
#else
        sat[0] = sqrt_sin_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (0), sat[0]);
        sat[1] = sqrt_sin_waveshaper_adaa<adaa_order>::tick (
          {}, get_waveshaper_states (1), sat[1]);
#endif
        break;
      default:
        break;
      }

      if constexpr (adaa_order == 1) {
        // half sample delay
        sat[0] = allpass_interpolator::tick<float> (
          make_crange (_allpass_coeff),
          make_crange (_allpass_states[0]),
          sat[0]);
        sat[1] = allpass_interpolator::tick<float> (
          make_crange (_allpass_coeff),
          make_crange (_allpass_states[1]),
          sat[1]);
      }
      sat += lo + hi;
      sat *= inv_cdrive;
      chnls[0][i] = sat[0];
      chnls[1][i] = sat[1];
    };
  }
  //----------------------------------------------------------------------------
private:
  enum sat_type {
    sat_sqrt,
    sat_sqrt_adaa,
    sat_tanh,
    sat_tanh_adaa,
    sat_hardclip,
    sat_hardclip_adaa,
    sat_sqrt_sin,
    sat_sqrt_sin_adaa,
    sat_type_count,
  };

  enum filter_indexes { lo_lp, hi_hp, n_filters };

  struct params {
    float drive             = 1.f;
    float compensated_drive = 1.f;
    float lo_cut_hz         = -1.f;
    float hi_cut_hz         = -1.f;
    char  type              = 0;
    char  type_prev         = 1;
  } _p;

  static constexpr uint adaa_order = 1;

  using shapers = mp_list<
    sqrt_waveshaper_adaa<adaa_order>,
    tanh_waveshaper_adaa<adaa_order>,
    hardclip_waveshaper_adaa<adaa_order>,
    sqrt_sin_waveshaper_adaa<adaa_order>>;

  template <class T>
  using wsh_to_n_states    = std::integral_constant<int, T::n_states>;
  using wsh_n_states_types = mp11::mp_transform<wsh_to_n_states, shapers>;

  static constexpr uint wsh_max_states
    = mp11::mp_max_element<wsh_n_states_types, mp11::mp_less>::value;

  using smoother = onepole_smoother;
  //----------------------------------------------------------------------------
  using butterworth_type = butterworth<1>;

  crange<double> get_filt_coeffs (uint filt_idx)
  {
    static constexpr uint n_coeffs = butterworth_type::n_coeffs;
    return {&_filt_coeffs[filt_idx * n_coeffs], n_coeffs};
  }
  crange<double> get_filt_states (uint filt_idx, uint channel)
  {
    static constexpr uint n_states = butterworth_type::n_states;
    return {&_filt_states[channel][filt_idx * n_states], n_states};
  }

  crange<double> get_waveshaper_states (uint channel)
  {
    return {&_wvsh_states[wsh_max_states * channel], wsh_max_states};
  }

  //----------------------------------------------------------------------------
  static constexpr uint n_channels = 2;
  using wsh_state_array
    = simd_array<double, wsh_max_states * n_channels, sse_bytes>;
  // using coeff_array = simd_array<double, max_coeffs, sse_bytes>;
  using crossv_state_array
    = simd_array<double, butterworth_type::n_states * n_filters, sse_bytes>;
  using crossv_coeff_array
    = simd_array<double, butterworth_type::n_coeffs * n_filters, sse_bytes>;
  using allpass_state_array
    = simd_array<float, allpass_interpolator::n_states, sse_bytes>;

  alignas (sse_bytes) wsh_state_array _wvsh_states;
  crossv_coeff_array                 _filt_coeffs;
  std::array<crossv_state_array, 2>  _filt_states;
  std::array<allpass_state_array, 2> _allpass_states;
  static_assert (allpass_interpolator::n_coeffs == 1, "");
  float _allpass_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
