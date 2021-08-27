#pragma once

#include <array>
#include <cstdint>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::modulation;
  //----------------------------------------------------------------------------
  struct stages_tag {};
  void set (stages_tag, int v)
  {
    static_assert (vec_traits<float_x4>().size == 4, "Refactor this...");
    auto stages                    = (v + 1) * vec_traits<float_x4>().size;
    _params.unsmoothed.n_allpasses = stages;
  }

  static constexpr auto get_parameter (stages_tag)
  {
    // TODO: before release, not so many choices make a difference, reduce to
    // e.g. 2, 4, 6, 8, 12, 16, 24, 32
    return choice_param (
      6,
      make_cstr_array (
        "2",
        "4",
        "6",
        "8",
        "10",
        "12",
        "14",
        "16",
        "18",
        "20",
        "22",
        "24",
        "26",
        "28",
        "30",
        "32"),
      32);
  }
  //----------------------------------------------------------------------------
  struct stages_mode_tag {};
  void set (stages_mode_tag, int v) { _params.unsmoothed.mode = v; }

  static constexpr auto get_parameter (stages_mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Exponential",
        "Exp Spread",
        "Exp Alternate",
        "Exp Stereo",
        "Linear",
        "Lin Spread",
        "Lin Alternate",
        "Lin Stereo"),
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
    _params.smooth_target.value.lfo_depth = v;
  }

  static constexpr auto get_parameter (lfo_depth_tag)
  {
    return float_param ("%", 0., 100., 25., 0.1);
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
    _params.smooth_target.value.lfo_start_phase = v;
  }

  static constexpr auto get_parameter (lfo_start_phase_tag)
  {
    return float_param ("deg", 0., 360., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct lfo_stereo_tag {};
  void set (lfo_stereo_tag, float v)
  {
    _params.smooth_target.value.lfo_stereo = v;
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
    _params.smooth_target.value.freq_lo = v;
  }

  static constexpr auto get_parameter (low_freq_tag)
  {
    return frequency_parameter (20., 10000., 225.);
  }
  //----------------------------------------------------------------------------
  struct high_freq_tag {};
  void set (high_freq_tag, float v)
  {
    v                                   = midi_note_to_hz (v);
    _params.smooth_target.value.freq_hi = v;
  }

  static constexpr auto get_parameter (high_freq_tag)
  {
    return frequency_parameter (20., 10000., 3300.);
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.00995;
    v = sqrt (abs (v));
    v = neg ? -v : v;

    _params.smooth_target.value.feedback = v;
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct q_tag {};
  void set (q_tag, float v) { _params.smooth_target.value.q = v; }

  static constexpr auto get_parameter (q_tag)
  {
    return float_param ("", 0.02, 3.5, 0.18, 0.001, 0.4);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    memset (&_allpass, 0, sizeof _allpass);
    memset (&_feedback_samples, 0, sizeof _feedback_samples);

    _lfos[0].reset();
    _lfos[1].reset();

    _params.smooth_target.value.freq_lo
      = get_parameter (low_freq_tag {}).defaultv;
    _params.smooth_target.value.freq_hi
      = get_parameter (high_freq_tag {}).defaultv;
    _params.smooth_target.value.lfo_stereo
      = get_parameter (lfo_stereo_tag {}).defaultv;
    _params.smooth_target.value.lfo_depth
      = get_parameter (lfo_depth_tag {}).defaultv;
    _params.smooth_target.value.lfo_start_phase
      = get_parameter (lfo_start_phase_tag {}).defaultv;
    _params.smooth_target.value.feedback
      = get_parameter (feedback_tag {}).defaultv;
    _params.smooth_target.value.q = get_parameter (q_tag {}).defaultv;

    _params.unsmoothed.lfo_hz_user = get_parameter (lfo_rate_tag {}).defaultv;
    _params.unsmoothed.lfo_eights
      = get_parameter (lfo_rate_sync_tag {}).defaultv;
    _params.unsmoothed.lfo_wave    = get_parameter (lfo_wave_tag {}).defaultv;
    _params.unsmoothed.n_allpasses = get_parameter (stages_tag {}).defaultv;
    _params.unsmoothed.mode = get_parameter (stages_mode_tag {}).defaultv;

    refresh_lfo_hz();

    memset (&_params.smooth_state, 0, sizeof _params.smooth_state);
    memset (&_dc_block_states, 0, sizeof _dc_block_states);

    mystran_dc_blocker::reset_coeffs (
      _dc_block_coeffs, vec_set<double_x2> (2.), pc.get_sample_rate());

    onepole_smoother::reset_coeffs (
      make_crange (_lp_smooth_coeff),
      make_vec_x1<decltype (_lp_smooth_coeff)> (1. / 0.01),
      pc.get_sample_rate());

    _n_processed_samples = 0;

    // Sample rates 44100 multiples update every 362.811us
    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _feedback_samples[0] = _feedback_samples[1] = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    for (uint i = 0; i < samples; ++i, ++_n_processed_samples) {

      constexpr uint vec_size = vec_traits<float_x4>().size;

      // parameter smoothing.
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
          make_crange (_params.smooth_state.arr, vec_size, vec_size * j),
          vec_load<float_x4> (in),
          single_coeff_set_tag {});
      }
      // from now on access parameters without caring if they are smoothed or
      // not. This copy also has the aditional advantage of telling the compiler
      // that it can keep parameters on registers.
      all_parameters pars;
      static_assert (
        onepole_smoother::n_states == 1,
        "For the assignment below to work only the Z-1 state has to exist");
      *((smoothed_parameters*) &pars)   = _params.smooth_state.value;
      *((unsmoothed_parameters*) &pars) = _params.unsmoothed;

      // control rate refresh block
      if ((_n_processed_samples & _control_rate_mask) == 0) {
        _lfos[0].set_freq (pars.lfo_hz_final, _plugcontext->get_sample_rate());
        _lfos[1].set_freq (pars.lfo_hz_final, _plugcontext->get_sample_rate());

        auto stereo_ph = phase {pars.lfo_stereo, phase::degrees {}};
        if (pars.lfo_hz_final != 0.f) {
          // updating stereo phase diff.
          _lfos[1].set_phase (_lfos[0].get_phase() + stereo_ph);
        }
        else {
          auto start_ph = phase {pars.lfo_start_phase, phase::degrees {}};
          _lfos[0].set_phase (start_ph);
          _lfos[1].set_phase (start_ph + stereo_ph);
        }
        auto                 n_samples = _control_rate_mask + 1;
        std::array<float, 2> lfov;
        switch (pars.lfo_wave) {
        case 0:
          lfov[0] = _lfos[0].tick_sine (n_samples);
          lfov[1] = _lfos[1].tick_sine (n_samples);
          break;
        case 1:
          lfov[0] = _lfos[0].tick_triangle (n_samples);
          lfov[1] = _lfos[1].tick_triangle (n_samples);
          break;
        case 2:
          lfov[0] = _lfos[0].tick_sample_hold (n_samples);
          lfov[1] = _lfos[1].tick_sample_hold (n_samples);
          break;
        case 3:
          lfov[0] = _lfos[0].tick_filtered_noise (n_samples);
          lfov[1] = _lfos[1].tick_filtered_noise (n_samples);
          break;
        case 4:
          lfov[0] = _lfos[0].tick_trapezoid (0.3, n_samples);
          lfov[1] = _lfos[1].tick_trapezoid (0.3, n_samples);
          break;
        case 5:
          lfov[0] = _lfos[0].tick_square (n_samples);
          lfov[1] = _lfos[1].tick_square (n_samples);
          break;
        case 6:
          lfov[0] = _lfos[0].tick_saw (n_samples);
          lfov[1] = _lfos[1].tick_saw (n_samples);
          break;
        }

        auto   interp_stages = (double) (pars.n_allpasses / 2) - 1;
        double fconstant_even;
        double fconstant_odd;

        auto freq_hi = pars.freq_hi;
        auto freq_lo = pars.freq_lo;

        switch (pars.mode) {
        case m_exp:
        case m_exp_spread:
        case m_exp_alternate:
          fconstant_even = pow (freq_hi / freq_lo, 1. / interp_stages);
          fconstant_odd  = 1.;
          break;
        case m_exp_stereo:
          interp_stages *= 2;
          fconstant_even = pow (freq_hi / freq_lo, 1. / interp_stages);
          fconstant_odd  = fconstant_even;
          break;
        case m_lin:
        case m_lin_spread:
        case m_lin_alternate:
          fconstant_even = abs (freq_hi - freq_lo) / interp_stages;
          fconstant_odd  = 0.;
          break;
        case m_lin_stereo:
          interp_stages *= 2;
          fconstant_even = abs (freq_hi - freq_lo) / interp_stages;
          fconstant_odd  = fconstant_even;
          break;
        default:
          fconstant_even = 0;
          fconstant_odd  = 0.;
          break;
        }
        auto f = freq_lo;

        assert ((pars.n_allpasses % vec_size) == 0 && "bug!");
        for (uint s = 0; s < (pars.n_allpasses / vec_size); ++s) {
          float_x4 freqs {};
          float_x4 qs = vec_set<float_x4> (pars.q); // Q's don't oscilate (yet)

          // fill freqs
          for (uint ap = 0; ap < (vec_size / 2); ++ap) {
            constexpr float max_depth = 0.6;
            float           k;
            uint            stage = (s * 4) + ap;

            switch (pars.mode) {
            case m_exp:
            case m_lin:
            case m_exp_stereo:
            case m_lin_stereo:
              k = max_depth;
              break;
            case m_exp_spread:
            case m_lin_spread:
              k = (stage < (pars.n_allpasses / 2)) ? max_depth : -max_depth;
              break;
            case m_exp_alternate:
            case m_lin_alternate:
              k = ((stage / 2) & 1) ? max_depth : -max_depth;
              break;
            default:
              break;
            }

            freqs[ap * 2] = f + (lfov[0] * pars.lfo_depth * f * k);

            if (mode_is_exponential (pars.mode)) {
              f *= fconstant_odd;
            }
            else {
              f += fconstant_odd;
            }

            freqs[(ap * 2) + 1] = f + (lfov[1] * pars.lfo_depth * f * k);

            if (mode_is_exponential (pars.mode)) {
              f *= fconstant_even;
            }
            else {
              f += fconstant_even;
            }
          }
          andy::svf::reset_coeffs (
            get_allpass_group_coeffs (s),
            freqs,
            qs,
            _plugcontext->get_sample_rate(),
            allpass_tag {});
        }
      }

      // regular processing
      float_x4 out {
        chnls[0][i] + (_feedback_samples[0] * pars.feedback),
        chnls[1][i] + (_feedback_samples[1] * pars.feedback),
        chnls[0][i] + (_feedback_samples[0] * pars.feedback),
        chnls[1][i] + (_feedback_samples[1] * pars.feedback)};

      for (uint g = 0; g < (pars.n_allpasses / vec_size); ++g) {
        // as of now this processes both in parallel and in series.
        out = andy::svf::tick (
          get_allpass_group_coeffs (g), get_allpass_group_states (g), out);
      }
      assert (
        (pars.n_allpasses % vec_size) == 0
        && "there are unprocessed allpasses");

      double_x2 outs = {
        (double) out[0] + (double) out[2], (double) out[1] + (double) out[3]};
      outs *= -0.5;
      chnls[0][i] = outs[0];
      chnls[1][i] = outs[1];

      outs
        = mystran_dc_blocker::tick (_dc_block_coeffs, _dc_block_states, outs);

      _feedback_samples[0] = outs[0];
      _feedback_samples[1] = outs[1];
    }
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
    q_tag>;
  //----------------------------------------------------------------------------
private:
  void refresh_lfo_hz()
  {
    float hz = 0.f;
    if (_params.unsmoothed.lfo_eights != 0.f) {
      hz = (8. * _plugcontext->get_play_state().bpm)
        / (60. * _params.unsmoothed.lfo_eights);
    }
    hz += _params.unsmoothed.lfo_hz_user;
    _params.smooth_target.value.lfo_hz_final = hz;
  }
  //----------------------------------------------------------------------------
  crange<float> get_allpass_group_coeffs (uint n)
  {
    return make_crange (
      &_allpass[n][0], andy::svf::n_coeffs * vec_traits<float_x4>().size);
  }
  //----------------------------------------------------------------------------
  crange<float> get_allpass_group_states (uint n)
  {
    return make_crange (
      &_allpass[n][andy::svf::n_coeffs * vec_traits<float_x4>().size],
      andy::svf::n_states * vec_traits<float_x4>().size);
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    float lfo_eights;
    float lfo_hz_user;
    uint  lfo_wave;
    uint  n_allpasses;
    uint  mode;
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
  static constexpr uint n_channels = 2;

  std::array<float, n_channels> _feedback_samples;
  using simd_allpass_group = std::array<
    float,
    (andy::svf::n_coeffs + andy::svf::n_states) * vec_traits_t<float_x4>::size>;

  using dc_block_coeff_array
    = simd_array<double, mystran_dc_blocker::n_coeffs * n_channels, sse_bytes>;
  using dc_block_state_array
    = simd_array<double, mystran_dc_blocker::n_states * n_channels, sse_bytes>;

  parameter_values _params;
  // 16 groups of 4 filters (SIMD), 32 allpasses, interleaved for L and R = 16
  alignas (sse_bytes) std::array<simd_allpass_group, 16> _allpass;
  alignas (sse_bytes) dc_block_coeff_array _dc_block_coeffs;
  alignas (sse_bytes) dc_block_state_array _dc_block_states;

  std::array<lfo, 2> _lfos; // 0 = L, 1 = R
  uint               _n_processed_samples;
  uint               _control_rate_mask;

  float _lp_smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv