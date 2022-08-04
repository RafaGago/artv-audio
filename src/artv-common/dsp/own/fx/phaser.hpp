#pragma once

#include <array>
#include <cstdint>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
class phaser {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::modulation;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
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
  void set (stages_mode_tag, int v)
  {
    bool lin                       = v < m_total_modes;
    _params.unsmoothed.lin_flo_mod = lin;
    _params.unsmoothed.mode        = v - (!lin * m_total_modes);
  }

  static constexpr auto get_parameter (stages_mode_tag)
  {
    return choice_param (
      8,
      make_cstr_array (
        "Lin Mod Exp F",
        "Lin Mod Exp F Spread",
        "Lin Mod Exp F Alternate",
        "Lin Mod Exp F Stereo",
        "Lin Mod Lin F",
        "Lin Mod Lin F Spread",
        "Lin Mod Lin F Alternate",
        "Lin Mod Lin F Stereo",
        "Exp Mod Exp F",
        "Exp Mod Exp F Spread",
        "Exp Mod Exp F Alternate",
        "Exp Mod Exp F Stereo",
        "Exp Mod Lin F",
        "Exp Mod Lin F Spread",
        "Exp Mod Lin F Alternate",
        "Exp Mod Lin F Stereo"),
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
    m_total_modes,
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
    _allpass.reset_states_cascade();
    memset (&_feedback_samples, 0, sizeof _feedback_samples);

    _lfos.reset();

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

    using x1_type = vec<decltype (_lp_smooth_coeff), 1>;

    onepole_smoother::reset_coeffs (
      make_crange (_lp_smooth_coeff).cast (x1_type {}),
      vec_set<x1_type> (1. / 0.01),
      pc.get_sample_rate());

    _dc_blocker.reset_states();
    _dc_blocker.reset_coeffs (vec_set<double_x2> (2.), pc.get_sample_rate());

    _n_processed_samples = 0;

    // Sample rates 44100 multiples update every 362.811us
    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _feedback_samples[0] = _feedback_samples[1] = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
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
          make_crange (_params.smooth_state.arr, vec_size, vec_size * j)
            .cast (float_x4 {}),
          vec_load<float_x4> (in));
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
        _lfos.set_freq (
          vec_set<2> (pars.lfo_hz_final), _plugcontext->get_sample_rate());
        auto stereo_ph
          = phase<1> {vec_set<1> (pars.lfo_stereo), phase<1>::degrees {}}
              .get_raw (0);
        if (pars.lfo_hz_final != 0.f) {
          // updating stereo phase diff.
          auto ph = _lfos.get_phase();
          ph.set_raw (ph.get_raw (0) + stereo_ph, 1);
          _lfos.set_phase (ph);
        }
        else {
          auto start_ph = phase<n_channels> {
            vec_set<n_channels> (pars.lfo_start_phase),
            phase<n_channels>::degrees {}};
          start_ph.set_raw (start_ph.get_raw (0) + stereo_ph, 1);
          _lfos.set_phase (start_ph);
        }
        auto          n_samples = _control_rate_mask + 1;
        vec<float, 2> lfov;
        switch (pars.lfo_wave) {
        case 0:
          lfov = _lfos.tick_sine (n_samples);
          break;
        case 1:
          lfov = _lfos.tick_triangle (n_samples);
          break;
        case 2:
          lfov = _lfos.tick_sample_hold (n_samples);
          break;
        case 3:
          lfov = _lfos.tick_filt_sample_and_hold (n_samples);
          break;
        case 4:
          lfov = _lfos.tick_trapezoid (vec_set<2> (0.3f), n_samples);
          break;
        case 5:
          lfov = _lfos.tick_square (n_samples);
          break;
        case 6:
          lfov = _lfos.tick_saw (n_samples);
          break;
        }

        auto interp_stages = (double) (pars.n_allpasses / 2) - 1;
        std::array<double, 2> fconstant;

        auto freq_hi = std::max (pars.freq_hi, pars.freq_lo);
        auto freq_lo = std::min (pars.freq_hi, pars.freq_lo);

        switch (pars.mode) {
        case m_exp:
        case m_exp_spread:
        case m_exp_alternate:
          fconstant[0] = pow (freq_hi / freq_lo, 1. / interp_stages);
          fconstant[1] = 1.;
          break;
        case m_exp_stereo:
          interp_stages *= 2;
          fconstant[0] = pow (freq_hi / freq_lo, 1. / interp_stages);
          fconstant[1] = fconstant[0];
          break;
        case m_lin:
        case m_lin_spread:
        case m_lin_alternate:
          fconstant[0] = (freq_hi - freq_lo) / interp_stages;
          fconstant[1] = 0.;
          break;
        case m_lin_stereo:
          interp_stages *= 2;
          fconstant[0] = (freq_hi - freq_lo) / interp_stages;
          fconstant[1] = fconstant[0];
          break;
        default:
          fconstant[0] = 0.;
          fconstant[1] = 0.;
          break;
        }
        auto f = freq_lo;

        assert ((pars.n_allpasses % vec_size) == 0 && "bug!");
        for (uint s = 0; s < (pars.n_allpasses / vec_size); ++s) {
          float_x4 freqs {};
          float_x4 qs = vec_set<float_x4> (pars.q); // Q's don't oscilate (yet)

          // fill freqs
          for (uint ap = 0; ap < (vec_size / 2); ++ap) {
            float max_depth = pars.lin_flo_mod ? 0.6 : 1.7;
            float k;
            uint  stage = (s * 4) + ap;

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

            std::array<decltype (f), 2> fs;
            for (uint iv = 0; iv < fs.size(); ++iv) {
              fs[iv] = f;
              if (mode_is_exponential (pars.mode)) {
                f *= fconstant[iv];
              }
              else {
                f += fconstant[iv];
              }
            }

            if (pars.lin_flo_mod) {
              freqs[ap * 2] = fs[0] + (lfov[0] * pars.lfo_depth * fs[0] * k);
              freqs[ap * 2 + 1]
                = fs[1] + (lfov[1] * pars.lfo_depth * fs[1] * k);
            }
            else {
              freqs[ap * 2]     = fs[0] * exp (lfov[0] * pars.lfo_depth * k);
              freqs[ap * 2 + 1] = fs[1] * exp (lfov[1] * pars.lfo_depth * k);
            }
          }
          freqs = vec_min (20000.f, freqs);
          _allpass.reset_coeffs_on_idx (
            s, freqs, qs, _plugcontext->get_sample_rate());
        }
      }

      // regular processing
      float_x4 out {
        ins[0][i] + (_feedback_samples[0] * pars.feedback),
        ins[1][i] + (_feedback_samples[1] * pars.feedback),
        ins[0][i] + (_feedback_samples[0] * pars.feedback),
        ins[1][i] + (_feedback_samples[1] * pars.feedback)};

      for (uint g = 0; g < (pars.n_allpasses / vec_size); ++g) {
        // as of now this processes both in parallel and in series.
        out = _allpass.tick_on_idx (g, out);
      }
      assert (
        (pars.n_allpasses % vec_size) == 0
        && "there are unprocessed allpasses");

      double_x2 outx2 = {
        (double) out[0] + (double) out[2], (double) out[1] + (double) out[3]};
      outx2 *= -0.5;
      outs[0][i] = outx2[0];
      outs[1][i] = outx2[1];

      outx2 = _dc_blocker.tick (outx2);

      _feedback_samples[0] = outx2[0];
      _feedback_samples[1] = outx2[1];
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
  struct unsmoothed_parameters {
    float lfo_eights;
    float lfo_hz_user;
    uint  lfo_wave;
    uint  n_allpasses;
    uint  mode;
    bool  lin_flo_mod;
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
  parameter_values              _params;

  part_class_array<andy::svf_allpass, float_x4, 16> _allpass;
  part_class_array<mystran_dc_blocker, double_x2>   _dc_blocker;

  lfo<n_channels> _lfos; // 0 = L, 1 = R
  uint            _n_processed_samples;
  uint            _control_rate_mask;

  float _lp_smooth_coeff;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
