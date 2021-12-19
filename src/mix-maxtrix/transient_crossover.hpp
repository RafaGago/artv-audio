#pragma once

#include <algorithm>
#include <cmath>
#include <gcem.hpp>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// Mod of Saike's transience for transient splitting.
//------------------------------------------------------------------------------
struct mixmaxtrix_transient_crossover {
public:
  using V = vec<double, 1>;
  using T = double;

  static constexpr dsp_types dsp_type  = dsp_types::dynamics;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 2;
  //----------------------------------------------------------------------------
  struct transient_output {};

  void set (transient_output, int)
  {
    // dummy, to be used by processor
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (transient_output)
  {
    constexpr auto str_array = make_cstr_array (
      "bus2", "bus3", "bus4", "bus5", "bus6", "bus7", "bus8");
    constexpr uint n_future_choices = 16;

    return choice_param (1, str_array, n_future_choices);
  }
  //----------------------------------------------------------------------------
  struct gainsmoothing_tag {};

  void set (gainsmoothing_tag, float v)
  {
    if (v == _gainsmoothingp) {
      return;
    }
    _gainsmoothingp = v;
    slider();
  }

  static constexpr auto get_parameter (gainsmoothing_tag)
  {
    return float_param ("", 0.0, 1.0, 0.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};

  void set (mode_tag, uint v)
  {
    _mode = v;
    // slider();
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "L+R Peak",
        "L+R Squared",
        "L-Peak",
        "L-Squared",
        "R-Peak",
        "R-Squared",
        "L-R Peak",
        "L-R Squared"),
      8);
  }
  //----------------------------------------------------------------------------
  enum mode {
    mode_lr_peak,
    mode_lr_squared,
    mode_l_peak,
    mode_l_squared,
    mode_r_peak,
    mode_r_squared,
    mode_lmr_peak,
    mode_lmr_squared,
  };
  //----------------------------------------------------------------------------
  struct sattack_tag {};

  void set (sattack_tag, float v)
  {
    if (v == _sattackp) {
      return;
    }
    _sattackp = v;
    slider();
  }

  static constexpr auto get_parameter (sattack_tag)
  {
    return float_param ("", 0.0, 1.0, 0.5, 1e-09);
  }
  //----------------------------------------------------------------------------
  struct sdecay_tag {};

  void set (sdecay_tag, float v)
  {
    if (v == _sdecayp) {
      return;
    }
    _sdecayp = v;
    slider();
  }

  static constexpr auto get_parameter (sdecay_tag)
  {
    return float_param ("", 0.0, 1.0, 0.5, 1e-09);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    sattack_tag,
    sdecay_tag,
    mode_tag,
    gainsmoothing_tag,
    transient_output>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _sample_rate   = pc.get_sample_rate(); // TODO: iterate every param
    _gain_prev     = 0;
    _gain_atk_prev = 0;
    _gain_dec_prev = 0;
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {

    std::array<std::array<double, 3>, 32> in;

    static constexpr uint l      = 0;
    static constexpr uint r      = 1;
    static constexpr uint detect = 2;

    size_t done = 0;
    while (done < samples) {
      uint blocksize = std::min (samples - done, in.size());

      switch (_mode) {
      case mode_lr_peak:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = std::abs (in[i][l]) + std::abs (in[i][r]) * 0.5;
        }
        break;
      case mode_lr_squared:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = (in[i][l] * in[i][l] + in[i][r] * in[i][r]) * 0.5;
        }
        break;
      case mode_l_peak:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = std::abs (in[i][l]);
        }
        break;
      case mode_l_squared:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = in[i][l] * in[i][l];
        }
        break;
      case mode_r_peak:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = std::abs (in[i][r]);
        }
        break;
      case mode_r_squared:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = in[i][r] * in[i][r];
        }
        break;
      case mode_lmr_peak:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l]      = ins[0][done + i];
          in[i][r]      = ins[1][done + i];
          in[i][detect] = std::abs (in[i][l] - in[i][r]) * 0.5;
        }
        break;
      case mode_lmr_squared:
        for (uint i = 0; i < blocksize; ++i) {
          in[i][l] = ins[0][done + i];
          in[i][r] = ins[1][done + i];
          in[i][detect]
            = std::abs ((in[i][l] * in[i][l] - in[i][r] * in[i][r]) * 0.5);
        }
        break;
      }

      for (uint i = 0; i < blocksize; ++i) {
        double in_l = in[i][l];
        double in_r = in[i][r];

        double in_gain = 20. * std::log10 (std::max (0.001, in[i][detect]));

        double c_env = _followers.tick<flw_main> (vec_set<V> (in_gain))[0];
        double c_target_atk
          = _followers.tick<flw_attack> (vec_set<V> (in_gain))[0];
        double c_target_decay
          = _followers.tick<flw_decay> (vec_set<V> (in_gain))[0];

        double gain = _alpha_gain * _gain_prev + (1. - _alpha_gain) * c_env;
        double gain_atk
          = _alpha_gain * _gain_atk_prev + (1. - _alpha_gain) * c_target_atk;
        double gain_dec
          = _alpha_gain * _gain_dec_prev + (1. - _alpha_gain) * c_target_decay;

        _gain_prev     = gain;
        _gain_atk_prev = gain_atk;
        _gain_dec_prev = gain_dec;

        double mix = (gain_atk - gain_dec) / gain;
        mix        = std::clamp (mix, 0., 1.);

        double out_l = in_l * mix;
        double out_r = in_r * mix;

        outs[0][done + i] = in_l - out_l;
        outs[1][done + i] = in_r - out_r;
        outs[2][done + i] = out_l;
        outs[3][done + i] = out_r;
      }
      done += blocksize;
    };
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void slider()
  {
    static constexpr double follow_attack = 1.;

    static constexpr double min_attack = 2.;
    static constexpr double max_attack = 120.;
    static constexpr double min_decay  = 130.;
    static constexpr double max_decay  = 1000.;

    static constexpr double atk_beta  = gcem::log (min_attack);
    static constexpr double atk_alpha = gcem::log (max_attack) - atk_beta;
    static constexpr double dec_beta  = gcem::log (min_decay);
    static constexpr double dec_alpha = gcem::log (max_decay) - dec_beta;

    // code on the "block" part of the original
    auto attack = exp (atk_alpha * _sattackp + atk_beta) - 1.f;
    auto decay  = exp (dec_alpha * _sdecayp + dec_beta) - 1.f;
    // Max gain smoothing is 15 ms;_sample_rate
    _alpha_gain = exp (-1 / (.5 * .015 * _gainsmoothingp * _sample_rate));

    // notice, all the envelope followers on the original multiply the sample
    // rate by 0.5 mine dont.
    //
    // To pass to msec a factor of 0.001 is added too,
    constexpr float sr_correct = 0.5 * 0.001;
    auto            srate      = _sample_rate * sr_correct;

    /* 20 Hz => ~50 ms period */
    _followers.reset_coeffs<flw_main> (
      vec_set<V> (follow_attack), vec_set<V> (120.), srate);

    _followers.reset_coeffs<flw_attack> (
      vec_set<V> (attack), vec_set<V> (150.), srate);

    _followers.reset_coeffs<flw_decay> (
      vec_set<V> (follow_attack), vec_set<V> (decay), srate);
  }
  //----------------------------------------------------------------------------
  // ported parameters from the old generated JSFX, names come from there too
  float _gainsmoothingp = 0.0f;
  float _sattackp       = 0.f;
  float _sdecayp        = 0.f;
  uint  _mode           = 0;

  // variables set on the block function
  float _alpha_gain    = 0.;
  float _gain_prev     = 0.;
  float _gain_atk_prev = 0.;
  float _gain_dec_prev = 0.;

  enum followers { flw_main, flw_attack, flw_decay, n_flw };
  part_to_class<V, slew_limiter, n_flw> _followers;
  float                                 _sample_rate = 0;
};

} // namespace artv
