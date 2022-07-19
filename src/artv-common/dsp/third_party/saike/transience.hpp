// Ported from https://github.com/JoepVanlier/JSFX.git
// commit sha: 5d54a165b806b2773792623376f6a51cc957368b
#pragma once

#include <cmath>
#include <gcem.hpp>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
namespace artv { namespace saike {
//------------------------------------------------------------------------------
struct transience {
public:
  using V = vec<double, 1>;
  using T = double;

  static constexpr dsp_types dsp_type  = dsp_types::dynamics;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct compensate_tag {};

  void set (compensate_tag, float v)
  {
    if (v == _compensatep) {
      return;
    }
    _compensatep = v;
    slider();
  }

  static constexpr auto get_parameter (compensate_tag)
  {
    return float_param ("", -20.0, 20.0, 0.0, 0.0005);
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
    if (!!v == _squared_mode) {
      return;
    }
    _squared_mode = !!v;
    // slider();
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (0, make_cstr_array ("Peak", "Squared"));
  }
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
  struct strength_tag {};

  void set (strength_tag, float v)
  {
    if (v == _strengthp) {
      return;
    }
    _strengthp = v;
    slider();
  }

  static constexpr auto get_parameter (strength_tag)
  {
    return float_param ("", -1.0, 1.0, 0.0, 1e-09);
  }
  //----------------------------------------------------------------------------
  struct strength2_tag {};

  void set (strength2_tag, float v)
  {
    if (v == _strength2p) {
      return;
    }
    _strength2p = v;
    slider();
  }

  static constexpr auto get_parameter (strength2_tag)
  {
    return float_param ("", -1.0, 1.0, 0.0, 1e-09);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    sattack_tag,
    strength_tag,
    sdecay_tag,
    strength2_tag,
    mode_tag,
    gainsmoothing_tag,
    compensate_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { reset (pc.get_sample_rate()); }
  //----------------------------------------------------------------------------
  void reset (double sample_rate)
  {
    _sample_rate  = sample_rate;
    _db_gain_prev = 0;
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    for (uint i = 0; i < samples; ++i) {
      T gain     = get_gain (ins[0][i], ins[1][i]);
      outs[0][i] = (ins[0][i] * gain);
      outs[1][i] = (ins[1][i] * gain);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<std::array<T, 2>> io)
  {
    for (uint i = 0; i < io.size(); ++i) {
      T gain   = get_gain (io[i][0], io[i][1]);
      io[i][0] = (io[i][0] * gain);
      io[i][1] = (io[i][1] * gain);
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  T get_gain (T in_l, T in_r)
  {
    double in_gain;
    if (_squared_mode) {
      in_gain = in_l * in_l + in_r * in_r;
    }
    else {
      in_gain = std::abs (in_l) + std::abs (in_r);
    }
    in_gain *= 0.5;
    in_gain = 20. * std::log10 (std::max (0.001, in_gain));

    double c_env        = _followers.tick<flw_main> (vec_set<V> (in_gain))[0];
    double c_target_atk = _followers.tick<flw_attack> (vec_set<V> (in_gain))[0];
    double c_target_decay
      = _followers.tick<flw_decay> (vec_set<V> (in_gain))[0];

    /* Gain changes in dB space */
    double l_diff_atk  = _strengthp * (c_target_atk - c_env);
    double l_diff_rel  = _strength2p * (c_target_decay - c_env);
    double db_gain_cur = -l_diff_atk + l_diff_rel + _compensatep;
    double db_gain
      = _alpha_gain * _db_gain_prev + (1. - _alpha_gain) * db_gain_cur;
    _db_gain_prev = db_gain;

    /* Convert to linear */
    return exp ((T) (M_LN10 * (0.05 * db_gain)));
  }
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
  float _compensatep    = 0.0f;
  float _gainsmoothingp = 0.0f;
  float _sattackp       = 0.5f;
  float _sdecayp        = 0.5f;
  float _strengthp      = 0.0f;
  float _strength2p     = 0.0f;
  bool  _squared_mode   = false;

  // variables set on the block function
  float _alpha_gain   = 0.;
  float _db_gain_prev = 0.;

  enum followers { flw_main, flw_attack, flw_decay, n_flw };
  part_class_array<slew_limiter, V, n_flw> _followers;
  float                                    _sample_rate = 0;
};
}} // namespace artv::saike
//------------------------------------------------------------------------------
