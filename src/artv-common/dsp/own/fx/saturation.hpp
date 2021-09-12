#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"
#include "artv-common/dsp/own/parts/filters/dc_blocker.hpp"
#include "artv-common/dsp/own/parts/filters/moving_average.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/misc/interpolators.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/compand1.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/hardclip.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/pow2.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/sqrt.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/sqrt_sigmoid.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/sqrt_sin.hpp"
#include "artv-common/dsp/own/parts/waveshapers/adaa/tanh.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

// A pretty naive but extreme saturation unit. Doing a subltler less featured
// variant based on Chebyshev polynomials could be fun.

// TODO: -parameter smoothing(?) (Gains, freq and filter amount(?).
namespace artv {

//------------------------------------------------------------------------------
class saturation {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::waveshaper;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct type_tag {};

  void set (type_tag, int v) { _p.type = (decltype (_p.type)) v; }

  static constexpr auto get_parameter (type_tag)
  {
    return choice_param (
      0, make_cstr_array ("Hardclip", "Tanh", "Sqrt", "SqrtSin"), 40);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};

  void set (mode_tag, int v)
  {
    _p.mode = (decltype (_p.mode)) v;
    if (_p.mode != _p.mode_prev) {
      uint delcomp;
      switch (_p.mode) {
      case mode_no_aa:
        delcomp = 0;
        break;
      case mode_normal:
        delcomp = 1;
        break;
      case mode_compand_1a:
      case mode_compand_sqrt:
        delcomp = 3;
        break;
      default:
        // assert (false);
        delcomp = 0;
        break;
      }
      _plugcontext->set_delay_compensation (delcomp);
    }
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "Normal no AA", "Normal", "Compand High", "Compand Squash"),
      40);
  }
  //----------------------------------------------------------------------------
  struct saturated_out_tag {};

  void set (saturated_out_tag, float v)
  {
    _sparams.set (db_to_gain (v), sm_out);
  }

  static constexpr auto get_parameter (saturated_out_tag)
  {
    return float_param ("dB", -20.0, 20., 0.0, 0.25, 0.6, true);
  }
  //----------------------------------------------------------------------------
  struct trim_tag {};

  void set (trim_tag, float v) { _sparams.set (db_to_gain (v), sm_trim); }

  static constexpr auto get_parameter (trim_tag)
  {
    return float_param ("dB", -20.0, 20., 0.0, 0.25, 0.6, true);
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};

  void set (drive_tag, float v)
  {
    _p.drive = db_to_gain (v);
    update_drive();
  }

  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -30.0, 30, 0.0, 0.25, 0.5, true);
  }
  //----------------------------------------------------------------------------
  struct drive_balance_tag {};

  void set (drive_balance_tag, float v)
  {
    _p.drive_bal = (v * 0.7 * 0.01) + 1.;
    update_drive();
  }

  static constexpr auto get_parameter (drive_balance_tag)
  {
    return float_param ("%", -100.0, 100., 0.0, 0.25, 0.6, true);
  }
  //----------------------------------------------------------------------------
  static constexpr float lo_cut_min_hz = 20.;

  struct lo_cut_tag {};

  void set (lo_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    if (v == _p.crossv_hz[lo_crossv]) {
      return;
    }
    _p.crossv_hz[lo_crossv]      = v;
    _p.crossv_enabled[lo_crossv] = (v > lo_cut_min_hz);
    update_crossover (lo_crossv, !_p.crossv_enabled[lo_crossv]);
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
    if (v == _p.crossv_hz[hi_crossv]) {
      return;
    }
    _p.crossv_hz[hi_crossv]      = v;
    _p.crossv_enabled[hi_crossv] = (v < hi_cut_max_hz);
    update_crossover (hi_crossv, !_p.crossv_enabled[hi_crossv]);
  }

  static constexpr auto get_parameter (hi_cut_tag)
  {
    return frequency_parameter (30., hi_cut_max_hz, hi_cut_max_hz);
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_crossv_order = 8;

  struct lo_mode_tag {};

  void set (lo_mode_tag, int v)
  {
    int  order = v + 1;
    bool is_lp = order <= max_crossv_order;
    order -= is_lp ? 0 : max_crossv_order;

    if (
      order == _p.crossv_order[lo_crossv]
      && is_lp == _p.crossv_is_lp[lo_crossv]) {
      return;
    }
    _p.crossv_order[lo_crossv] = order;
    _p.crossv_is_lp[lo_crossv] = is_lp;

    if (_p.crossv_enabled[lo_crossv]) {
      update_crossover (lo_crossv, true);
    }
  }

  static constexpr auto get_parameter (lo_mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
        "6dB/Oct",
        "12dB/Oct",
        "18dB/Oct",
        "24dB/Oct",
        "30dB/Oct",
        "36dB/Oct",
        "42dB/Oct",
        "48dB/Oct",
        "drop 6dB/Oct",
        "drop 12dB/Oct",
        "drop 18dB/Oct",
        "drop 24dB/Oct",
        "drop 30dB/Oct",
        "drop 36dB/Oct",
        "drop 42dB/Oct",
        "drop 48dB/Oct"),
      32);
  }
  //----------------------------------------------------------------------------
  struct hi_mode_tag {};

  void set (hi_mode_tag, int v)
  {
    int  order = v + 1;
    bool is_lp = (order > max_crossv_order);
    order -= is_lp ? max_crossv_order : 0;

    if (
      order == _p.crossv_order[hi_crossv]
      && is_lp == _p.crossv_is_lp[hi_crossv]) {
      return;
    }
    _p.crossv_order[hi_crossv] = order;
    _p.crossv_is_lp[hi_crossv] = is_lp;

    if (_p.crossv_enabled[hi_crossv]) {
      update_crossover (hi_crossv, true);
    }
  }

  static constexpr auto get_parameter (hi_mode_tag)
  {
    return get_parameter (lo_mode_tag {});
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};
  void set (feedback_tag, float v)
  {
    bool neg = v < 0.f;
    v *= 0.01;
    v = sqrtf (fabs (v));
    v *= 0.9f;
    v = neg ? -v : v;
    _sparams.set (v, sm_feedback);
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    return float_param ("%", -100., 100, 0., 0.1, 0.5, true);
  }
  //----------------------------------------------------------------------------
  struct emphasis_freq_tag {};

  void set (emphasis_freq_tag, float v)
  {
    _sparams.set (midi_note_to_hz (v), sm_emphasis_freq);
  }

  static constexpr auto get_parameter (emphasis_freq_tag)
  {
    return frequency_parameter (lo_cut_min_hz, 6000., 200.);
  }
  //----------------------------------------------------------------------------
  struct emphasis_amount_tag {};

  void set (emphasis_amount_tag, float v)
  {
    _sparams.set (midi_note_to_hz (v), sm_emphasis_amt);
  }

  static constexpr auto get_parameter (emphasis_amount_tag)
  {
    return float_param ("dB", -30.0, 30, 0.0, 0.25, 0.4, true);
  }
  //----------------------------------------------------------------------------
  struct emphasis_q_tag {};

  void set (emphasis_q_tag, float v) { _sparams.set (v, sm_emphasis_q); }

  static constexpr auto get_parameter (emphasis_q_tag)
  {
    return float_param ("", 0.01, 10., 0.5, 0.01, 0.6);
  }
  //----------------------------------------------------------------------------
  struct envfollow_attack_tag {};

  void set (envfollow_attack_tag, float v)
  {
    v *= 0.001f; // to seconds
    if (v != _p.ef_attack) {
      _p.ef_attack = v;
      update_envelope_follower();
    }
  }

  static constexpr auto get_parameter (envfollow_attack_tag)
  {
    return float_param ("ms", 1, 180., 20., 0.1);
  }
  //----------------------------------------------------------------------------
  struct envfollow_mode_tag {};

  void set (envfollow_mode_tag, int v) { _p.ef_mode = v; }

  static constexpr auto get_parameter (envfollow_mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("sqrt(exp)", "exp", "exp^2", "exp^3"), 20);
  }
  //----------------------------------------------------------------------------
  struct envfollow_release_tag {};

  void set (envfollow_release_tag, float v)
  {
    v *= 0.001f; // to seconds
    if (v != _p.ef_release) {
      _p.ef_release = v;
      update_envelope_follower();
    }
  }

  static constexpr auto get_parameter (envfollow_release_tag)
  {
    return float_param ("ms", 40., 1000., 150., 1.);
  }
  //----------------------------------------------------------------------------
  struct envfollow_sensitivity_tag {};

  void set (envfollow_sensitivity_tag, float v) { _p.ef_gain = db_to_gain (v); }

  static constexpr auto get_parameter (envfollow_sensitivity_tag)
  {
    return float_param ("dB", -20., 20., 0.);
  }
  //----------------------------------------------------------------------------
  static constexpr float envfollow_max_db = 30.;

  struct envfollow_to_drive_tag {};

  void set (envfollow_to_drive_tag, float v)
  {
    // -1 because the signal of the envelope follower has "1" added.
    constexpr auto max_db = (constexpr_db_to_gain (envfollow_max_db) - 1.f);

    auto efdrv = sqrt (abs (v) * 0.01); // max 1
    efdrv *= sgn_no_zero (v, -max_db, max_db);
    _sparams.set (efdrv, sm_ef_to_drive);
  }

  static constexpr auto get_parameter (envfollow_to_drive_tag)
  {
    return float_param ("%", -100., 100., 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct envfollow_to_emphasis_freq_tag {};

  void set (envfollow_to_emphasis_freq_tag, float v)
  {
    // max = 2, 2 octaves up and down
    _sparams.set (v * (0.01 * 2.), sm_ef_to_emphasis_freq);
  }

  static constexpr auto get_parameter (envfollow_to_emphasis_freq_tag)
  {
    return float_param ("%", -100., 100., 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct envfollow_to_emphasis_amount_tag {};

  void set (envfollow_to_emphasis_amount_tag, float v)
  {
    auto amt = v * 0.01f; // max 1
    amt      = sqrt (fabs (amt));
    amt *= sgn_no_zero (v, -30.f, 30.f);
    _sparams.set (amt, sm_ef_to_emphasis_amt);
  }

  static constexpr auto get_parameter (envfollow_to_emphasis_amount_tag)
  {
    return float_param ("%", -100., 100., 0., 0.01);
  }
  //----------------------------------------------------------------------------
  struct envfollow_to_dc_tag {};

  void set (envfollow_to_dc_tag, float v)
  {
    auto amt = sqrt (fabs (v * 0.01));
    amt *= sgn_no_zero (v);
    _sparams.set (amt, sm_ef_to_dc);
  }

  static constexpr auto get_parameter (envfollow_to_dc_tag)
  {
    return float_param ("%", -100., 100., 0., 0.01);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    type_tag,
    mode_tag,
    emphasis_amount_tag,
    emphasis_freq_tag,
    emphasis_q_tag,
    saturated_out_tag,
    drive_balance_tag,
    drive_tag,
    trim_tag,
    lo_cut_tag,
    hi_cut_tag,
    envfollow_attack_tag,
    envfollow_release_tag,
    envfollow_sensitivity_tag,
    envfollow_mode_tag,
    envfollow_to_drive_tag,
    envfollow_to_emphasis_freq_tag,
    envfollow_to_emphasis_amount_tag,
    envfollow_to_dc_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;

    memset (&_envfollow_states, 0, sizeof _envfollow_states);
    memset (&_compressor_states, 0, sizeof _compressor_states);
    memset (&_expander_states, 0, sizeof _expander_states);
    memset (&_wvsh_states, 0, sizeof _wvsh_states);
    memset (&_crossv_states, 0, sizeof _crossv_states);
    memset (&_dc_block_states, 0, sizeof _dc_block_states);
    memset (&_crossv_coeffs, 0, sizeof _crossv_coeffs);
    memset (&_pre_emphasis_states, 0, sizeof _pre_emphasis_states);
    memset (&_pre_emphasis_coeffs, 0, sizeof _pre_emphasis_coeffs);
    memset (&_post_emphasis_states, 0, sizeof _post_emphasis_states);
    memset (&_post_emphasis_coeffs, 0, sizeof _post_emphasis_coeffs);

    adaa::fix_eq_and_delay_coeff_initialization<adaa_order>::reset_coeffs<
      double_x2> (_adaa_fix_eq_delay_coeffs);

    // this DC blocker at a very low frequency is critical for the sound of
    // the companded modes to be acceptable. Unfortunately it causes DC itself
    // when the sound is muted. I didn't find a solution.
    mystran_dc_blocker::reset_coeffs (
      _dc_block_coeffs, make_vec_x1 (1.), pc.get_sample_rate());

    _p = params {};

    _sparams.reset ((float) pc.get_sample_rate());

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _p.type_prev = sat_type_count; // force initial click removal routine.

    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _n_processed_samples = 0;
    _sat_prev = _crossv_prev[lo_crossv] = _crossv_prev[hi_crossv]
      = double_x2 {};

    update_emphasis();
    _sparams.set_to_target();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    static constexpr double pow2compand_clip = 40.;
    params                  p                = _p;

    if (unlikely (p.type_prev != p.type || p.mode_prev != p.mode)) {
      memset (&_dc_block_states, 0, sizeof _dc_block_states);
      memset (&_compressor_states, 0, sizeof _compressor_states);
      memset (&_expander_states, 0, sizeof _expander_states);
      // some waveshapers will create peaks, as the integral on 0 might not be
      // 0. Running them for some samples of silence to initialize. This avoids
      // too having to run the "reset_states" functions on the waveshapers.
      static constexpr uint n_samples = 32;
      _p.type_prev                    = _p.type;
      _p.mode_prev                    = _p.mode;
      std::array<T, 2> value_increment {
        (T) (ins[0][0] * (1. / n_samples)), (T) (ins[1][0] * (1. / n_samples))};
      std::array<T, n_samples * 2> in {};
      for (uint i = 1; i < n_samples; ++i) {
        in[i] = in[i - 1] + value_increment[0];
      }
      for (uint i = n_samples + 1; i < (n_samples * 2); ++i) {
        in[i] = in[i - 1] + value_increment[1];
      }
      std::array<T*, 2> io  = {&in[0], &in[n_samples]};
      auto              cio = array_const_cast<T const*> (io);
      process<T> (io, cio, n_samples);
    }

    static constexpr uint                 max_block_size    = 32;
    uint                                  samples_processed = 0;
    std::array<double_x2, max_block_size> sat, lo, hi;

    while (samples_processed < block_samples) {
      uint subblock_size
        = std::min (block_samples - samples_processed, max_block_size);

      // Input DC blocking
      for (uint i = 0; i < subblock_size; ++i) {
        sat[i][0] = ins[0][samples_processed + i];
        sat[i][1] = ins[1][samples_processed + i];

        sat[i] = mystran_dc_blocker::tick (
          _dc_block_coeffs,
          _dc_block_states[dc_block_in],
          sat[i],
          single_coeff_set_tag {});
      }

      // pre process/Companding + data initialization block
      //
      // Using "compand_1b_aa" first required extremely good DC blocking
      // unfortunately. The sound easily broke.
      switch (p.mode) {
      case mode_compand_1a:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = compand_1a_aa::tick (
            _adaa_fix_eq_delay_coeffs, _compressor_states, sat[i]);
        }
        break;
      case mode_compand_sqrt:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = sqrt_aa::tick (
            _adaa_fix_eq_delay_coeffs, _compressor_states, sat[i]);
        }
        break;
      case mode_no_aa:
      case mode_normal:
      default:
        break;
      }

      // Crossover section. Notice that this is a crossover that sums to flat
      // frequency and phase response but has ripples on the passed band.
      for (uint i = 0; i < subblock_size; ++i) {
        if (_p.crossv_enabled[lo_crossv]) {
          lo[i] = butterworth_any_order::tick (
            _crossv_coeffs[lo_crossv],
            _crossv_states[lo_crossv],
            sat[i],
            p.crossv_order[lo_crossv]);
          sat[i] = p.crossv_is_lp[lo_crossv] ? sat[i] - lo[i] : lo[i];
          lo[i]  = p.crossv_is_lp[lo_crossv] ? lo[i] : double_x2 {};
        }
        else {
          lo[i] = double_x2 {};
        }

        if (_p.crossv_enabled[hi_crossv]) {
          hi[i] = butterworth_any_order::tick (
            _crossv_coeffs[hi_crossv],
            _crossv_states[hi_crossv],
            sat[i],
            p.crossv_order[hi_crossv]);
          sat[i] = !p.crossv_is_lp[hi_crossv] ? sat[i] - hi[i] : hi[i];
          hi[i]  = !p.crossv_is_lp[hi_crossv] ? hi[i] : double_x2 {};
        }
        else {
          hi[i] = double_x2 {};
        }
      }

      // Main block with audio-rate modulations. Done sample-wise for
      // simplicity
      for (uint i = 0; i < subblock_size; ++i, ++_n_processed_samples) {
        if (adaa_order == 1 && waveshaper_type_is_adaa (p.mode)) {
          // One sample delay for hi and lo, as the ADAA chain will add 1
          // sample delay, we mix with the previous crossover outputs.
          std::swap (_crossv_prev[lo_crossv], lo[i]);
          std::swap (_crossv_prev[hi_crossv], hi[i]);
        }

        // smoothing some params
        _sparams.tick();

        double_x2 drive {_sparams.get (sm_drive_l), _sparams.get (sm_drive_r)};
        double_x2 inv_drive = 1. / drive;
        inv_drive *= _sparams.get (sm_trim);

        // Envelope follower/modulation
        double_x2 follow
          = slew_limiter::tick (_envfollow_coeffs, _envfollow_states, sat[i]);
        // Make unipolar and clip at 1.
        follow = vec_min (vec_abs (follow) * p.ef_gain, 1.);

        switch (_p.ef_mode) {
        case ef_sqrt:
          follow = vec_sqrt (follow);
          break;
        case ef_linear:
          break;
        case ef_pow2:
          follow *= follow;
          break;
        case ef_pow3:
          follow *= follow * follow;
          break;
        default:
          assert (false);
          break;
        }

        if ((_n_processed_samples & _control_rate_mask) == 0) {
          // expensive stuff at lower than audio rate
          update_emphasis (
            vec_exp2 (follow * _sparams.get (sm_ef_to_emphasis_freq))
              * _sparams.get (sm_emphasis_freq),
            follow * _sparams.get (sm_ef_to_emphasis_amt));
        }

        float ef_to_drive = _sparams.get (sm_ef_to_drive);

        if (ef_to_drive != 0.f) {
          // gain modulation. Audio rate.
          double_x2 gainfollow
            = (follow * abs (ef_to_drive)) + 1.; // 1 to N range

          if (ef_to_drive > 0.f) {
            drive *= gainfollow;
            inv_drive /= gainfollow;
          }
          else {
            // inverse modulation
            inv_drive *= gainfollow;
            drive /= gainfollow;
          }
        }

        sat[i] *= drive;

        // Feedback section
        static constexpr auto fb_att = constexpr_db_to_gain (
          -constexpr_db_to_gain (envfollow_max_db) - 1.);

        auto feedback = _sat_prev * _sparams.get (sm_feedback);

        double_x2 feedback_follow = follow * ef_to_drive * fb_att;

        // inverted channels on feedback_follow, just for fun...
        feedback_follow = vec_shuffle (feedback_follow, feedback_follow, 1, 0);
        feedback += feedback * feedback_follow;

        sat[i] += feedback;
        auto dcmod = (0.5 * follow * _sparams.get (sm_ef_to_dc));
        sat[i] += dcmod;

        // pre emphasis
        sat[i] = andy::svf::tick (
          _pre_emphasis_coeffs, _pre_emphasis_states, sat[i]);

        // waveshaping
        auto wsh_type
          = p.type + (sat_type_count * (uint) waveshaper_type_is_adaa (p.mode));

        switch (wsh_type) {
        case sat_hardclip:
          sat[i] = wavesh_tick<hardclip_adaa<0>> (sat[i]);
          break;
        case sat_tanh:
          sat[i] = wavesh_tick<tanh_adaa<0>> (sat[i]);
          break;
        case sat_sqrt:
          sat[i] = wavesh_tick<sqrt_sigmoid_adaa<0>> (sat[i]);
          break;
        case sat_sqrt_sin:
          sat[i] = wavesh_tick<sqrt_sin_sigmoid_adaa<0>> (sat[i]);
          break;
        case sat_hardclip_adaa:
          sat[i] = wavesh_tick<hardclip_aa> (sat[i]);
          break;
        case sat_tanh_adaa:
          sat[i] = wavesh_tick<tanh_aa> (sat[i]);
          break;
        case sat_sqrt_adaa:
          sat[i] = wavesh_tick<sqrt_sigmoid_aa> (sat[i]);
          break;
        case sat_sqrt_sin_adaa:
          sat[i] = wavesh_tick<sqrt_sin_sigmoid_aa> (sat[i]);
          break;
        default:
          assert (false);
          break;
        }

        // post emphasis
        sat[i] = andy::svf::tick (
          _post_emphasis_coeffs, _post_emphasis_states, sat[i]);

        // not a DC blocker, but reduces DC on the FB path. Notice that adding
        // a third DC blocker here had negative effect on the sound when using
        // "mode_compand_1b".
        sat[i] -= dcmod;
        _sat_prev = sat[i];

        // gain and crossover join
        sat[i] *= inv_drive * _sparams.get (sm_out);
        sat[i] += lo[i] + hi[i];
      }

      // Post DC-blocker.
      for (uint i = 0; i < subblock_size; ++i) {
        sat[i] = mystran_dc_blocker::tick (
          _dc_block_coeffs,
          _dc_block_states[dc_block_out],
          sat[i],
          single_coeff_set_tag {});

        outs[0][samples_processed + i] = sat[i][0];
        outs[1][samples_processed + i] = sat[i][1];
      }

      // post process / restore
      switch (p.mode) {
      case mode_compand_1a:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = compand_1b_aa::tick (
            _adaa_fix_eq_delay_coeffs, _expander_states, sat[i]);
        }
        break;
      case mode_compand_sqrt:
        for (uint i = 0; i < subblock_size; ++i) {
          // sat[i] = vec_clamp (sat[i], -pow2compand_clip, pow2compand_clip);
          sat[i] = pow2_aa::tick (
            _adaa_fix_eq_delay_coeffs, _expander_states, sat[i]);
        }
        break;
      case mode_no_aa:
      case mode_normal:
      default:
        break;
      }
      samples_processed += subblock_size;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  enum smoothed_state {
    sm_drive_l,
    sm_drive_r,
    sm_feedback,
    sm_out,

    sm_trim,
    sm_emphasis_freq,
    sm_emphasis_amt,
    sm_emphasis_q,

    sm_ef_to_drive,
    sm_ef_to_emphasis_freq,
    sm_ef_to_emphasis_amt,
    sm_ef_to_dc,

    sm_count,
  };
  //----------------------------------------------------------------------------
  void update_crossover (uint idx, bool reset_states)
  {
    butterworth_any_order::reset_coeffs (
      _crossv_coeffs[idx],
      double_x2 {_p.crossv_hz[idx], _p.crossv_hz[idx] * 1.009},
      _plugcontext->get_sample_rate(),
      _p.crossv_order[idx],
      _p.crossv_is_lp[idx]);

    if (reset_states) {
      _crossv_prev[idx] = double_x2 {};
      memset (&_crossv_states[idx], 0, sizeof _crossv_states[idx]);
    }
  }
  //----------------------------------------------------------------------------
  void update_drive()
  {
    _sparams.set (_p.drive * (2. - _p.drive_bal), sm_drive_l);
    _sparams.set (_p.drive * (_p.drive_bal), sm_drive_r);
  }
  //----------------------------------------------------------------------------
  static constexpr bool waveshaper_type_is_adaa (char mode)
  {
    return mode != mode_no_aa;
  }
  //----------------------------------------------------------------------------
  static constexpr bool has_compander (char mode)
  {
    return mode >= mode_compander_first;
  }
  //----------------------------------------------------------------------------
  void update_emphasis (
    double_x2 freq_offset = double_x2 {},
    double_x2 amt_offset  = double_x2 {})
  {
    double_x2 f  = vec_set<double_x2> (_sparams.get (sm_emphasis_freq));
    double_x2 q  = vec_set<double_x2> (_sparams.get (sm_emphasis_q));
    double_x2 db = vec_set<double_x2> (_sparams.get (sm_emphasis_amt));

    f += freq_offset;
    db += amt_offset;

    andy::svf::reset_coeffs (
      _pre_emphasis_coeffs,
      f,
      q,
      db,
      _plugcontext->get_sample_rate(),
      bell_tag {});

    andy::svf::reset_coeffs (
      _post_emphasis_coeffs,
      f,
      q,
      -db,
      _plugcontext->get_sample_rate(),
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void update_envelope_follower()
  {
    slew_limiter::reset_coeffs (
      _envfollow_coeffs,
      vec_set<double_x2> (_p.ef_attack),
      vec_set<double_x2> (_p.ef_release),
      _plugcontext->get_sample_rate());
  }
  //----------------------------------------------------------------------------
  template <class wsh>
  double_x2 wavesh_tick (double_x2 x)
  {
    return wsh::template tick (_adaa_fix_eq_delay_coeffs, _wvsh_states, x);
  }
  //----------------------------------------------------------------------------
  enum sat_type {

    sat_hardclip,
    sat_tanh,
    sat_sqrt,
    sat_sqrt_sin,
    // adaa
    sat_hardclip_adaa,
    sat_tanh_adaa,
    sat_sqrt_adaa,
    sat_sqrt_sin_adaa,
    sat_type_count = (sat_sqrt_sin_adaa + 1) / 2,
  };

  enum mode_type {
    mode_no_aa,
    mode_normal,
    mode_compander_first,
    mode_compand_1a = mode_compander_first,
    mode_compand_sqrt,
    mode_count,
  };

  enum ef_type {
    ef_sqrt,
    ef_linear,
    ef_pow2,
    ef_pow3,
  };

  enum ef_block_type {
    dc_block_in,
    dc_block_out,
    dc_block_count,
  };

  enum crossv_indexes { lo_crossv, hi_crossv, n_crossv };

  struct params {
    params()
    {
      for (uint i = 0; i < n_crossv; ++i) {
        crossv_hz[i]      = -1.f;
        crossv_order[i]   = 1;
        crossv_enabled[i] = false;
      }
      crossv_is_lp[lo_crossv] = true;
      crossv_is_lp[hi_crossv] = false;
    }

    std::array<int, n_crossv>   crossv_order;
    std::array<float, n_crossv> crossv_hz;
    float                       drive      = 1.f;
    float                       drive_bal  = 0.f;
    float                       ef_attack  = 0.;
    float                       ef_release = 0.f;
    float                       ef_gain    = 1.f;
    char                        type       = 0;
    char                        type_prev  = 1;
    char                        mode       = 1;
    char                        mode_prev  = 0;
    char                        ef_mode    = ef_linear;
    std::array<bool, n_crossv>  crossv_is_lp;
    std::array<bool, n_crossv>  crossv_enabled;
  };

  template <class T>
  using to_n_states = std::integral_constant<int, T::n_states>;

  template <class T>
  using to_n_coeffs = std::integral_constant<int, T::n_coeffs>;

  static constexpr uint adaa_order = 1;

  using sqrt_sigmoid_aa = adaa::fix_eq_and_delay<adaa_order, sqrt_sigmoid_adaa>;
  using tanh_aa         = adaa::fix_eq_and_delay<adaa_order, tanh_adaa>;
  using hardclip_aa     = adaa::fix_eq_and_delay<adaa_order, hardclip_adaa>;
  using sqrt_sin_sigmoid_aa
    = adaa::fix_eq_and_delay<adaa_order, sqrt_sin_sigmoid_adaa>;

  using shapers
    = mp_list<sqrt_sigmoid_aa, tanh_aa, hardclip_aa, sqrt_sin_sigmoid_aa>;

  static constexpr uint wsh_max_states = mp11::mp_max_element<
    mp11::mp_transform<to_n_states, shapers>,
    mp11::mp_less>::value;

  using sqrt_aa       = adaa::fix_eq_and_delay<adaa_order, sqrt_adaa>;
  using pow2_aa       = adaa::fix_eq_and_delay<adaa_order, pow2_adaa>;
  using compand_1a_aa = adaa::fix_eq_and_delay<adaa_order, compand_1a_adaa>;
  using compand_1b_aa = adaa::fix_eq_and_delay<adaa_order, compand_1b_adaa>;

  using companders = mp_list<sqrt_aa, pow2_aa, compand_1a_aa, compand_1b_aa>;

  static constexpr uint compander_max_states = mp11::mp_max_element<
    mp11::mp_transform<to_n_states, companders>,
    mp11::mp_less>::value;

  using smoother = onepole_smoother;

  static constexpr uint n_channels = 2;
  using wsh_state_array
    = simd_array<double, wsh_max_states * n_channels, sse_bytes>;
  using compander_state_array
    = simd_array<double, compander_max_states * n_channels, sse_bytes>;
  using fix_eq_and_delay_coeff_array = simd_array<
    double,
    adaa::fix_eq_and_delay_coeff_initialization<adaa_order>::n_coeffs
      * n_channels,
    sse_bytes>;

  using crossv_coeff_array = simd_array<
    double,
    butterworth<max_crossv_order>::n_coeffs * n_channels,
    sse_bytes>;
  using crossv_state_array = simd_array<
    double,
    butterworth<max_crossv_order>::n_states * n_channels,
    sse_bytes>;

  using emphasis_coeff_array
    = simd_array<double, andy::svf::n_coeffs * n_channels, sse_bytes>;
  using emphasis_state_array
    = simd_array<double, andy::svf::n_states * n_channels, sse_bytes>;

  using envfollow_coeff_array
    = simd_array<double, slew_limiter::n_coeffs * n_channels, sse_bytes>;
  using envfollow_state_array
    = simd_array<double, slew_limiter::n_states * n_channels, sse_bytes>;

  using dc_block_coeff_array
    = simd_array<double, mystran_dc_blocker::n_coeffs * n_channels, sse_bytes>;
  using dc_block_state_array
    = simd_array<double, mystran_dc_blocker::n_states * n_channels, sse_bytes>;

  // all arrays are multiples of the simd size, no need to alignas on
  // everything.
  params _p;
  alignas (sse_bytes) envfollow_coeff_array _envfollow_coeffs;
  envfollow_state_array                    _envfollow_states;
  fix_eq_and_delay_coeff_array             _adaa_fix_eq_delay_coeffs;
  compander_state_array                    _compressor_states;
  std::array<crossv_coeff_array, n_crossv> _crossv_coeffs;
  std::array<crossv_state_array, n_crossv> _crossv_states;
  value_smoother<float, sm_count>          _sparams;

  alignas (sse_bytes) emphasis_coeff_array _pre_emphasis_coeffs;
  emphasis_state_array                             _pre_emphasis_states;
  wsh_state_array                                  _wvsh_states;
  dc_block_coeff_array                             _dc_block_coeffs;
  std::array<dc_block_state_array, dc_block_count> _dc_block_states;
  emphasis_coeff_array                             _post_emphasis_coeffs;
  emphasis_state_array                             _post_emphasis_states;
  compander_state_array                            _expander_states;
  double_x2                                        _sat_prev;
  std::array<double_x2, n_crossv>                  _crossv_prev;
  uint                                             _n_processed_samples;
  uint                                             _control_rate_mask;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
