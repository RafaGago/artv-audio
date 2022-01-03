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
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
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
class waveshaper {
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
      1, make_cstr_array ("Hardclip", "Tanh", "Sqrt", "SqrtSin"), 40);
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
      case mode_wvsh_1b:
      case mode_wvsh_1a:
      case mode_wvsh_pow:
      case mode_wvsh_sqrt:
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
        "Normal no AA",
        "Normal",
        "Low Levels",
        "High Levels",
        "Very Low Levels",
        "Very High Levels"),
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
    _p.crossv_hz[lo_crossv] = v;
    if (_p.crossv_enabled[lo_crossv]) {
      update_crossover (lo_crossv, false);
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
    if (v == _p.crossv_hz[hi_crossv]) {
      return;
    }
    _p.crossv_hz[hi_crossv] = v;
    if (_p.crossv_enabled[hi_crossv]) {
      update_crossover (hi_crossv, false);
    }
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
    _p.crossv_enabled[lo_crossv] = (v != 0);
    int  order                   = v;
    bool is_lp                   = order <= max_crossv_order;
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
        "Off",
        "HP 6dB/Oct",
        "HP 12dB/Oct",
        "HP 18dB/Oct",
        "HP 24dB/Oct",
        "HP 30dB/Oct",
        "HP 36dB/Oct",
        "HP 42dB/Oct",
        "HP 48dB/Oct",
        "LR 6dB/Oct",
        "LR 12dB/Oct",
        "LR 18dB/Oct",
        "LR 24dB/Oct",
        "LR 30dB/Oct",
        "LR 36dB/Oct",
        "LR 42dB/Oct",
        "LR 48dB/Oct"),
      32);
  }
  //----------------------------------------------------------------------------
  struct hi_mode_tag {};

  void set (hi_mode_tag, int v)
  {
    _p.crossv_enabled[hi_crossv] = (v != 0);
    int  order                   = v;
    bool is_lp                   = (order > max_crossv_order);
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
    return choice_param (
      0,
      make_cstr_array (
        "Off",
        "LP 6dB/Oct",
        "LP 12dB/Oct",
        "LP 18dB/Oct",
        "LP 24dB/Oct",
        "LP 30dB/Oct",
        "LP 36dB/Oct",
        "LP 42dB/Oct",
        "LP 48dB/Oct",
        "HR 6dB/Oct",
        "HR 12dB/Oct",
        "HR 18dB/Oct",
        "HR 24dB/Oct",
        "HR 30dB/Oct",
        "HR 36dB/Oct",
        "HR 42dB/Oct",
        "HR 48dB/Oct"),
      32);
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
    return float_param ("ms", 1, 180., 20., 0.02);
  }
  //----------------------------------------------------------------------------
  struct envfollow_mode_tag {};

  void set (envfollow_mode_tag, int v) { _p.ef_mode = v; }

  static constexpr auto get_parameter (envfollow_mode_tag)
  {
    return choice_param (
      2, make_cstr_array ("sqrt(exp)", "exp", "exp^2", "exp^3"), 20);
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
    return float_param ("ms", 1., 300., 56., 0.02);
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

    _envfollow.reset_states();
    _pre_emphasis.reset_states();
    _post_emphasis.reset_states();
    _crossv.zero_all_states();
    _wvsh.zero_all_states();

    // all instances have the same coefficients
    _wvsh.reset_coeffs_on_idx<sqrt_aa> (compand_in);
    _wvsh.reset_coeffs_on_idx<sqrt_aa> (sigmoid_wsh);
    _wvsh.reset_coeffs_on_idx<sqrt_aa> (compand_out);

    // this DC blocker at a very low frequency is critical for the sound of
    // the companded modes to be acceptable. Unfortunately it causes DC itself
    // when the sound is muted. I didn't find a solution.
    _dc_block.reset_coeffs (make_vec_x1 (1.), pc.get_sample_rate());

    _p = params {};

    _sparams.reset ((float) pc.get_sample_rate());

    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
    _p.type_prev = sat_type_count; // force initial click removal routine.

    uint sr_order        = get_samplerate_order (pc.get_sample_rate()) + 3;
    _control_rate_mask   = lsb_mask<uint> (sr_order);
    _n_processed_samples = 0;
    _sat_prev = _dcmod_prev = _crossv_prev[lo_crossv] = _crossv_prev[hi_crossv]
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
      _dc_block.reset_states();
      _wvsh.zero_all_states();
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

      // interleaving
      for (uint i = 0; i < subblock_size; ++i) {
        sat[i][0] = ins[0][samples_processed + i];
        sat[i][1] = ins[1][samples_processed + i];
      }
      // Input DC blocking
      for (uint i = 0; i < subblock_size; ++i) {
        sat[i] = _dc_block.tick_on_idx (dc_block_in, sat[i]);
      }
      // pre process/Companding + data initialization block
      //
      // Using "compand_1b_aa" first required extremely good DC blocking
      // unfortunately. The sound easily broke.
      switch (p.mode) {
      case mode_wvsh_1b:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<compand_1b_aa> (compand_in, sat[i]);
        }
        break;
      case mode_wvsh_1a:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<compand_1a_aa> (compand_in, sat[i]);
        }
        break;
      case mode_wvsh_pow:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<pow2_aa> (compand_in, sat[i]);
        }
        break;
      case mode_wvsh_sqrt:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<sqrt_aa> (compand_in, sat[i]);
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
          if (p.crossv_is_lp[lo_crossv]) {
            lo[i] = _crossv.tick_on_idx<lp_type> (
              lo_crossv, sat[i], p.crossv_order[lo_crossv]);
            sat[i] = sat[i] - lo[i];
          }
          else {
            lo[i] = _crossv.tick_on_idx<hp_type> (
              lo_crossv, sat[i], p.crossv_order[lo_crossv]);
            sat[i] = lo[i];
            lo[i]  = double_x2 {};
          }
        }
        else {
          lo[i] = double_x2 {};
        }

        if (_p.crossv_enabled[hi_crossv]) {
          if (!p.crossv_is_lp[hi_crossv]) {
            hi[i] = _crossv.tick_on_idx<hp_type> (
              hi_crossv, sat[i], p.crossv_order[hi_crossv]);
            sat[i] = sat[i] - hi[i];
          }
          else {
            hi[i] = _crossv.tick_on_idx<lp_type> (
              hi_crossv, sat[i], p.crossv_order[hi_crossv]);
            sat[i] = hi[i];
            hi[i]  = double_x2 {};
          }
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
        double_x2 follow = _envfollow.tick (sat[i]);
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
        auto ef2dc = _sparams.get (sm_ef_to_dc);
        auto amt   = fabs (ef2dc) * 0.40;
        auto dcmod = follow * amt;
        dcmod += _dcmod_prev;
        dcmod *= 0.5; // Moving average LP
        sat[i] += dcmod;

        sat[i] = _pre_emphasis.tick (sat[i]);
        // waveshaping
        auto wsh_type
          = p.type + (sat_type_count * (uint) waveshaper_type_is_adaa (p.mode));

        switch (wsh_type) {
        case sat_hardclip:
          sat[i] = _wvsh.tick_on_idx<hardclip_adaa<0>> (sigmoid_wsh, sat[i]);
          break;
        case sat_tanh:
          sat[i] = _wvsh.tick_on_idx<tanh_adaa<0>> (sigmoid_wsh, sat[i]);
          break;
        case sat_sqrt:
          sat[i]
            = _wvsh.tick_on_idx<sqrt_sigmoid_adaa<0>> (sigmoid_wsh, sat[i]);
          break;
        case sat_sqrt_sin:
          sat[i]
            = _wvsh.tick_on_idx<sqrt_sin_sigmoid_adaa<0>> (sigmoid_wsh, sat[i]);
          break;
        case sat_hardclip_adaa:
          sat[i] = _wvsh.tick_on_idx<hardclip_aa> (sigmoid_wsh, sat[i]);
          break;
        case sat_tanh_adaa:
          sat[i] = _wvsh.tick_on_idx<tanh_aa> (sigmoid_wsh, sat[i]);
          break;
        case sat_sqrt_adaa:
          sat[i] = _wvsh.tick_on_idx<sqrt_sigmoid_aa> (sigmoid_wsh, sat[i]);
          break;
        case sat_sqrt_sin_adaa:
          sat[i] = _wvsh.tick_on_idx<sqrt_sin_sigmoid_aa> (sigmoid_wsh, sat[i]);
          break;
        default:
          assert (false);
          break;
        }
        sat[i] = _post_emphasis.tick (sat[i]);

        sat[i] -= dcmod;
        _sat_prev = sat[i];

        _dcmod_prev = dcmod;
        // gain and crossover join
        sat[i] *= inv_drive * _sparams.get (sm_out);
        sat[i] += lo[i] + hi[i];
      }
      // post process / restore if required
      switch (p.mode) {
      case mode_wvsh_1b:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<compand_1a_aa> (compand_out, sat[i]);
        }
        break;
      case mode_wvsh_1a:
        for (uint i = 0; i < subblock_size; ++i) {
          sat[i] = _wvsh.tick_on_idx<compand_1b_aa> (compand_out, sat[i]);
        }
        break;
      case mode_wvsh_pow:
        for (uint i = 0; i < subblock_size; ++i) {
          // sat[i] = vec_clamp (sat[i], -pow2compand_clip, pow2compand_clip);
          sat[i] = _wvsh.tick_on_idx<sqrt_aa> (compand_out, sat[i]);
        }
        break;
      case mode_wvsh_sqrt:
        for (uint i = 0; i < subblock_size; ++i) {
          // sat[i] = vec_clamp (sat[i], -pow2compand_clip, pow2compand_clip);
          sat[i] = _wvsh.tick_on_idx<pow2_aa> (compand_out, sat[i]);
        }
        break;
      case mode_no_aa:
      case mode_normal:
      default:
        break;
      }
      // Post DC-blocker.
      for (uint i = 0; i < subblock_size; ++i) {
        sat[i] = _dc_block.tick_on_idx (dc_block_out, sat[i]);
      }
      // Deinterleaving
      for (uint i = 0; i < subblock_size; ++i) {
        outs[0][samples_processed + i] = sat[i][0];
        outs[1][samples_processed + i] = sat[i][1];
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
    auto f = double_x2 {_p.crossv_hz[idx], _p.crossv_hz[idx] * 1.009};
    if (_p.crossv_is_lp[idx]) {
      _crossv.reset_coeffs_on_idx<lp_type> (
        idx, f, _plugcontext->get_sample_rate(), _p.crossv_order[idx]);
    }
    else {
      _crossv.reset_coeffs_on_idx<hp_type> (
        idx, f, _plugcontext->get_sample_rate(), _p.crossv_order[idx]);
    }

    if (reset_states) {
      _crossv_prev[idx] = double_x2 {};
      _crossv.reset_states_on_idx (idx);
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
  static constexpr bool has_wvsher (char mode)
  {
    return mode >= mode_wvsher_first;
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

    _pre_emphasis.reset_coeffs (
      f, q, db, _plugcontext->get_sample_rate(), bell_tag {});

    _post_emphasis.reset_coeffs (
      f, q, -db, _plugcontext->get_sample_rate(), bell_tag {});
  }
  //----------------------------------------------------------------------------
  void update_envelope_follower()
  {
    _envfollow.reset_coeffs (
      vec_set<double_x2> (_p.ef_attack),
      vec_set<double_x2> (_p.ef_release),
      _plugcontext->get_sample_rate());
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

  enum crossv_indexes { lo_crossv, hi_crossv, n_crossv };

  enum mode_type {
    mode_no_aa,
    mode_normal,
    mode_wvsher_first,
    mode_wvsh_1b = mode_wvsher_first,
    mode_wvsh_1a,
    mode_wvsh_pow,
    mode_wvsh_sqrt,
    mode_count,
  };

  enum waveshaper_idx { compand_in, sigmoid_wsh, compand_out, n_waveshapers };

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

  static constexpr uint adaa_order = 1;

  using sqrt_sigmoid_aa = adaa::fix_eq_and_delay<adaa_order, sqrt_sigmoid_adaa>;
  using tanh_aa         = adaa::fix_eq_and_delay<adaa_order, tanh_adaa>;
  using hardclip_aa     = adaa::fix_eq_and_delay<adaa_order, hardclip_adaa>;
  using sqrt_sin_sigmoid_aa
    = adaa::fix_eq_and_delay<adaa_order, sqrt_sin_sigmoid_adaa>;

  using sqrt_aa       = adaa::fix_eq_and_delay<adaa_order, sqrt_adaa>;
  using pow2_aa       = adaa::fix_eq_and_delay<adaa_order, pow2_adaa>;
  using compand_1a_aa = adaa::fix_eq_and_delay<adaa_order, compand_1a_adaa>;
  using compand_1b_aa = adaa::fix_eq_and_delay<adaa_order, compand_1b_adaa>;
  using waveshapers   = parts_union_array<
    mp_list<
      sqrt_sigmoid_adaa<0>,
      sqrt_sigmoid_aa,
      tanh_adaa<0>,
      tanh_aa,
      hardclip_adaa<0>,
      hardclip_aa,
      sqrt_sin_sigmoid_adaa<0>,
      sqrt_sin_sigmoid_aa,
      sqrt_aa,
      pow2_aa,
      compand_1a_aa,
      compand_1b_aa>,
    double_x2,
    n_waveshapers>;

  using lp_type = butterworth_lowpass_any_order;
  using hp_type = butterworth_highpass_any_order;

  using crossv
    = parts_union_array<mp_list<lp_type, hp_type>, double_x2, n_crossv>;

  using dc_blockers = part_class_array_coeffs_global<
    mystran_dc_blocker,
    double_x2,
    dc_block_count>;

  params                                    _p;
  value_smoother<float, sm_count>           _sparams;
  part_class_array<andy::svf, double_x2>    _pre_emphasis;
  crossv                                    _crossv;
  waveshapers                               _wvsh;
  part_class_array<slew_limiter, double_x2> _envfollow;
  dc_blockers                               _dc_block;
  part_class_array<andy::svf, double_x2>    _post_emphasis;
  double_x2                                 _sat_prev;
  double_x2                                 _dcmod_prev;
  std::array<double_x2, n_crossv>           _crossv_prev;
  uint                                      _n_processed_samples;
  uint                                      _control_rate_mask;

  plugin_context* _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
} // namespace artv
