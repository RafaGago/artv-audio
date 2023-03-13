#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/ducker.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/own/parts/waveshapers/sigmoid.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

#define LOFIVERB_ADD_DEBUG_ALGO   1
#define LOFIVERB_GAIN_CALIBRATION 1

namespace artv {
namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
static constexpr uint max_block_size = 32;

#ifdef LOFIVERB_ADD_DEBUG_ALGO
//------------------------------------------------------------------------------
static constexpr auto get_debug_algo_spec()
{
  return make_array<stage_data> (
    make_ap (147, -0.707), // 0
    make_ap (183, 0.707), // 1
    make_ap (389, -0.6), // 2
    make_ap (401, 0.6), // 3
    make_quantizer(), // 4
    make_delay (max_block_size + 1) // 5 delay block (== blocksz + 1)
  );
}

struct debug_algo_spec {
  static constexpr auto values {get_debug_algo_spec()};
};
#endif
//------------------------------------------------------------------------------
static constexpr auto get_algo1_spec()
{
  return make_array<stage_data> (
    // diffusors
    make_ap (147, -0.707), // 0
    make_ap (183, 0.707), // 1
    make_ap (389, -0.6), // 2
    make_ap (401, 0.6), // 3
    // er (not ER at the end, as this is now another reverb...)
    make_ap (1367, 0.35, 71 + 70), // 4
    make_lp(), // 5
    make_quantizer(), // 6
    make_ap (1787, 0.5, 261), // 7
    make_delay (max_block_size + 1), // 8 to allow block processing
    // loop1
    make_ap (977, 0.5 /*overridden*/, 51), // 9
    make_delay (2819), // 10
    make_lp(), // 11
    make_quantizer(), // 12
    make_ap (863, -0.5 /*overridden*/), // 13
    make_delay (1021), // 14
    make_quantizer(), // 15
    make_ap (1453, 0.618), // 16
    make_delay (787), // 17 delay (allows block processing) (> blocksz + 1)
    // loop2
    make_ap (947, 0.5 /*overridden*/, 67), // 18
    make_delay (3191), // 19
    make_lp(), // 20
    make_quantizer(), // 21
    make_ap (887, -0.5 /*overridden*/), // 22
    make_delay (1049), // 23
    make_quantizer(), // 24
    make_ap (1367, 0.618), // 25
    make_hp (0.98), // 26
    make_delay (647)); // 27 delay (allows block processing) (> blocksz + 1)
}

struct algo1_spec {
  static constexpr auto values {get_algo1_spec()};
};

//------------------------------------------------------------------------------
static constexpr auto get_midifex49_spec()
{
  return make_array<stage_data> (
    // diffusors
    make_ap (321, 0.5), // 0 PreAP
    make_ap (431, 0.5), // 1 PreAP
    make_ap (968, 0.5), // 2 PreAP
    make_ap (1620, 0.5), // 3 PreAP

    make_delay (21), // 4 L1
    make_delay (1010), // 5 R1

    make_delay (1624), // 6 Loop
    make_ap (1992, 0.5, 17), // 7 Loop

    make_delay (1891), // 8 L2
    make_delay (890), // 9 R2

    make_delay (2110), // 10 Loop
    make_quantizer(), // 11
    make_ap (2371, 0.5), // 12 Loop nested allpass 1
    make_ap (1378, 0.2), // 13 Loop nested allpass 2

    make_delay (2003), // 14 L3
    make_delay (671), // 15 R3

    make_delay (2157), // 16 Loop
    make_quantizer(), // 17
    make_lp(), // 18
    make_ap (2712, 0.5, 22), // 19 Loop nested allpass 1
    make_ap (1783, 0.2), // 20 Loop nested allpass 2

    make_hp (0.995), // 21 HP

    make_delay (max_block_size + 1) // 21 delay block (== blocksz + 1)
  );
}

struct midifex49_spec {
  static constexpr auto values {get_midifex49_spec()};
};

//------------------------------------------------------------------------------
static constexpr auto get_midifex50_spec()
{
  return make_array<stage_data> (
    make_quantizer(), // 0
    make_ap (13, 0.5), // 1
    make_ap (83, 0.5), // 2
    make_ap (116, 0.5), // 3
    make_ap (239, 0.5), // 4
    make_delay (2, 32), // 5 TODO: FIX or STUDY. min of 2 spls are required
    make_ap (339, 0.5), // 6
    make_ap (481, 0.5), // 7
    make_ap (555, 0.5), // 8
    make_ap (823, 0.5), // 9
    make_delay (2, 64), // 10 TODO: FIX or STUDY. min of 2 spls are required
    make_ap (999, 0.5), // 11
    make_ap (1100, 0.5), // 12
    make_ap (1347, 0.5), // 13
    make_ap (1563, 0.5), // 14
    make_delay (2, 64), // 15 TODO: FIX or STUDY. min of 2 spls are required
    make_ap (1841, 0.5), // 16
    make_ap (2001, 0.5, 67), // 17
    make_ap (2083, 0.5, 127), // 18
    make_lp(), // 19
    make_delay (2, 96), // 20 TODO: FIX or STUDY. min of 2 spls are required
    make_hp (0.995), // 21 HP
    make_delay (max_block_size + 1), // 22 (FB point) (== blocksz + 1)
    make_ap (147, 0.5), // 23 L diff
    make_ap (43, 0.5), // 24 L diff
    make_ap (55, 0.5), // 25 L diff
    make_ap (249, 0.5), // 26 R diff
    make_ap (48, 0.5), // 27 R diff
    make_ap (21, 0.5) // 28 R diff
  );
}

struct midifex50_spec {
  static constexpr auto values {get_midifex50_spec()};
};

}} // namespace detail::lofiverb
//------------------------------------------------------------------------------
// A reverb using 16-bit fixed-point arithmetic on the main loop. One design
// criteria has been for it to be extremely CPU friendly.
// Notice: __fp16 (half-precision floating point storage) could achieve the
// same memory savings with more dynamic range. This is still kept as a 16-bit
// fixed-point reverb.
class lofiverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    if (v == _param.mode) {
      return;
    }
    switch (v) {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo_flt:
    case mode::debug_algo: {
      _param.pre_gain  = db_to_gain (-29.f);
      _param.post_gain = db_to_gain (29.f - 6.f); // match all algos
      auto& rev        = _modes.emplace<debug_algo_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    } break;
#endif
    case mode::algo1_flt:
    case mode::algo1: {
      _param.pre_gain  = db_to_gain (-29.f);
      _param.post_gain = db_to_gain (29.f - 6.f); // match all algos
      auto& rev        = _modes.emplace<algo1_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    } break;
    case mode::midifex49_flt:
    case mode::midifex49: {
      _param.pre_gain  = db_to_gain (-24.f);
      _param.post_gain = db_to_gain (24.f - 10.6f);
      auto& rev        = _modes.emplace<midifex49_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    case mode::midifex50_flt:
    case mode::midifex50: {
      _param.pre_gain  = db_to_gain (-23.f);
      _param.post_gain = db_to_gain (23.f - 6.f);
      auto& rev        = _modes.emplace<midifex50_type>();
      rev.reset_memory (_mem_reverb);
      _lfo.set_phase (
        phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    } break;
    default:
      return;
    }
    _param.mode          = v;
    _n_processed_samples = 0; // trigger the control block on first sample
    xspan_memset (_mem_reverb, 0);
    update_mod();
  }
  struct mode {
    enum {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
      debug_algo_flt,
      debug_algo,
#endif
      algo1_flt,
      algo1,
      midifex49_flt,
      midifex49,
      midifex50_flt,
      midifex50
    };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
#ifdef LOFIVERB_ADD_DEBUG_ALGO
        "Debug ",
        "Debug 16-bit",
#endif
        "Long 1",
        "Long 1 16-bit",
        "Midifex 49",
        "Midifex 49 16-bit",
        "Midifex 50",
        "Midifex 50 16-bit"),
      128);
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};
  void set (decay_tag, float v)
  {
    _param_smooth.target().decay = v * 0.01f * fixpt_max_flt;
  }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mod_tag {};
  void set (mod_tag, float v)
  {
    v *= 0.01f * fixpt_max_flt;
    if (v == _param_smooth.target().mod) {
      return;
    }
    _param_smooth.target().mod = v;
    update_mod();
  }

  static constexpr auto get_parameter (mod_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct character_tag {};
  void set (character_tag, float v)
  {
    _param_smooth.target().character = v * 0.01f * fixpt_max_flt;
  }

  static constexpr auto get_parameter (character_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct damp_tag {};
  void set (damp_tag, float v) { _param.damp = v * 0.01f * fixpt_max_flt; }

  static constexpr auto get_parameter (damp_tag)
  {
    return float_param ("%", 0., 100., 30., 0.001);
  }
  //----------------------------------------------------------------------------
  struct freq_balace_tag {};
  void set (freq_balace_tag, float v)
  {
    v *= 0.01f;
    if (v == _param.tilt) {
      return;
    }
    _param.tilt = v;
    auto db     = vec_set<2> ((float) v * -14.f);
    _filt.reset_coeffs<0> (vec_set<2> (350.f), vec_set<2> (0.5f), db, t_spl);
    _filt.reset_coeffs<1> (vec_set<2> (2200.f), db * -0.5f);
  }

  static constexpr auto get_parameter (freq_balace_tag)
  {
    return float_param ("%", -100, 100., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_predelay_qb = 4;
  struct predelay_tag {};
  void set (predelay_tag, float v) { _param.predelay = v; }

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("quarters", 0., max_predelay_qb, 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct clip_level_tag {};
  void set (clip_level_tag, float v) { _param.gain = db_to_gain (v); }

  static constexpr auto get_parameter (clip_level_tag)
  {
    return float_param ("dBFS", -36., 36., 0., 0.1, 0.5, true);
  }
  //----------------------------------------------------------------------------
  struct stereo_tag {};
  void set (stereo_tag, float v) { _param_smooth.target().stereo = v * 0.01f; }

  static constexpr auto get_parameter (stereo_tag)
  {
    return float_param ("%", -100., 100., 100., 0.01);
  }
  //----------------------------------------------------------------------------
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v)
  {
    if (v == _param.ducking_threshold) {
      return;
    }
    _param.ducking_threshold = v;
    _ducker.set_threshold (vec_set<2> (v));
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -60.f, 0.0f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    if (v == _param.ducking_speed) {
      return;
    }
    _param.ducking_speed = v;
    v *= 0.01f;
    _ducker.set_speed (vec_set<2> (v * v), t_spl);
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    auto beat_hz         = pc.get_play_state().bpm * (1.f / 60.f);
    _1_4beat_spls        = (0.25f / beat_hz) * srate;
    _n_processed_samples = 0;

    _resampler.reset (
      srate,
      pc.get_sample_rate(),
      10500,
      10500,
      32,
      16,
      210,
      true,
      max_block_size,
      6 * 1024);

    _filt.reset_states_cascade();
    _ducker.reset();
    _lfo.reset();

    // resize memory
    uint predelay_spls = std::ceil (_1_4beat_spls * max_predelay_qb) * 2;
    uint rev_spls      = 0;
    mp11::mp_for_each<decltype (_modes)> ([&] (auto mode) {
      rev_spls = std::max (mode.get_required_size(), rev_spls);
    });
    uint predelay_flt = predelay_spls * (sizeof (float) / sizeof _mem[0]);

    _mem.clear();
    _mem.resize (predelay_spls + predelay_flt + rev_spls);

    _mem_reverb = xspan {_mem};
    // reminder, hack, _mem_reverb is span<s16> and here it is casting to
    // float. Do keep the predelay (float) first in the memory block to avoid
    // alignment issues.
    _predelay.reset (_mem_reverb.cut_head (predelay_flt).cast<float>(), 2);

    // set defaults
    mp11::mp_for_each<parameters> ([&] (auto param) {
      set (param, get_parameter (param).min);
      if constexpr (!is_choice<decltype (get_parameter (param))>) {
        set (param, get_parameter (param).max);
      }
      else {
        // max might be not yet impl.
        set (param, get_parameter (param).min + 1);
      }
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
#if MOD_DBG_DENORMALS
    feenableexcept (FE_INVALID);
#endif
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    _resampler.process (outs, ins, samples, [=] (auto io) {
      if ((_param.mode % 2) == 0) {
        process_float (io);
      }
      else {
        process_fixed_point (io);
      }
    });
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    mode_tag,
    character_tag,
    damp_tag,
    decay_tag,
    predelay_tag,
    freq_balace_tag,
    clip_level_tag,
    mod_tag,
    stereo_tag,
    ducking_threshold_tag,
    ducking_speed_tag>;
  //----------------------------------------------------------------------------
private:
  // fixpt_t truncates when dropping fractional bits (leaky).
  // fixpt_tr rounds to the nearest (never reaches full zero).
  //
  // Using fixpt_t as default, with fixpt_tr at some points compensate the
  // truncating leakage.
  using fixpt_t  = detail::lofiverb::fixpt_t;
  using fixpt_tr = detail::lofiverb::fixpt_tr;
  //----------------------------------------------------------------------------
  static constexpr uint  max_block_size = detail::lofiverb::max_block_size;
  static constexpr uint  n_channels     = 2;
  static constexpr uint  srate          = 23400;
  static constexpr float t_spl          = (float) (1. / srate);
  // this is using 16 bits fixed-point arithmetic, positive values can't
  // represent one, so instead of correcting everywhere the parameters are
  // scaled instead to never reach 1.
  static constexpr auto fixpt_max_flt = fixpt_t::max_float();
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr auto one()
  {
    if constexpr (std::is_floating_point_v<T>) {
      return fixpt_max_flt;
    }
    else {
      return fixpt_t::max();
    }
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  void process_fixed_point (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    array2d<fixpt_t, 2, max_block_size> wet;
    std::array<f32_x2, max_block_size>  ducker_gain;
    loop_parameters                     pars;

    auto pre_gain = (1.f / _param.gain) * fixpt_max_flt;

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wetv = _filt.tick_cascade (vec_from_array (io[i]));
#ifndef LOFIVERB_GAIN_CALIBRATION
      wetv *= pre_gain;
      wetv           = vec_clamp (wetv, -fixpt_max_flt, fixpt_max_flt);
      ducker_gain[i] = _ducker.tick (wetv);
      wetv *= _param.pre_gain;
#else // pre_gain calibration. (shouldn't be on when Released).
      wetv *= pre_gain;
      wetv = vec_clamp (wetv, -fixpt_max_flt, fixpt_max_flt);
#endif
      io[i] = vec_to_array (wetv);

      _param_smooth.tick();
      pars.stereo[i] = _param_smooth.get().stereo;
      pars.decay[i].load_float (_param_smooth.get().decay);
      pars.character[i].load_float (_param_smooth.get().character);
      pars.mod[i].load_float (_param_smooth.get().mod);
    }
    // predelay + int conversion
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0].load_float (_predelay.get (predelay_spls, 0));
        wet[i][1].load_float (_predelay.get (predelay_spls, 1));
        _predelay.push (io[i]);
      }
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0].load_float (io[i][0]);
        wet[i][1].load_float (io[i][1]);
      }
    }
    // main loop
    switch (_param.mode) {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo:
      process_debug_algo (xspan {wet.data(), io.size()}, pars);
      break;
#endif
    case mode::algo1:
      process_algo1 (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::midifex49:
      process_midifex49 (xspan {wet.data(), io.size()}, pars);
      break;
    case mode::midifex50:
      process_midifex50 (xspan {wet.data(), io.size()}, pars);
      break;
    default:
      assert (false);
    }
    // float conversion
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0].to_floatp();
      auto r = wet[i][1].to_floatp();
#ifndef LOFIVERB_GAIN_CALIBRATION
      l *= ducker_gain[i][0] * _param.gain * _param.post_gain;
      r *= ducker_gain[i][1] * _param.gain * _param.post_gain;
#else // pre_gain calibration (shouldn't be on when Released).
      l *= _param.gain;
      r *= _param.gain;
#endif

      l = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
  //----------------------------------------------------------------------------
  void process_float (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    std::array<f32_x2, max_block_size> ducker_gain;
    loop_parameters_flt                pars;

    auto pre_gain = 1.f / _param.gain;

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wet = _filt.tick_cascade (vec_from_array (io[i]));
      wet *= pre_gain;
      wet            = vec_clamp (wet, -1.f, 1.f);
      ducker_gain[i] = _ducker.tick (wet);
      wet *= _param.pre_gain;
      io[i] = vec_to_array (wet);

      _param_smooth.tick();
      pars.stereo[i]    = _param_smooth.get().stereo;
      pars.decay[i]     = _param_smooth.get().decay;
      pars.character[i] = _param_smooth.get().character;
      pars.mod[i]       = _param_smooth.get().mod;
    }
    // predelay
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        // TODO: block fetch?
        std::array<float, 2> spl;
        spl[0] = _predelay.get (predelay_spls, 0);
        spl[1] = _predelay.get (predelay_spls, 1);
        _predelay.push (xspan {io[i]});
        io[i] = spl;
      }
    }
    // main loop
    switch (_param.mode) {
#ifdef LOFIVERB_ADD_DEBUG_ALGO
    case mode::debug_algo_flt:
      process_debug_algo (io, pars);
      break;
#endif
    case mode::algo1_flt:
      process_algo1 (io, pars);
      break;
    case mode::midifex49_flt:
      process_midifex49 (io, pars);
      break;
    case mode::midifex50_flt:
      process_midifex50 (io, pars);
      break;
    default:
      assert (false);
    }
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = io[i][0] * ducker_gain[i][0] * _param.gain * _param.post_gain;
      auto r = io[i][1] * ducker_gain[i][1] * _param.gain * _param.post_gain;
      l      = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
//----------------------------------------------------------------------------
#ifdef LOFIVERB_ADD_DEBUG_ALGO
  template <class T, class Params>
  void process_debug_algo (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<debug_algo_type> (_modes);

    using arr    = std::array<T, max_block_size>;
    using arr_fb = std::array<T, max_block_size + 1>;

    arr    in_arr;
    arr_fb out_arr;

    auto in  = xspan {in_arr.data(), io.size()};
    auto out = xspan {out_arr.data(), io.size() + 1};

    rev.fetch_block<5> (out, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      in[i] = (T) (((io[i][0] + io[i][1]) * 0.25_r) + out[i]);
    }
    out.cut_head (1); // drop feedback sample from previous block

    rev.run<0> (in);
    rev.run<1> (in);
    rev.run<2> (in);
    rev.run<3> (in);
    rev.run<4> (in, [&] (auto v, uint i) { // decay fixup
      auto decay = (T) (one<T>() - par.decay[i]);
      decay      = (T) (one<T>() - decay * decay);
      decay      = (T) (0.6_r + decay * 0.39_r);
      return v * decay;
    });
    rev.push<5> (in.to_const()); // feedback point

    // Mixdown (not a lot going on...)
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      io[i][0] = out[i];
      io[i][1] = out[i];
    }
  }
#endif
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_algo1 (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<algo1_type> (_modes);

    using arr    = std::array<T, max_block_size>;
    using arr_fb = std::array<T, max_block_size + 1>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (T) ((io[i][0] + io[i][1]) * 0.5_r);
      auto mod   = (T) (0.25_r + (1_r - par.mod[i]) * 0.75_r);
      // ER + late lfo
      auto lfo = tick_lfo<T>();
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = (T) (T {lfo[1]} * (0.5_r + par.character[i] * 0.5_r));
      lfo3[i]  = (T) (T {lfo[2]} * mod);
      lfo4[i]  = (T) (T {lfo[3]} * mod);

      // decay fixup
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.6_r + decay * 0.38_r);
    }
    // diffusion -----------------------------
    rev.run<0> (late_in);
    rev.run<1> (late_in);
    rev.run<2> (late_in);
    rev.run<3> (late_in);

    // ER (first reverb, not exactly ER at the end...) -----------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    rev.fetch_block<8> (er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = (T) (late_in[i] * 0.5_r + er2[i]);
    });
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run<4> (er1, xspan {lfo2}, nullptr);
    rev.run<5> (er1, load_float<T> (0.0001f + 0.17f * _param.damp));
    xspan_memcpy (er1b, er1);
    rev.run<6> (er1b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<7> (er1b, xspan {lfo1}, nullptr);
    rev.push<8> (er1b.to_const()); // feedback point

    // Late -----------------------------
    arr    late_arr;
    arr_fb l_arr;
    arr_fb r_arr;
    arr    g_arr;

    auto late = xspan {late_arr.data(), io.size()};
    auto g    = xspan {g_arr.data(), io.size()};
    auto l    = xspan {l_arr.data(), io.size() + 1}; // +1: Feedback on head
    auto r    = xspan {r_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    rev.fetch_block<17> (l, 1); // feedback, fetching block + 1 samples
    rev.fetch_block<27> (r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
      g[i] = (T) (0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r));
    }
    r.cut_head (1); // drop feedback sample from previous block

    float late_dampf = (0.9f - _param.damp * 0.9f);
    late_dampf       = 1.f - late_dampf * late_dampf;
    late_dampf *= 0.4f;
    auto late_damp = load_float<T> (late_dampf);
    rev.run<9> (late, xspan {lfo3}, g);
    rev.run<10> (late);
    rev.run<11> (late, late_damp);
    rev.run<12> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<13> (late, nullptr, [g] (uint i) { return -g[i]; });
    rev.run<14> (late);
    rev.run<15> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<16> (late);
    rev.push<17> (late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
    }
    l.cut_head (1); // drop feedback sample from previous block
    rev.run<18> (late, xspan {lfo4}, g);
    rev.run<19> (late);
    rev.run<20> (late, late_damp);
    rev.run<21> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<22> (late, nullptr, [g] (uint i) { return -g[i]; });
    rev.run<23> (late);
    rev.run<24> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<25> (late);
    rev.run<26> (late);
    rev.push<27> (late.to_const()); // feedback point

    // Mixdown
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto er_amt     = par.character[i] * 0.2_r;
      auto e_l        = (-er1[i] * 0.66_r - er2[i] * 0.34_r) * er_amt;
      auto e_r        = (-er1[i] * 0.66_r + er2[i] * 0.34_r) * er_amt;
      auto direct_amt = one<T>() - er_amt;
      io[i][0]        = (T) (l[i] * direct_amt + e_l);
      io[i][1]        = (T) (r[i] * direct_amt + e_r);
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_midifex49 (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<midifex49_type> (_modes);

    using arr = std::array<T, max_block_size>;

    arr   tmp_arr;
    xspan tmp {tmp_arr.data(), io.size()};
    arr   sig_arr;
    xspan sig {sig_arr.data(), io.size()};
    arr   l_arr;
    xspan l {l_arr.data(), io.size()};
    arr   r_arr;
    xspan r {r_arr.data(), io.size()};
    arr   lfo1_arr;
    xspan lfo1 {lfo1_arr.data(), io.size()};
    arr   lfo2_arr;
    xspan lfo2 {lfo2_arr.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      sig[i] = (T) (((spl[0] + spl[1]) * 0.25_r)); // gain = 0.5 + feedback = 1
      // LFO
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      // decay fixup
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.1_r + decay * 0.8375_r);
    });

    rev.run<0> (sig);
    rev.run<1> (sig);
    rev.run<2> (sig);
    rev.run<3> (sig);

    // feedback handling, fetching the block with a negative offset of 1
    rev.fetch_block<22> (tmp, 1);
    span_visit (sig, [&] (auto& v, uint i) { v = (T) (v + tmp[i]); });

    // 1st output point for L and R signal
    xspan_memcpy (l, sig); // delay LT a block -> might overlap, requires copy
    rev.run<4> (l);
    rev.fetch_block<5> (r); // delay GT a block will never overlap
    rev.push<5> (sig.to_const());

    // continuing the loop
    rev.run<6> (sig);
    rev.run<7> (sig, lfo1, nullptr);

    // 2nd output point for L and R signal
    rev.fetch_block<8> (tmp); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<8> (sig.to_const());
    rev.fetch_block<9> (tmp); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<9> (sig.to_const());

    // continuing the loop
    rev.run<10> (sig);
    rev.run<11> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.14_r);
    });
    rev.run<12, 13> (sig, nullptr, nullptr, nullptr, tmp);

    // 3rd output point for L and R signal
    rev.fetch_block<14> (tmp); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<14> (sig.to_const());
    rev.fetch_block<15> (tmp); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<15> (sig.to_const());

    // continuing the loop
    rev.run<16> (sig);
    rev.run<17> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<18> (sig, load_float<T> (_param.damp * 0.8f));
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.2_r);
    });
    rev.run<19, 20> (sig, lfo2, nullptr, nullptr, tmp);
    rev.run<21> (sig);
    // push to delay feedback
    rev.push<22> (sig.to_const());

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class Params>
  void process_midifex50 (xspan<std::array<T, 2>> io, Params& par)
  {
    auto& rev = std::get<midifex50_type> (_modes);

    using arr = std::array<T, max_block_size>;

    arr   l_arr;
    xspan l {l_arr.data(), io.size()};
    arr   r_arr;
    xspan r {r_arr.data(), io.size()};
    arr   lfo1_arr;
    xspan lfo1 {lfo1_arr.data(), io.size()};
    arr   lfo2_arr;
    xspan lfo2 {lfo2_arr.data(), io.size()};

    // fetch feedback values
    rev.fetch_block<22> (l, 1);
    rev.run<0> (l, [&] (T fb_spl, uint i) {
      // running lfo (unrelated to quantizing)
      auto lfo = tick_lfo<T>();
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);

      // decay fixup
      auto decay = (T) (one<T>() - par.decay[i]);
      decay      = (T) (one<T>() - decay * decay);
      decay      = (T) - (0.3_r + decay * 0.45_r);

      return (fb_spl * decay) + ((io[i][0] + io[i][1]) * 0.25_r); // gain = 1
    }); // feedback + input summing with quantizer

    rev.run<1> (l);
    rev.run<2> (l);
    rev.run<3> (l);
    rev.run<4> (l);
    rev.run<5> (l, xspan {par.character}, nullptr); // variable delay

    rev.run<6> (l);
    rev.run<7> (l);
    rev.run<8> (l);
    rev.run<9> (l);
    rev.run<10> (l, xspan {par.character}, nullptr); // variable delay

    rev.run<11> (l);
    rev.run<12> (l);
    rev.run<13> (l);
    rev.run<14> (l);
    rev.run<15> (l, xspan {par.character}, nullptr); // variable delay

    rev.run<16> (l);
    rev.run<17> (l, xspan {lfo1}, nullptr);
    rev.run<18> (l, xspan {lfo2}, nullptr);
    float dampv = _param.damp;
    dampv       = dampv * dampv * 0.4f;
    dampv       = 0.05f + dampv;
    rev.run<19> (l, load_float<T> (dampv));
    rev.run<20> (l, xspan {par.character}, nullptr); // variable delay
    rev.run<21> (l);

    // push to delay feedback
    rev.push<22> (l.to_const());

    xspan_memcpy (r, l);
    rev.run<23> (l);
    rev.run<24> (l);
    rev.run<25> (l);

    rev.run<26> (r);
    rev.run<27> (r);
    rev.run<28> (r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T load_float (float v)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      return fixpt_t::from_float (v);
    }
    else {
      return v;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  auto tick_lfo()
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      auto ret = _lfo.tick_sine_fixpt().spec_cast<fixpt_t>().value();
      // At this point "ret" has fixed traits, but different fixed point
      // conversions configured
      return vec_cast<fixpt_t::value_type> (ret);
    }
    else {
      return _lfo.tick_sine();
    }
  }
  //----------------------------------------------------------------------------
  void update_mod()
  {
    // reminder, the phase relations are set on "void set (mode_tag, int v)"
    auto mod = _param_smooth.target().mod;
    switch (_param.mode) {
    case mode::algo1:
    case mode::algo1_flt: {
      auto f_er   = 0.3f + mod * 0.3f;
      auto f_late = 0.1f + mod * 1.2f;
      _lfo.set_freq (f32_x4 {f_er, f_er, f_late, f_late}, t_spl);
    } break;
    case mode::midifex49_flt:
    case mode::midifex49: {
      auto f_late = 0.2f + mod * 0.2f;
      _lfo.set_freq (f32_x4 {f_late, f_late, f_late, f_late}, t_spl);
    } break;
    case mode::midifex50_flt:
    case mode::midifex50: {
      auto f_late = 0.3f + mod * 0.1f;
      _lfo.set_freq (f32_x4 {f_late, f_late, f_late, f_late}, t_spl);
    } break;
    default:
      break;
    }
  }
  //----------------------------------------------------------------------------
  // just a convenience function for iterating block loops while not bloating
  // the code with more loop unroll hints than necessary
  template <class T, class F>
  void span_visit (xspan<T> block, F&& visitor)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block.size(); ++i) {
      visitor (block[i], i);
    }
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   mode;
    float post_gain; // set to rougly match loudness for all algorithms and to
                     // compensate for pre_gain
    float pre_gain; // using a 440Hz 0dbFS square wave in looking at which
                    // level breaks with decay set to the maximum value and
                    // moving abruptly all parameters.
    float gain; // a parameter
    float tilt;
    float predelay;
    float damp;
    float ducking_threshold;
    float ducking_speed;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float mod;
    float stereo;
    float er;
    float decay;
    float character;
    // dry, wet, ducker/gate
  };
  //----------------------------------------------------------------------------
  struct loop_parameters {
    std::array<fixpt_t, max_block_size> decay;
    std::array<fixpt_t, max_block_size> character;
    std::array<fixpt_t, max_block_size> mod;
    std::array<float, max_block_size>   stereo;
  };
  struct loop_parameters_flt {
    std::array<float, max_block_size> decay;
    std::array<float, max_block_size> character;
    std::array<float, max_block_size> mod;
    std::array<float, max_block_size> stereo;
  };
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>                                       _resampler {};
  part_classes<mp_list<tilt_eq, onepole_naive_highshelf>, f32_x2> _filt {};
  static_delay_line<float, true, false>                           _predelay;

  lfo<4> _lfo;

  using algo1_type
    = detail::lofiverb::engine<detail::lofiverb::algo1_spec, max_block_size>;
  using midifex49_type = detail::lofiverb::
    engine<detail::lofiverb::midifex49_spec, max_block_size>;
  using midifex50_type = detail::lofiverb::
    engine<detail::lofiverb::midifex50_spec, max_block_size>;

#ifndef LOFIVERB_ADD_DEBUG_ALGO
  using modes_type = std::variant<algo1_type, midifex49_type, midifex50_type>;
#else
  using debug_algo_type = detail::lofiverb::
    engine<detail::lofiverb::debug_algo_spec, max_block_size>;

  using modes_type
    = std::variant<debug_algo_type, algo1_type, midifex49_type, midifex50_type>;
#endif

  modes_type     _modes;
  ducker<f32_x2> _ducker;

  uint  _n_processed_samples;
  float _1_4beat_spls;

  std::vector<s16, overaligned_allocator<s16, 16>> _mem;
  xspan<s16>                                       _mem_reverb;
};
//------------------------------------------------------------------------------
} // namespace artv
