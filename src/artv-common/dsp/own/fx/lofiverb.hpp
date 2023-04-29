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
#include "artv-common/dsp/own/fx/lofiverb-algorithms.hpp"
#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"
#include "artv-common/dsp/own/parts/filters/composite/tilt.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/own/parts/waveshapers/sigmoid.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

// #define LOFIVERB_DEBUG_ALGO 1

namespace artv {

// TODO list
// - predelay broken
// - rename room to ambience
// - dre algorithms are too dark now
// - add pre and post hooks for eq and gain
// - companding on saving fixed point? e.g from a Q0.20
// - wider fixed point type
// - dither?
// - multitaps with lerp, to do modulations on output taps

//------------------------------------------------------------------------------
// A reverb supporting 16-bit fixed-point arithmetic on the main loop. One
// design criteria has been for it to be extremely CPU friendly. Notice: __fp16
// (half-precision floating point storage) could achieve the same memory savings
// with more dynamic range. This is still kept as a 16-bit fixed-point reverb.
class lofiverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct algorithm_tag {};
  void set (algorithm_tag, int v)
  {
    u32 curr = _param.mode / n_bit_formats;
    if (v == curr) {
      return;
    }
    _param.mode = (_param.mode % n_bit_formats) + (v * n_bit_formats);
    update_operation_mode();
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (algorithm_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
#ifndef LOFIVERB_DEBUG_ALGO
        "Artv Abyss",
        "Artv Plate",
        "Artv Room",
        "Acreil Midifex 49",
        "Acreil Midifex 50",
        "Acreil Dre-2000 A",
        "Acreil Dre-2000 B",
        "Acreil Dre-2000 C",
        "Acreil Dre-2000 D",
        "Acreil REV5 L Hall"
#else
        "Debug "
#endif
        ),
      64);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    u32 curr = _param.mode % n_bit_formats;
    if (v == curr) {
      return;
    }
    _param.mode = round_floor (_param.mode, n_bit_formats) + v;
    update_operation_mode();
  }

  static constexpr uint n_bit_formats = 3;

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      1, make_cstr_array ("16-bit fixed", "16-bit float", "32-bit float"), 16);
  }
  //----------------------------------------------------------------------------
  struct clock_tag {};
  void set (clock_tag, int v)
  {
    if (v == _param.srateid) {
      return;
    }
    _param.srateid = v;
    update_operation_mode();
  }

  static constexpr auto get_parameter (clock_tag)
  {
    return choice_param (
      3,
      make_cstr_array (
        "Downclocked 3",
        "Downclocked 2",
        "Downclocked 1",
        "Base Clock",
        "Overclocked 1",
        "Overclocked 2",
        "Overclocked 3"));
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
    return float_param ("%", 0., 100., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lf_amt_tag {};
  void set (lf_amt_tag, float v)
  {
    _param.algo.lf_amt = v * 0.01f * fixpt_max_flt;
  }

  static constexpr auto get_parameter (lf_amt_tag)
  {
    return float_param ("%", 0., 100., 80., 0.01);
  }
  //----------------------------------------------------------------------------
  struct hf_amt_tag {};
  void set (hf_amt_tag, float v)
  {
    _param.algo.hf_amt = v * 0.01f * fixpt_max_flt;
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (hf_amt_tag)
  {
    return float_param ("%", 0., 100., 80., 0.01);
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
  void set (clip_level_tag, float v)
  {
    _param.gain = db_to_gain (v) * clip_value();
  }

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
    update_ducker();
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _pc                  = &pc;
    _n_processed_samples = 0;
    _srate               = 0; // ensure triggering a resample reset
    _param.mode
      = (decltype (_param.mode)) -1ull; // trigger a mode and resampler reset
    _param.srateid = 3;

    _param_smooth.reset (_t_spl);
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
      if ((_param.mode % n_bit_formats) == 0) {
        process_fixed_point (io);
      }
      else {
        process_float (io);
      }
    });
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    algorithm_tag, // keep first! important for reset to work
    mode_tag,
    clock_tag,
    character_tag,
    lf_amt_tag,
    decay_tag,
    predelay_tag,
    hf_amt_tag,
    clip_level_tag,
    mod_tag,
    stereo_tag,
    ducking_threshold_tag,
    ducking_speed_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = detail::lofiverb::max_block_size;
  static constexpr uint n_channels     = 2;
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  void process_fixed_point (xspan<std::array<float, 2>> io)
  {
    using algo = detail::lofiverb::algorithm;
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    array2d<fixpt_t, 2, max_block_size>                       wet;
    std::array<f32_x2, max_block_size>                        ducker_gain;
    detail::lofiverb::algorithm::smoothed_parameters<fixpt_t> pars;

    auto pre_gain = (1.f / _param.gain);

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto wetv      = vec_from_array (io[i]) * pre_gain;
      wetv           = vec_clamp (wetv, -clip_value(), clip_value());
      ducker_gain[i] = _ducker.tick (wetv);
      io[i]          = vec_to_array (wetv);
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

    std::visit (
      [&, this] (auto& algo_engine) {
        using engine = typename std::decay_t<decltype (algo_engine)>;
        using algo   = typename engine::algorithm;
        if constexpr (std::is_same_v<fixpt_t, typename engine::value_type>) {
          algo::process_block (
            algo_engine,
            _lfo,
            xspan {wet.data(), io.size()},
            pars,
            _param.algo,
            _srate);
        }
      },
      _algorithms);

    // float conversion
    auto postgain = _param.gain;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0].to_floatp();
      auto r = wet[i][1].to_floatp();
      l *= ducker_gain[i][0] * postgain;
      r *= ducker_gain[i][1] * postgain;
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
    std::array<f32_x2, max_block_size>                      ducker_gain;
    detail::lofiverb::algorithm::smoothed_parameters<float> pars;
    pars;

    auto pre_gain = 1.f / _param.gain;

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto wet       = vec_from_array (io[i]) * pre_gain;
      wet            = vec_clamp (wet, -clip_value(), clip_value());
      ducker_gain[i] = _ducker.tick (wet);
      io[i]          = vec_to_array (wet);
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

    std::visit (
      [&, this] (auto& algo_engine) {
        using engine = typename std::decay_t<decltype (algo_engine)>;
        using algo   = typename engine::algorithm;
        if constexpr (std::is_same_v<float, typename engine::value_type>) {
          algo::process_block (
            algo_engine, _lfo, io, pars, _param.algo, _srate);
        }
      },
      _algorithms);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = io[i][0] * ducker_gain[i][0] * _param.gain;
      auto r = io[i][1] * ducker_gain[i][1] * _param.gain;
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
  struct mode {
    enum {
#ifndef LOFIVERB_DEBUG_ALGO
      abyss_fix16,
      abyss_flt16,
      abyss_flt32,
      plate1_fix16,
      plate1_flt16,
      plate1_flt32,
      room_fix16,
      room_flt16,
      room_flt32,
      midifex49_fix16,
      midifex49_flt16,
      midifex49_flt32,
      midifex50_fix16,
      midifex50_flt16,
      midifex50_flt32,
      dre2000a_fix16,
      dre2000a_flt16,
      dre2000a_flt32,
      dre2000b_fix16,
      dre2000b_flt16,
      dre2000b_flt32,
      dre2000c_fix16,
      dre2000c_flt16,
      dre2000c_flt32,
      dre2000d_fix16,
      dre2000d_flt16,
      dre2000d_flt32,
      rev5_l_hall_fix16,
      rev5_l_hall_flt16,
      rev5_l_hall_flt32,
#else
      debug_algo_fix16,
      debug_algo_flt16,
      debug_algo_flt32,
#endif
    };
  };
  //----------------------------------------------------------------------------
  void update_operation_mode()
  {
    /*
      srate t44k t48k
      7200 49   20
      8400 21   40
      9000 49   16
      10500 21   32
      10800 49   40
      13500 49   32
      14400 49   10
      16800 21   20
      18000 49   08
      20160 35   50
      21000 21   16
      21600 49   20
      22500 49   32
      25200 07   40
      27000 49   16
      28800 49   05
      31500 07   32
      32400 49   40
      33600 21   10
      36000 49   04
      39200 09   60
      39600 49   40
      40320 35   25
      40500 49   32
      42000 21   08
      43200 49   10
      45000 50   16
      46800 52   40
      47040 16   50
      49000 10   49
      49500 55   33
      50400 08   21
      52500 25   35
      54000 60   09
      57600 64   06
      58500 65   39
      58800 04   49
      60480 48   63
      61200 68   51
      63000 10   21
      64800 72   27
      67200 32   07
      67500 75   45
      68400 76   57
    */
    uint srate {};

    switch (_param.mode) {
#ifndef LOFIVERB_DEBUG_ALGO
    case mode::abyss_fix16:
    case mode::abyss_flt16:
    case mode::abyss_flt32:
    case mode::plate1_fix16:
    case mode::plate1_flt16:
    case mode::plate1_flt32: {
      constexpr auto srates
        = make_array (10500, 16800, 21000, 25200, 31500, 40320, 50400);
      srate = srates[_param.srateid];
    } break;
    case mode::midifex49_fix16:
    case mode::midifex49_flt16:
    case mode::midifex49_flt32:
    case mode::midifex50_fix16:
    case mode::midifex50_flt16:
    case mode::midifex50_flt32:
    case mode::room_fix16:
    case mode::room_flt16:
    case mode::room_flt32: {
      constexpr auto srates
        = make_array (10500, 16800, 21000, 23400, 33600, 42000, 50400);
      srate = srates[_param.srateid];
    } break;
    case mode::dre2000a_fix16:
    case mode::dre2000a_flt16:
    case mode::dre2000a_flt32:
    case mode::dre2000b_fix16:
    case mode::dre2000b_flt16:
    case mode::dre2000b_flt32:
    case mode::dre2000c_fix16:
    case mode::dre2000c_flt16:
    case mode::dre2000c_flt32:
    case mode::dre2000d_fix16:
    case mode::dre2000d_flt16:
    case mode::dre2000d_flt32: {
      constexpr auto srates
        = make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
      srate = srates[_param.srateid];
    } break;
    case mode::rev5_l_hall_fix16:
    case mode::rev5_l_hall_flt16:
    case mode::rev5_l_hall_flt32: {
      constexpr auto srates
        = make_array (25200, 33600, 40320, 44100, 45000, 52500, 57600);
      srate = srates[_param.srateid];
    } break;
#else
    case mode::debug_algo_fix16:
    case mode::debug_algo_flt16:
    case mode::debug_algo_flt32: {
      constexpr auto srates
        = make_array (10500, 16800, 21000, 23400, 33600, 42000, 50400);
      srate = srates[_param.srateid];
    } break;
#endif
    default:
      return;
    }
    if (_srate != srate) {
      _srate = srate;
      update_internal_srate (srate, (srate / 20) * 9);
    }

    switch (_param.mode) {
#ifndef LOFIVERB_DEBUG_ALGO
    case mode::abyss_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::abyss::engine<dt_fix16>>();
      detail::lofiverb::abyss::reset_lfo_phase (_lfo);
    } break;
    case mode::abyss_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::abyss::engine<dt_flt16>>();
      detail::lofiverb::abyss::reset_lfo_phase (_lfo);
    } break;
    case mode::abyss_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::abyss::engine<dt_flt32>>();
      detail::lofiverb::abyss::reset_lfo_phase (_lfo);
    } break;
    case mode::plate1_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::plate1::engine<dt_fix16>>();
      detail::lofiverb::plate1::reset_lfo_phase (_lfo);
    } break;
    case mode::plate1_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::plate1::engine<dt_flt16>>();
      detail::lofiverb::plate1::reset_lfo_phase (_lfo);
    } break;
    case mode::plate1_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::plate1::engine<dt_flt32>>();
      detail::lofiverb::plate1::reset_lfo_phase (_lfo);
    } break;
    case mode::room_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::room::engine<dt_fix16>>();
      detail::lofiverb::room::reset_lfo_phase (_lfo);
    } break;
    case mode::room_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::room::engine<dt_flt16>>();
      detail::lofiverb::room::reset_lfo_phase (_lfo);
    } break;
    case mode::room_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::room::engine<dt_flt32>>();
      detail::lofiverb::room::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex49_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex49::engine<dt_fix16>>();
      detail::lofiverb::midifex49::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex49_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex49::engine<dt_flt16>>();
      detail::lofiverb::midifex49::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex49_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex49::engine<dt_flt32>>();
      detail::lofiverb::midifex49::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex50_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex50::engine<dt_fix16>>();
      detail::lofiverb::midifex50::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex50_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex50::engine<dt_flt16>>();
      detail::lofiverb::midifex50::reset_lfo_phase (_lfo);
    } break;
    case mode::midifex50_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::midifex50::engine<dt_flt32>>();
      detail::lofiverb::midifex50::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000a_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000a::engine<dt_fix16>>();
      detail::lofiverb::dre2000a::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000a_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000a::engine<dt_flt16>>();
      detail::lofiverb::dre2000a::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000a_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000a::engine<dt_flt32>>();
      detail::lofiverb::dre2000a::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000b_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000b::engine<dt_fix16>>();
      detail::lofiverb::dre2000b::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000b_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000b::engine<dt_flt16>>();
      detail::lofiverb::dre2000b::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000b_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000b::engine<dt_flt32>>();
      detail::lofiverb::dre2000b::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000c_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000c::engine<dt_fix16>>();
      detail::lofiverb::dre2000c::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000c_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000c::engine<dt_flt16>>();
      detail::lofiverb::dre2000c::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000c_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000c::engine<dt_flt32>>();
      detail::lofiverb::dre2000c::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000d_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000d::engine<dt_fix16>>();
      detail::lofiverb::dre2000d::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000d_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000d::engine<dt_flt16>>();
      detail::lofiverb::dre2000d::reset_lfo_phase (_lfo);
    } break;
    case mode::dre2000d_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::dre2000d::engine<dt_flt32>>();
      detail::lofiverb::dre2000d::reset_lfo_phase (_lfo);
    } break;
    case mode::rev5_l_hall_fix16: {
      auto& rev = _algorithms
                    .emplace<detail::lofiverb::rev5_l_hall::engine<dt_fix16>>();
      detail::lofiverb::rev5_l_hall::reset_lfo_phase (_lfo);
    } break;
    case mode::rev5_l_hall_flt16: {
      auto& rev = _algorithms
                    .emplace<detail::lofiverb::rev5_l_hall::engine<dt_flt16>>();
      detail::lofiverb::rev5_l_hall::reset_lfo_phase (_lfo);
    } break;
    case mode::rev5_l_hall_flt32: {
      auto& rev = _algorithms
                    .emplace<detail::lofiverb::rev5_l_hall::engine<dt_flt32>>();
      detail::lofiverb::rev5_l_hall::reset_lfo_phase (_lfo);
    } break;
#else
    case mode::debug_algo_fix16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::debug::engine<dt_fix16>>();
      detail::lofiverb::debug::reset_lfo_phase (_lfo);
    } break;
    case mode::debug_algo_flt16: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::debug::engine<dt_flt16>>();
      detail::lofiverb::debug::reset_lfo_phase (_lfo);
    } break;
    case mode::debug_algo_flt32: {
      auto& rev
        = _algorithms.emplace<detail::lofiverb::debug::engine<dt_flt32>>();
      detail::lofiverb::debug::reset_lfo_phase (_lfo);
    } break;
#endif
    default:
      return;
    }
    std::visit (
      [&] (auto& rev) { rev.reset_memory (_mem_reverb); }, _algorithms);
    _n_processed_samples = 0; // trigger the control block on first sample
    update_mod();
  }
  //----------------------------------------------------------------------------
  void update_internal_srate (uint srate, uint fc)
  {
    auto daw_srate = _pc->get_sample_rate();
    auto target_fc = std::min (fc, (daw_srate / 20) * 9);
    _resampler.reset (
      srate,
      daw_srate,
      fc,
      target_fc,
      32,
      16,
      210,
      true,
      max_block_size,
      6 * 1024);

    auto  state   = _pc->get_play_state();
    float beat_hz = 120.f;
    if (state.is_valid) {
      // playhead may not exist, eg. when scanning the plugin
      beat_hz = state.bpm * (1.f / 60.f);
    }
    _1_4beat_spls = (0.25f / beat_hz) * srate;
    _t_spl        = 1.f / srate;

    _ducker.reset();
    _lfo.reset();
    _param_smooth.reset_srate (_t_spl);

    // resize memory
    uint predelay_spls = std::ceil (_1_4beat_spls * max_predelay_qb) * 2;
    uint rev_bytes     = 0;
    mp11::mp_for_each<decltype (_algorithms)> ([&] (auto mode) {
      rev_bytes = std::max (mode.get_required_bytes(), rev_bytes);
    });
    uint predelay_flt = predelay_spls * sizeof (float);
    uint padding      = round_ceil (predelay_flt, 16u) - predelay_flt;

    _mem.clear();
    _mem.resize (predelay_flt + padding + rev_bytes);

    _mem_reverb = xspan {_mem};
    _predelay.reset (_mem_reverb.cut_head (predelay_flt).cast<float>(), 2);
    _mem_reverb.cut_head (padding);

    update_mod();
    update_ducker();
  }
  //----------------------------------------------------------------------------
  void update_ducker()
  {
    auto v = _param.ducking_speed * 0.01f;
    _ducker.set_speed (vec_set<2> (v * v), _t_spl);
  }
  //----------------------------------------------------------------------------
  void update_mod()
  {
    std::visit (
      [this] (auto& algo_engine) {
        using algo = typename std::decay_t<decltype (algo_engine)>::algorithm;
        algo::reset_lfo_freq (_lfo, _param_smooth.target().mod, _t_spl);
      },
      _algorithms);
  }
  //----------------------------------------------------------------------------
  static float clip_value() { return std::min (fixpt_max_flt, float16::max()); }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32                                                mode;
    u32                                                srateid;
    float                                              gain; // a parameter
    float                                              predelay;
    detail::lofiverb::algorithm::unsmoothed_parameters algo;
    float                                              ducking_threshold;
    float                                              ducking_speed;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float mod;
    float stereo;
    float decay;
    float character;
    // dry, wet, ducker/gate
  };
  //----------------------------------------------------------------------------
  using fixpt_t  = detail::lofiverb::fixpt_t;
  using fixpt_tr = detail::lofiverb::fixpt_tr;
  using float16  = detail::lofiverb::float16;

  static constexpr auto dt_fix16 = detail::lofiverb::delay::data_type::fixpt16;
  static constexpr auto dt_flt16 = detail::lofiverb::delay::data_type::float16;
  static constexpr auto dt_flt32 = detail::lofiverb::delay::data_type::float32;

  static constexpr auto fixpt_max_flt = detail::lofiverb::fixpt_max_flt;
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>             _resampler {};
  static_delay_line<float, true, false> _predelay;

  lfo<4> _lfo;

  using algorithms_type = std::variant<
#ifndef LOFIVERB_DEBUG_ALGO
    detail::lofiverb::abyss::engine<dt_fix16>,
    detail::lofiverb::abyss::engine<dt_flt16>,
    detail::lofiverb::abyss::engine<dt_flt32>,
    detail::lofiverb::plate1::engine<dt_fix16>,
    detail::lofiverb::plate1::engine<dt_flt16>,
    detail::lofiverb::plate1::engine<dt_flt32>,
    detail::lofiverb::room::engine<dt_fix16>,
    detail::lofiverb::room::engine<dt_flt16>,
    detail::lofiverb::room::engine<dt_flt32>,
    detail::lofiverb::midifex49::engine<dt_fix16>,
    detail::lofiverb::midifex49::engine<dt_flt16>,
    detail::lofiverb::midifex49::engine<dt_flt32>,
    detail::lofiverb::midifex50::engine<dt_fix16>,
    detail::lofiverb::midifex50::engine<dt_flt16>,
    detail::lofiverb::midifex50::engine<dt_flt32>,
    detail::lofiverb::dre2000a::engine<dt_fix16>,
    detail::lofiverb::dre2000a::engine<dt_flt16>,
    detail::lofiverb::dre2000a::engine<dt_flt32>,
    detail::lofiverb::dre2000b::engine<dt_fix16>,
    detail::lofiverb::dre2000b::engine<dt_flt16>,
    detail::lofiverb::dre2000b::engine<dt_flt32>,
    detail::lofiverb::dre2000c::engine<dt_fix16>,
    detail::lofiverb::dre2000c::engine<dt_flt16>,
    detail::lofiverb::dre2000c::engine<dt_flt32>,
    detail::lofiverb::dre2000d::engine<dt_fix16>,
    detail::lofiverb::dre2000d::engine<dt_flt16>,
    detail::lofiverb::dre2000d::engine<dt_flt32>,
    detail::lofiverb::rev5_l_hall::engine<dt_fix16>,
    detail::lofiverb::rev5_l_hall::engine<dt_flt16>,
    detail::lofiverb::rev5_l_hall::engine<dt_flt32>
#else
    detail::lofiverb::debug::engine<dt_fix16>,
    detail::lofiverb::debug::engine<dt_flt16>,
    detail::lofiverb::debug::engine<dt_flt32>
#endif
    >;
  algorithms_type _algorithms;
  ducker<f32_x2>  _ducker;

  uint  _n_processed_samples;
  float _1_4beat_spls;
  float _t_spl;
  uint  _srate;

  plugin_context* _pc;

  std::vector<u8> _mem;
  xspan<u8>       _mem_reverb;
};
//------------------------------------------------------------------------------
} // namespace artv
