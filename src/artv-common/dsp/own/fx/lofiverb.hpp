#pragma once

// #define LOFIVERB_DEBUG_ALGO 1

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
#ifndef LOFIVERB_DEBUG_ALGO
#include "artv-common/dsp/own/fx/lofiverb/acreil-algos.hpp"
#include "artv-common/dsp/own/fx/lofiverb/artv-algos.hpp"
#else
#include "artv-common/dsp/own/fx/lofiverb/debug-algo.hpp"
#endif
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

namespace artv {

// TODO list
// - dre algorithms are too dark now
// - dynamics at the output compressor/expander.

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
    update_algorithm();
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (algorithm_tag)
  {
    return choice_param (
      0,
      make_cstr_array (
#ifndef LOFIVERB_DEBUG_ALGO
        "Artv Ambience",
        "Artv Room",
        "Artv Hall",
        "Artv Broken Hall",
        "Artv Arena",
        "Artv Palace",
        "Artv FlutterPlate",
        "Artv Abyss",
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
      128);
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
    update_algorithm();
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
    update_algorithm();
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
  void set (decay_tag, float v) { _param_smooth.target().decay = v * 0.01f; }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mod_tag {};
  void set (mod_tag, float v)
  {
    v *= 0.01f;
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
    _param_smooth.target().character = v * 0.01f;
  }

  static constexpr auto get_parameter (character_tag)
  {
    return float_param ("%", 0., 100., 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct lf_amt_tag {};
  void set (lf_amt_tag, float v) { _param.algo.lf_amt = v * 0.01f; }

  static constexpr auto get_parameter (lf_amt_tag)
  {
    return float_param ("%", 0., 100., 80., 0.01);
  }
  //----------------------------------------------------------------------------
  struct hf_amt_tag {};
  void set (hf_amt_tag, float v) { _param.algo.hf_amt = v * 0.01f; }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (hf_amt_tag)
  {
    return float_param ("%", 0., 100., 80., 0.01);
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_predelay_msec = 1000;
  struct predelay_tag {};
  void set (predelay_tag, float v) { _param.predelay = v * 0.001f; }

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("msec", 0., max_predelay_msec, 0., 0.1, 0.6);
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
  struct dyn_threshold_tag {};
  void set (dyn_threshold_tag, float v)
  {
    if (v == _param.ducking_threshold) {
      return;
    }
    _param.ducking_threshold = v;
    _ducker.set_threshold (vec_set<2> (v));
  }

  static constexpr auto get_parameter (dyn_threshold_tag)
  {
    return float_param ("dB", -60.f, 0.0f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct dyn_speed_tag {};
  void set (dyn_speed_tag, float v)
  {
    if (v == _param.ducking_speed) {
      return;
    }
    _param.ducking_speed = v;
    update_ducker();
  }

  static constexpr auto get_parameter (dyn_speed_tag)
  {
    return float_param ("%", -100.f, 100.f, 10.f, 0.01f);
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
        process<fixpt_t> (io);
      }
      else {
        process<float> (io);
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
    dyn_threshold_tag,
    dyn_speed_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = detail::lofiverb::max_block_size;
  static constexpr uint n_channels     = 2;
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<std::array<float, 2>> io)
  {
    using algo = detail::lofiverb::algorithm;
    assert (io.size() <= max_block_size);

    static constexpr bool is_fixpt_t = std::is_same_v<T, fixpt_t>;

    std::array<f32_x2, max_block_size>                  ducker_gain;
    detail::lofiverb::algorithm::smoothed_parameters<T> pars;
    bool gating = _param.ducking_speed < 0.f;

    auto pre_gain = (1.f / _param.gain);
    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto wetv      = vec_from_array (io[i]) * pre_gain;
      wetv           = vec_clamp (wetv, -clip_value(), clip_value());
      ducker_gain[i] = _ducker.tick (wetv);
      if (gating) {
        ducker_gain[i] = 1.f - ducker_gain[i];
      }
      io[i] = vec_to_array (wetv);
      _param_smooth.tick();
      pars.stereo[i] = _param_smooth.get().stereo;
      if constexpr (is_fixpt_t) {
        pars.decay[i].load_float (_param_smooth.get().decay);
        pars.character[i].load_float (_param_smooth.get().character);
        pars.mod[i].load_float (_param_smooth.get().mod);
      }
      else {
        pars.decay[i]     = _param_smooth.get().decay;
        pars.character[i] = _param_smooth.get().character;
        pars.mod[i]       = _param_smooth.get().mod;
      }
    }
    // predelay
    uint predelay_spls = _srate * _param.predelay;
    if (predelay_spls != 0) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        std::array<float, 2> spl;
        spl[0] = _predelay.get (predelay_spls, 0);
        spl[1] = _predelay.get (predelay_spls, 1);
        _predelay.push (xspan {io[i]});
        io[i] = spl;
      }
    }
    std::visit (
      [&, this] (auto& algo) {
        using algo_type = typename std::decay_t<decltype (algo)>;
        if constexpr (std::is_same_v<T, typename algo_type::value_type>) {
          if constexpr (is_fixpt_t) {
            array2d<fixpt_t, 2, max_block_size> wet;
            ARTV_LOOP_UNROLL_SIZE_HINT (16)
            for (uint i = 0; i < io.size(); ++i) {
              wet[i][0].load_float (io[i][0]);
              wet[i][1].load_float (io[i][1]);
            }
            algo.process_block (
              xspan {wet.data(), io.size()}, pars, _param.algo, _srate);
            ARTV_LOOP_UNROLL_SIZE_HINT (16)
            for (uint i = 0; i < io.size(); ++i) {
              io[i][0] = wet[i][0].to_floatp();
              io[i][1] = wet[i][1].to_floatp();
            }
          }
          else {
            algo.process_block (io, pars, _param.algo, _srate);
          }
          algo.post_process_block (io);
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
  void update_algorithm()
  {
    mp11::mp_for_each<mp11::mp_iota<mp11::mp_size<algorithms_type>>> (
      [&] (auto i) {
        if (_param.mode == i) {
          using algo = mp11::mp_at_c<algorithms_type, i>;
          _algorithms.emplace<algo>();
        }
      });
    std::visit (
      [&] (auto& algo) {
        uint srate = algo.get_sample_rates()[_param.srateid];
        if (_srate != srate) {
          _srate = srate;
          update_internal_srate (srate, (srate / 20) * 9);
        }
        algo.reset (_mem_reverb, _t_spl);
      },
      _algorithms);
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
    float beat_hz = 120.f * (1.f / 60.f);
    if (state.is_valid) {
      // playhead may not exist, eg. when scanning the plugin
      beat_hz = state.bpm * (1.f / 60.f);
    }
    _1_4beat_spls = (0.25f / beat_hz) * srate;
    _t_spl        = 1.f / srate;

    _ducker.reset();
    _param_smooth.reset_srate (_t_spl);

    // resize memory
    uint predelay_spls = std::ceil (srate * max_predelay_msec * 0.001f) * 2;
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
    float v = abs (_param.ducking_speed * 0.01f);
    _ducker.set_speed (vec_set<2> (v * v), _t_spl);
  }
  //----------------------------------------------------------------------------
  void update_mod()
  {
    std::visit (
      [this] (auto& algo) {
        algo.mod_changed (_param_smooth.target().mod, _t_spl);
      },
      _algorithms);
  }
  //----------------------------------------------------------------------------
  static float clip_value()
  {
    return std::min (fixpt_sto16::max_float(), float16::max());
  }
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
  using fixpt_t     = detail::lofiverb::fixpt_t;
  using fixpt_sto16 = detail::lofiverb::fixpt_sto16;
  using float16     = detail::lofiverb::float16;

  static constexpr auto dt_fix16 = detail::lofiverb::delay::data_type::fixpt16;
  static constexpr auto dt_flt16 = detail::lofiverb::delay::data_type::float16;
  static constexpr auto dt_flt32 = detail::lofiverb::delay::data_type::float32;
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>             _resampler {};
  static_delay_line<float, true, false> _predelay;

  using algorithms_type = std::variant<
#ifndef LOFIVERB_DEBUG_ALGO
    detail::lofiverb::ambience<dt_fix16>,
    detail::lofiverb::ambience<dt_flt16>,
    detail::lofiverb::ambience<dt_flt32>,
    detail::lofiverb::room<dt_fix16>,
    detail::lofiverb::room<dt_flt16>,
    detail::lofiverb::room<dt_flt32>,
    detail::lofiverb::hall<dt_fix16>,
    detail::lofiverb::hall<dt_flt16>,
    detail::lofiverb::hall<dt_flt32>,
    detail::lofiverb::broken_hall<dt_fix16>,
    detail::lofiverb::broken_hall<dt_flt16>,
    detail::lofiverb::broken_hall<dt_flt32>,
    detail::lofiverb::arena<dt_fix16>,
    detail::lofiverb::arena<dt_flt16>,
    detail::lofiverb::arena<dt_flt32>,
    detail::lofiverb::palace<dt_fix16>,
    detail::lofiverb::palace<dt_flt16>,
    detail::lofiverb::palace<dt_flt32>,
    detail::lofiverb::plate1<dt_fix16>,
    detail::lofiverb::plate1<dt_flt16>,
    detail::lofiverb::plate1<dt_flt32>,
    detail::lofiverb::abyss<dt_fix16>,
    detail::lofiverb::abyss<dt_flt16>,
    detail::lofiverb::abyss<dt_flt32>,
    detail::lofiverb::midifex49<dt_fix16>,
    detail::lofiverb::midifex49<dt_flt16>,
    detail::lofiverb::midifex49<dt_flt32>,
    detail::lofiverb::midifex50<dt_fix16>,
    detail::lofiverb::midifex50<dt_flt16>,
    detail::lofiverb::midifex50<dt_flt32>,
    detail::lofiverb::dre2000a<dt_fix16>,
    detail::lofiverb::dre2000a<dt_flt16>,
    detail::lofiverb::dre2000a<dt_flt32>,
    detail::lofiverb::dre2000b<dt_fix16>,
    detail::lofiverb::dre2000b<dt_flt16>,
    detail::lofiverb::dre2000b<dt_flt32>,
    detail::lofiverb::dre2000c<dt_fix16>,
    detail::lofiverb::dre2000c<dt_flt16>,
    detail::lofiverb::dre2000c<dt_flt32>,
    detail::lofiverb::dre2000d<dt_fix16>,
    detail::lofiverb::dre2000d<dt_flt16>,
    detail::lofiverb::dre2000d<dt_flt32>,
    detail::lofiverb::rev5_l_hall<dt_fix16>,
    detail::lofiverb::rev5_l_hall<dt_flt16>,
    detail::lofiverb::rev5_l_hall<dt_flt32>
#else
    detail::lofiverb::debug<dt_fix16>,
    detail::lofiverb::debug<dt_flt16>,
    detail::lofiverb::debug<dt_flt32>
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
