#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <utility>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

#define SOUND_DELAY_LEAN 1

//------------------------------------------------------------------------------
class sound_delay {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::delay;
  //----------------------------------------------------------------------------
  struct delay_ms_tag {};

  void set (delay_ms_tag, float v) { _delay_ms = v; }

  static constexpr auto get_parameter (delay_ms_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_ms_l_offset_tag {};

  void set (delay_ms_l_offset_tag, float v) { _delay_ms_l = v; }

  static constexpr auto get_parameter (delay_ms_l_offset_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_ms_r_offset_tag {};

  void set (delay_ms_r_offset_tag, float v) { _delay_ms_r = v; }

  static constexpr auto get_parameter (delay_ms_r_offset_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_beats_l_tag {};

  void set (delay_beats_l_tag, float v) { _delay_beats_l = v * (1. / 16.); }

  static constexpr auto get_parameter (delay_beats_l_tag)
  {
    return float_param ("Sixteenths", 0.0, 32, 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct delay_beats_r_tag {};

  void set (delay_beats_r_tag, float v) { _delay_beats_r = v * (1. / 16.); }

  static constexpr auto get_parameter (delay_beats_r_tag)
  {
    return float_param ("Sixteenths", 0.0, 32, 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct delay_samples_tag {};

  void set (delay_samples_tag, int v) { _delay_samples = v; }

  static constexpr auto get_parameter (delay_samples_tag)
  {
    constexpr auto max_delay = 1024.f;
    static_assert (
      ((get_parameter (delay_ms_tag {}).max * 44100.)) >= max_delay, "");
    return float_param ("Spls", 0., max_delay, 0., 1.);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    delay_ms_tag,
    delay_ms_l_offset_tag,
    delay_ms_r_offset_tag,
    delay_beats_l_tag,
    delay_beats_r_tag,
    delay_samples_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;

    uint _max_delay_samples = (uint) std::ceil (
      ((float) pc.get_sample_rate()) * get_parameter (delay_ms_tag {}).max);

    uint _max_tempo_delay = (uint) std::ceil (std::ceil (
      pc.get_samples_per_beat() * get_parameter (delay_beats_l_tag {}).max
      * (1. / 16.)));

    _delay.reset (2, std::max (_max_delay_samples, _max_tempo_delay));

#if !SOUND_DELAY_LEAN
    onepole_smoother<double>::lowpass (
      make_crange (_smooth_coeff), 1 / 0.015, pc.get_sample_rate());
    _smooth_state[0] = _smooth_state[1] = 0.;
#endif
    _delay_times_samples[0] = _delay_times_samples[1] = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint block_samples)
  {
    update_delay_times();

    for (uint i = 0; i < block_samples; ++i) {
#if SOUND_DELAY_LEAN
      auto t_samples = _delay_times_samples;
#else
      auto t_samples = onepole_smoother<double>::tick (
        make_crange (_smooth_coeff),
        make_array (
          make_crange (_smooth_state[0]), make_crange (_smooth_state[1])),
        _delay_times_samples);
#endif
      auto in = make_array (chnls[0][i], chnls[1][i]);
      _delay.push (in);

      chnls[0][i] = _delay.get (0, t_samples[0]);
      chnls[1][i] = _delay.get (1, t_samples[1]);
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  float get_delay_samples (float ms)
  {
    return _plugcontext->get_sample_rate() * ms * 0.001;
  }
  //----------------------------------------------------------------------------
  float get_synced_delay_samples (float beats)
  {
    return _plugcontext->get_samples_per_beat() * beats;
  }
  //----------------------------------------------------------------------------
  void update_delay_times()
  {
    _delay_times_samples[0]
      = get_delay_samples (_delay_ms) + (float) _delay_samples;
    _delay_times_samples[1] = _delay_times_samples[0];

    _delay_times_samples[0] += get_delay_samples (_delay_ms_l);
    _delay_times_samples[0] += get_synced_delay_samples (_delay_beats_l);
    _delay_times_samples[0]
      = std::min (_delay_times_samples[0], (double) _delay.max_delay_samples());

    _delay_times_samples[1] += get_delay_samples (_delay_ms_r);
    _delay_times_samples[1] += get_synced_delay_samples (_delay_beats_r);
    _delay_times_samples[1]
      = std::min (_delay_times_samples[1], (double) _delay.max_delay_samples());
  }
  //----------------------------------------------------------------------------
  float _delay_ms;
  float _delay_ms_l;
  float _delay_ms_r;
  float _delay_beats_l;
  float _delay_beats_r;
  float _delay_samples;

  std::array<double, 2> _delay_times_samples;
#if !SOUND_DELAY_LEAN
  std::array<double, 2> _smooth_state;
  double                _smooth_coeff;
#endif
  delay_line<float> _delay;
  plugin_context*   _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
}; // namespace artv
