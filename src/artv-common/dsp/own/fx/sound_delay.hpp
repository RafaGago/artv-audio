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
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
class sound_delay {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::delay;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct delay_ms_tag {};

  void set (delay_ms_tag, float v)
  {
    if (_delay_ms == v) {
      return;
    }
    _delay_ms = v;
    update_delay_times();
  }

  static constexpr auto get_parameter (delay_ms_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_ms_l_offset_tag {};

  void set (delay_ms_l_offset_tag, float v)
  {
    if (_delay_ms_l == v) {
      return;
    }
    _delay_ms_l = v;
    update_delay_times();
  }

  static constexpr auto get_parameter (delay_ms_l_offset_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_ms_r_offset_tag {};

  void set (delay_ms_r_offset_tag, float v)
  {
    if (_delay_ms_r == v) {
      return;
    }
    _delay_ms_r = v;
    update_delay_times();
  }

  static constexpr auto get_parameter (delay_ms_r_offset_tag)
  {
    return float_param ("ms", 0.0, 500, 0.0, 0.001);
  }
  //----------------------------------------------------------------------------
  struct delay_beats_l_tag {};

  void set (delay_beats_l_tag, float v)
  {
    v *= (1. / 16.);
    if (_delay_beats_l == v) {
      return;
    }
    _delay_beats_l = v;
    update_delay_times();
  }

  static constexpr auto get_parameter (delay_beats_l_tag)
  {
    return float_param ("Sixteenths", 0.0, 32, 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct delay_beats_r_tag {};

  void set (delay_beats_r_tag, float v)
  {
    v *= (1. / 16.);
    if (_delay_beats_r == v) {
      return;
    }
    _delay_beats_r = v;
    update_delay_times();
  }

  static constexpr auto get_parameter (delay_beats_r_tag)
  {
    return float_param ("Sixteenths", 0.0, 32, 0., 1.);
  }
  //----------------------------------------------------------------------------
  struct delay_samples_tag {};

  void set (delay_samples_tag, int v)
  {
    if (_delay_samples == v) {
      return;
    }
    _delay_samples = v;
    update_delay_times();
  }

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
      get_delay_samples (get_parameter (delay_ms_tag {}).max));

    uint _max_tempo_delay = (uint) std::ceil (get_synced_delay_samples (
      get_parameter (delay_beats_l_tag {}).max * (1. / 16.)));

    _delay.reset (2, std::max (_max_delay_samples, _max_tempo_delay));
    _delay_times_samples[0] = _delay_times_samples[1] = 0.;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < block_samples; ++i) {
      auto t_samples = _delay_times_samples;
      auto in        = make_array (ins[0][i], ins[1][i]);
      _delay.push (in);
      outs[0][i] = _delay.get (t_samples[0], 0);
      outs[1][i] = _delay.get (t_samples[1], 1);
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
      = std::min (_delay_times_samples[0], (double) _delay.size());

    _delay_times_samples[1] += get_delay_samples (_delay_ms_r);
    _delay_times_samples[1] += get_synced_delay_samples (_delay_beats_r);
    _delay_times_samples[1]
      = std::min (_delay_times_samples[1], (double) _delay.size());
  }
  //----------------------------------------------------------------------------
  float _delay_ms;
  float _delay_ms_l;
  float _delay_ms_r;
  float _delay_beats_l;
  float _delay_beats_r;
  float _delay_samples;

  std::array<double, 2>           _delay_times_samples;
  dynamic_delay_line<float, true> _delay;
  plugin_context*                 _plugcontext = nullptr;
};
//------------------------------------------------------------------------------
}; // namespace artv
