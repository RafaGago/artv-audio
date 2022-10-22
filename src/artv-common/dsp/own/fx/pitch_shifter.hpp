#pragma once

#include <algorithm>

#include <gcem.hpp>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/pitch_shift.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
class pitch_shifter {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::pitch;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct mode_tag {};

  void set (mode_tag, uint v)
  {
    if (v == _mode) {
      return;
    }
    _mode = v;
    reset_delay();
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      1,
      make_cstr_array (
        "SR/2",
        "SR/4",
        "SR/8",
        "SR/16",
        "SR/32",
        "SR/64",
        "SR/128",
        "SR/256",
        "SR/512"),
      32);
  }
  //----------------------------------------------------------------------------
  struct semitones_tag {};

  void set (semitones_tag, float v)
  {
    if (v == _amt) {
      return;
    }
    _amt = v;
    reset_amt();
  }

  static constexpr auto get_parameter (semitones_tag)
  {
    return float_param ("Semitones", -24.0, 24.0, 0.0, 1.);
  }
  //----------------------------------------------------------------------------
  static constexpr double detune_resolution = 0.001;

  struct detune_tag {};

  void set (detune_tag, float v)
  {
    v = (v < detune_resolution && v > -detune_resolution) ? 0.f : v;
    if (v == _detune) {
      return;
    }
    _detune = v;
    reset_amt (false);
  }

  static constexpr auto get_parameter (detune_tag)
  {
    return float_param ("Semitones", -12.0, 12.0, 0.0, detune_resolution);
  }
  //----------------------------------------------------------------------------
  struct count_tag {};

  void set (count_tag, uint v)
  {
    if (v == _n_shift) {
      return;
    }
    _n_shift = v;
    reset_amt();
    reset_factors();
  }

  static constexpr auto get_parameter (count_tag)
  {
    return int_param ("", 1, max_readers, 1);
  }
  //----------------------------------------------------------------------------
  struct width_tag {};

  void set (width_tag, float v)
  {
    decltype (_width) det = v * 0.01;
    if (det == _width) {
      return;
    }
    _width = det;
    reset_factors();
  }

  static constexpr auto get_parameter (width_tag)
  {
    return float_param ("%", -100.0, 100.0, 0.0, 0.1);
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<mode_tag, semitones_tag, count_tag, detune_tag, width_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _mem.clear();
    _mem.resize (get_delay_size (0));

    _mode    = 1;
    _n_shift = 1;
    _width = _amt = _detune = 0.f;

    reset_delay();
    reset_factors();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < samples; ++i) {
      _shift.push (vec<float, 2> {ins[0][i], ins[1][i]});
      vec<float, 2> out {};
      for (uint p = 0; p < _n_shift; ++p) {
        out += _shift.read (_readers[p]) * _factors[p];
      }
      outs[0][i] = out[0];
      outs[1][i] = out[1];
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  uint get_delay_size (uint mode)
  {
    uint delay_samples = _plugcontext->get_sample_rate() / 2;
    uint delay_size = 1u << (last_bit_set (delay_samples) - 1); // floor (log2)
    delay_size      = delay_size < delay_samples ? delay_size << 1 : delay_size;
    return delay_size >> mode;
  }
  //----------------------------------------------------------------------------
  void reset_delay()
  {
    uint size = get_delay_size (_mode);
    _plugcontext->set_delay_compensation (size / 2);
    assert (size <= _mem.size());
    _shift.reset (xspan {_mem.data(), size});
    reset_amt(); // reset all readers
  }
  //----------------------------------------------------------------------------
  void reset_amt (bool resync = false)
  {
    float val = _amt + _detune;
    float dec = _detune / ((float) _n_shift - 1.f);

    for (uint i = 0; i < _n_shift; ++i) {
      _shift.set_reader (_readers[i], val, resync);
      val -= dec;
    }
  }
  //----------------------------------------------------------------------------
  void reset_factors()
  {
    // A power based pan law skewing radically towards the extremes
    static constexpr float pan_law_power = gcem::log (1.f / 16.f);

    float n_recip = 1.f / (float) _n_shift;
    float step    = _width / (float) (_n_shift - 1);
    float pan     = 0.5f - _width * 0.5f;

    for (uint i = 0; i < _n_shift; ++i) {
      _factors[i][0] = n_recip * expf (pan * pan_law_power);
      _factors[i][1] = n_recip * expf ((1.f - pan) * pan_law_power);
      pan += step;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_readers = 12;

  float _width;
  float _detune;
  float _amt;
  uint  _n_shift;
  uint  _mode;

  pitch_shift_sin<vec<float, 2>, true>               _shift;
  std::array<decltype (_shift)::reader, max_readers> _readers;
  std::array<vec<float, 2>, max_readers>             _factors;

  std::vector<vec<float, 2>> _mem;

  plugin_context* _plugcontext;
};
//------------------------------------------------------------------------------

} // namespace artv
