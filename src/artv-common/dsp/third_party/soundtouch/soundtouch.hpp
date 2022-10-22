#pragma once

#include <algorithm>
#include <cmath>
#include <optional>

#include "SoundTouch.h"

#include "artv-common/dsp/own/classes/misc.hpp"
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
class soundtouch_fx {
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
    switch (_mode) {
    case 0:
      _shift->setSetting (SETTING_USE_AA_FILTER, 0);
      break;
    case 1:
      _shift->setSetting (SETTING_USE_AA_FILTER, 1);
      _shift->setSetting (SETTING_AA_FILTER_LENGTH, 8);
      break;
    case 2:
      _shift->setSetting (SETTING_USE_AA_FILTER, 1);
      _shift->setSetting (SETTING_AA_FILTER_LENGTH, 32);
      break;
    case 3:
      _shift->setSetting (SETTING_USE_AA_FILTER, 1);
      _shift->setSetting (SETTING_AA_FILTER_LENGTH, 128);
      break;
    default:
      return;
    }
    fix_latency();
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Normal", "AA-8", "AA-32", "AA-128"), 16);
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
    return float_param ("Semitones", -36.0, 36.0, 0.0, 1.);
  }
  //----------------------------------------------------------------------------
  struct detune_tag {};

  void set (detune_tag, float v)
  {
    if (v == _detune) {
      return;
    }
    _detune = v;
    reset_amt();
  }

  static constexpr auto get_parameter (detune_tag)
  {
    return float_param ("Semitones", -1.0, 1.0, 0.0, 0.01);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<semitones_tag, detune_tag, mode_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;

    _shift.emplace();
    _shift->setChannels (2);
    _shift->setSampleRate (pc.get_sample_rate());

    _latency = 0;
    _mode = _overlap = 1000;
    _ratio           = 1.f;
    _amt = _detune = -10.f;

    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void process (xspan<float*> outs, xspan<float const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    std::array<std::array<float, 2>, 128> buff;

    for (uint offset = 0; offset < samples; offset += buff.size()) {
      uint blocksize = std::min<uint> (buff.size(), samples - offset);

      // interleave
      for (uint i = 0; i < blocksize; ++i) {
        buff[i][0] = ins[0][offset + i];
        buff[i][1] = ins[1][offset + i];
      }

      _shift->putSamples (buff.data()->data(), blocksize);
      uint got = _shift->receiveSamples (buff.data()->data(), blocksize);

      // deinterleave
      if (likely (got == blocksize)) {
        for (uint i = 0; i < blocksize; ++i) {
          outs[0][offset + i] = buff[i][0];
          outs[1][offset + i] = buff[i][1];
        }
      }
      else {
        uint blank = blocksize - got;
        memset (&outs[0][offset], 0, sizeof outs[0][offset] * blank);
        memset (&outs[1][offset], 0, sizeof outs[1][offset] * blank);
        for (uint i = blank; i < blocksize; ++i) {
          outs[0][offset + i] = buff[i][0];
          outs[1][offset + i] = buff[i][1];
        }
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void reset_amt()
  {
    _ratio = exp ((_amt + _detune) * (1. / 12.) * M_LN2);
    _shift->setPitch (_ratio);
    fix_latency();
  }
  //----------------------------------------------------------------------------
  void fix_latency()
  {
    // The soundtouch author(s) should really add a prefix to its macros, or if
    // the library is CPP only use them for things that make sense (not
    // constants).
    auto lat = _shift->getSetting (SETTING_INITIAL_LATENCY);
    if (lat != _latency) {
      _latency = lat;
      _plugcontext->set_delay_compensation (lat);
    }
  }
  //----------------------------------------------------------------------------
  float _detune;
  float _amt;
  float _ratio;

  uint _mode;
  uint _overlap;
  uint _latency;

  // it is default constructible, but it doesn't take samplerate changes well
  std::optional<soundtouch::SoundTouch> _shift;
  plugin_context*                       _plugcontext;
};
//------------------------------------------------------------------------------

} // namespace artv
