#pragma once

#include <algorithm>
#include <cmath>
#include <optional>

#include "rubberband/RubberBandStretcher.h"

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
class rubberband {
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
      _shift->setPhaseOption (rbs::OptionPhaseIndependent);
      _shift->setTransientsOption (rbs::OptionTransientsSmooth);
      break;
    case 1:
      _shift->setPhaseOption (rbs::OptionPhaseLaminar);
      _shift->setTransientsOption (rbs::OptionTransientsSmooth);
      break;
    case 2:
      _shift->setPhaseOption (rbs::OptionPhaseLaminar);
      _shift->setTransientsOption (rbs::OptionTransientsMixed);
      break;
    case 3:
      _shift->setPhaseOption (rbs::OptionPhaseLaminar);
      _shift->setTransientsOption (rbs::OptionTransientsCrisp);
      break;
    default:
      return;
    }
    fix_latency();
  }

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0,
      make_cstr_array ("Smooth", "Multitimbral", "Two sources", "Standard"),
      16);
  }
  //----------------------------------------------------------------------------
  struct formant_mode_tag {};

  void set (formant_mode_tag, uint v)
  {
    if (v == _formant) {
      return;
    }
    _formant = v;
    auto opt = (_formant == 0) ? rbs::OptionFormantShifted
                               : rbs::OptionFormantPreserved;
    _shift->setFormantOption (opt);
    fix_latency();
  }

  static constexpr auto get_parameter (formant_mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Formant Shifted", "Formant Preserved"), 16);
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
  using parameters
    = mp_list<mode_tag, formant_mode_tag, semitones_tag, detune_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;

    auto flags = rbs::OptionProcessRealTime | rbs::OptionPitchHighConsistency;
    _shift.emplace (pc.get_sample_rate(), 2, flags);
    _shift->reset();

    _latency = 0;
    _mode = _formant = 1000;
    _ratio           = 1.f;
    _amt = _detune = -10.f;

    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void process (crange<float*> outs, crange<float const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    uint blocksize;
    for (uint offset = 0; offset < samples; offset += blocksize) {
      blocksize
        = std::min<uint> (_shift->getSamplesRequired(), samples - offset);

      std::array<float const*, 2> in  = {ins[0] + offset, ins[1] + offset};
      std::array<float*, 2>       out = {outs[0] + offset, outs[1] + offset};

      _shift->process (in.data(), blocksize, false);

      // maybe add blank samples while waiting for the latency to pass.
      uint got = std::min<uint> (_shift->available(), blocksize);

      if (unlikely (got != blocksize)) {
        uint blank = blocksize - got;
        for (uint i = 0; i < 2; ++i) {
          crange_memset (make_crange (out[i], blank), 0);
          out[i] += blank;
        }
      }
      _shift->retrieve (out.data(), got);
    }
  }
  //----------------------------------------------------------------------------
private:
  using rbs = RubberBand::RubberBandStretcher;
  //----------------------------------------------------------------------------
  void reset_amt()
  {
    _ratio = exp ((_amt + _detune) * (1. / 12.) * M_LN2);
    _shift->setPitchScale (_ratio);
    fix_latency();
  }
  //----------------------------------------------------------------------------
  void fix_latency()
  {
    auto lat = _shift->getLatency();
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
  uint _formant;
  uint _latency;

  std::optional<rbs> _shift;
  plugin_context*    _plugcontext;
};
//------------------------------------------------------------------------------

} // namespace artv
