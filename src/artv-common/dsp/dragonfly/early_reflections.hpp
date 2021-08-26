#pragma once

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/dragonfly/compilation_firewalls.hpp"

namespace artv { namespace dragonfly {

class early_reflections {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::reverb;
  //----------------------------------------------------------------------------
  struct dry_tag {};
  void set (dry_tag, float v)
  {
    if (v > -60.f) {
      _dsp.set_parameter (paramDry, db_to_gain (v) * 100.f);
    }
    else {
      _dsp.set_parameter (paramDry, 0.f);
    }
  }

  static constexpr auto get_parameter (dry_tag)
  {
    // {paramDry, "Dry Level", "dry_level", 0.0f, 100.0f, "%"},
    return float_param ("dB", -60., 0., -60., 0.1);
  }
  //----------------------------------------------------------------------------
  struct wet_tag {};
  void set (wet_tag, float v)
  {
    _dsp.set_parameter (paramWet, db_to_gain (v) * 100.f);
  }

  static constexpr auto get_parameter (wet_tag)
  {
    // {paramWet, "Wet Level", "early_level", 0.0f, 100.0f, "%"},
    return float_param ("dB", -60., 0., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  struct size_tag {};
  void set (size_tag, float v) { _dsp.set_parameter (paramSize, v); }

  static constexpr auto get_parameter (size_tag)
  {
    // {paramSize, "Size", "size", 10.0f, 60.0f, "m"},
    return float_param ("m", 10., 60., 20., 0.1);
  }
  //----------------------------------------------------------------------------
  struct width_tag {};
  void set (width_tag, float v) { _dsp.set_parameter (paramWidth, v); }

  static constexpr auto get_parameter (width_tag)
  {
    // {paramWidth, "Width", "width", 50.0f, 150.0f, "%"},
    return float_param ("%", 50., 150., 100., 0.1);
  }
  //----------------------------------------------------------------------------
  struct low_cut_tag {};
  void set (low_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramLowCut, v);
  }

  static constexpr auto get_parameter (low_cut_tag)
  {
    // {paramLowCut, "Low Cut", "low_cut", 0.0f, 200.0f, "Hz"},
    return frequency_parameter_from_zero (200.0, 0.);
  }
  //----------------------------------------------------------------------------
  struct high_cut_tag {};
  void set (high_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramHighCut, v);
  }

  static constexpr auto get_parameter (high_cut_tag)
  {
    // {paramHighCut, "High Cut", "high_cut", 1000.0f, 16000.0f, "Hz"}
    return frequency_parameter (1000.0, 16000.0, 16000.0);
  }
  //----------------------------------------------------------------------------
  using parameters = mp11::
    mp_list<dry_tag, wet_tag, size_tag, width_tag, low_cut_tag, high_cut_tag>;
  //----------------------------------------------------------------------------
  early_reflections()
  {
    mp11::mp_for_each<parameters> ([=] (auto p) {
      set (p, get_parameter (p).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { _dsp.reset (pc.get_sample_rate()); }
  //----------------------------------------------------------------------------
  void process_block_replacing (std::array<float*, 2> chnls, uint samples)
  {
    _dsp.process (chnls.data(), samples);
  }
  //----------------------------------------------------------------------------
private:
  // from dragonfly's sources
  enum Parameters {
    paramDry = 0,
    paramWet,
    paramProgram,
    paramSize,
    paramWidth,
    paramLowCut,
    paramHighCut,
    paramCount
  };

  er_compile_firewall _dsp;
};

}} // namespace artv::dragonfly
