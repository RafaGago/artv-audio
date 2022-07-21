#pragma once

#include "artv-common/dsp/own/classes/add_ducker.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/stereo_processor.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "artv-common/dsp/third_party/dragonfly/compilation_firewalls.hpp"

namespace artv { namespace dragonfly {

class plate : private add_ducker<double_x2> {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
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
    return float_param ("dB", -60., 0., -18., 0.1);
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
    return frequency_parameter_from_zero (200.0);
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
  struct algorithm_tag {};
  void set (algorithm_tag, int v) { _dsp.set_parameter (paramAlgorithm, v); }

  static constexpr auto get_parameter (algorithm_tag)
  {
    // {paramAlgorithm, "Algorithm", "algorithm", 0.0f, 2.0f, ""},
    return choice_param (0, make_cstr_array ("Simple", "Nested", "Tank"));
  }
  //----------------------------------------------------------------------------
  struct predelay_tag {};
  void set (predelay_tag, float v) { _dsp.set_parameter (paramPredelay, v); }

  static constexpr auto get_parameter (predelay_tag)
  {
    // {paramPredelay, "Predelay", "predelay", 0.0f, 100.0f, "ms"},
    return float_param ("ms", 0., 100., 10., 0.1, 0.8);
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};
  void set (decay_tag, float v) { _dsp.set_parameter (paramDecay, v); }

  static constexpr auto get_parameter (decay_tag)
  {
    // {paramDecay, "Decay", "decay", 0.1f, 10.0f, "s"},
    return float_param ("s", 0.1, 10., 2., 0.1, 0.8);
  }
  //----------------------------------------------------------------------------
  struct early_damp_tag {};
  void set (early_damp_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramDamp, v);
  }

  static constexpr auto get_parameter (early_damp_tag)
  {
    // {paramDamp, "Dampen", "early_damp", 1000.0f, 16000.0f, "Hz"}};
    return frequency_parameter (1000.0, 16000.0, 16000.0);
  }
  //----------------------------------------------------------------------------
  using add_ducker::get_parameter;
  using add_ducker::set;
  using ducking_speed_tag     = add_ducker::ducking_speed_tag;
  using ducking_threshold_tag = add_ducker::ducking_threshold_tag;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    dry_tag,
    wet_tag,
    width_tag,
    low_cut_tag,
    high_cut_tag,
    algorithm_tag,
    predelay_tag,
    decay_tag,
    early_damp_tag,
    ducking_speed_tag,
    ducking_threshold_tag>;
  //----------------------------------------------------------------------------
  plate()
  {
    mp11::mp_for_each<parameters> ([=] (auto p) {
      set (p, get_parameter (p).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _dsp.reset (pc.get_sample_rate());
    add_ducker::reset (pc.get_sample_rate());
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    add_ducker::process (
      outs,
      ins,
      samples,
      [=] (crange<T*> outs_fw, crange<T const*> ins_fw, uint samples_fw) {
        this->process_intern (outs_fw, ins_fw, samples_fw);
      });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  void process_intern (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    _dsp.process (outs.data(), ins.data(), samples);
  }
  // from dragonfly's sources
  enum Parameters {
    paramDry = 0,
    paramWet,
    paramAlgorithm,
    paramWidth,
    paramPredelay,
    paramDecay,
    paramLowCut,
    paramHighCut,
    paramDamp,
    paramCount
  };

  plate_compile_firewall _dsp;
};

}} // namespace artv::dragonfly
