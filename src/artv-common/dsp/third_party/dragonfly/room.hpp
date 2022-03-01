#pragma once

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

class room {
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
  struct early_tag {};
  void set (early_tag, float v)
  {
    _dsp.set_parameter (paramEarly, db_to_gain (v) * 100.f);
  }

  static constexpr auto get_parameter (early_tag)
  {
    // {paramEarly,      "Early Level", "early_level",  0.0f,   100.0f,   "%"},
    return float_param ("dB", -60., 0., -35., 0.1);
  }
  //----------------------------------------------------------------------------
  struct early_send_tag {};
  void set (early_send_tag, float v)
  {
    _dsp.set_parameter (paramEarlySend, db_to_gain (v) * 100.f);
  }

  static constexpr auto get_parameter (early_send_tag)
  {
    // {paramEarlySend,  "Early Send",  "early_send",   0.0f,   100.0f,   "%"},
    return float_param ("dB", -60., 0., -35., 0.1);
  }
  //----------------------------------------------------------------------------
  struct late_tag {};
  void set (late_tag, float v)
  {
    _dsp.set_parameter (paramLate, db_to_gain (v) * 100.f);
  }

  static constexpr auto get_parameter (late_tag)
  {
    // {paramLate,       "Late Level",  "late_level",   0.0f,   100.0f,   "%"},
    return float_param ("dB", -60., 0., -15., 0.1);
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
  struct predelay_tag {};
  void set (predelay_tag, float v) { _dsp.set_parameter (paramPredelay, v); }

  static constexpr auto get_parameter (predelay_tag)
  {
    // {paramPredelay, "Predelay", "predelay", 0.0f, 100.0f, "ms"},
    return float_param ("ms", 0., 100., 10., 0.1, 0.8);
  }
  //----------------------------------------------------------------------------
  struct size_tag {};
  void set (size_tag, float v) { _dsp.set_parameter (paramSize, v); }

  static constexpr auto get_parameter (size_tag)
  {
    // {paramSize, "Size", "size", 10.0f, 60.0f, "m"},
    return float_param ("m", 8., 32., 12., 0.1);
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
  struct diffuse_tag {};
  void set (diffuse_tag, float v) { _dsp.set_parameter (paramDiffuse, v); }

  static constexpr auto get_parameter (diffuse_tag)
  {
    // {paramDiffuse,    "Diffuse",     "diffuse",      0.0f,   100.0f,   "%"},
    return float_param ("%", 0, 100., 20., 0.1);
  }
  //----------------------------------------------------------------------------
  struct low_cut_tag {};
  void set (low_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramInLowCut, v);
  }

  static constexpr auto get_parameter (low_cut_tag)
  {
    // {paramLowCut,     "Low Cut",     "low_cut",      0.0f,   200.0f,  "Hz"},
    return frequency_parameter_from_zero (200.0, 60.0);
  }
  //----------------------------------------------------------------------------
  struct low_boost_freq_tag {};
  void set (low_boost_freq_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramBoostLPF, v);
  }

  static constexpr auto get_parameter (low_boost_freq_tag)
  {
    // {paramBoostLPF, "Boost Freq", "boost_freq", 50.0f, 1050.0f, "Hz"},
    return frequency_parameter (50.0, 1050.0, 440.0);
  }
  //----------------------------------------------------------------------------
  struct low_boost_tag {};
  void set (low_boost_tag, float v) { _dsp.set_parameter (paramBoost, v); }

  static constexpr auto get_parameter (low_boost_tag)
  {
    // {paramBoost, "Low Boost", "low_boost", 0.0f, 100.0f, "%"},
    return float_param ("%", 0., 100., 5., 0.01);
  }
  //----------------------------------------------------------------------------
  struct high_cut_tag {};
  void set (high_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramInHighCut, v);
  }

  static constexpr auto get_parameter (high_cut_tag)
  {
    // {paramInHighCut, "High Cut", "in_high_cut", 1000.0f, 16000.0f, "Hz"},
    return frequency_parameter (1000.0, 16000.0, 16000.0);
  }
  //----------------------------------------------------------------------------
  struct early_damp_tag {};
  void set (early_damp_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramEarlyDamp, v);
  }

  static constexpr auto get_parameter (early_damp_tag)
  {
    // {paramEarlyDamp, "Early Damp", "early_damp", 1000.0f, 16000.0f, "Hz"}
    return frequency_parameter (1000.0, 16000.0, 3000.0);
  }
  //----------------------------------------------------------------------------
  struct late_damp_tag {};
  void set (late_damp_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramLateDamp, v);
  }

  static constexpr auto get_parameter (late_damp_tag)
  {
    // {paramLateDamp, "Late Damp", "late_damp", 1000.0f, 16000.0f, "Hz"},
    return frequency_parameter (1000.0, 16000.0, 5000.0);
  }
  //----------------------------------------------------------------------------
  struct spin_tag {};
  void set (spin_tag, float v)
  {
    v = midi_note_to_hz (v);
    _dsp.set_parameter (paramSpin, v);
  }

  static constexpr auto get_parameter (spin_tag)
  {
    // {paramSpin, "Spin", "spin", 0.0f, 10.0f, "Hz"},
    return frequency_parameter_from_zero (5.0, 1.5);
  }
  //----------------------------------------------------------------------------
  struct wander_tag {};
  void set (wander_tag, float v) { _dsp.set_parameter (paramWander, v); }

  static constexpr auto get_parameter (wander_tag)
  {
    // {paramWander, "Wander", "wander", 0.0f, 40.0f, "ms"},
    return float_param ("%", 0., 100., 20., 0.01);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    dry_tag,
    early_tag,
    early_send_tag,
    late_tag,
    size_tag,
    width_tag,
    predelay_tag,
    decay_tag,
    diffuse_tag,
    low_cut_tag,
    high_cut_tag,
    low_boost_freq_tag,
    late_damp_tag,
    low_boost_tag,
    early_damp_tag,
    spin_tag,
    wander_tag>;
  //----------------------------------------------------------------------------
  room()
  {
    mp11::mp_for_each<parameters> ([=] (auto p) {
      set (p, get_parameter (p).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { _dsp.reset (pc.get_sample_rate()); }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    _dsp.process (outs.data(), ins.data(), samples);
  }
  //----------------------------------------------------------------------------
private:
  // from dragonfly's sources
  enum Parameters {
    paramDry = 0,
    paramEarly,
    paramEarlySend,
    paramLate,
    paramSize,
    paramWidth,
    paramPredelay,
    paramDecay,
    paramDiffuse,
    paramSpin,
    paramWander,
    paramInHighCut,
    paramEarlyDamp,
    paramLateDamp,
    paramBoost,
    paramBoostLPF,
    paramInLowCut,
    paramCount
  };

  room_compile_firewall _dsp;
};

}} // namespace artv::dragonfly
