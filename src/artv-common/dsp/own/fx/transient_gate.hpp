#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/transient_gate.hpp"

namespace artv {

//------------------------------------------------------------------------------
struct transient_gate_fx {
public:
  using V = vec<double, 1>;
  using T = double;

  static constexpr dsp_types dsp_type  = dsp_types::dynamics;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct output_tag {};

  void set (output_tag, uint v) { _is_tonal = (v == 1); }

  static constexpr auto get_parameter (output_tag)
  {
    return choice_param (0, make_cstr_array ("Transient", "Tonal"), 8);
  }
  //----------------------------------------------------------------------------
  struct detector_hipass_tag {};

  void set (detector_hipass_tag, float v)
  {
    v = midi_note_to_hz (v);
    if (v == _hipass) {
      return;
    }
    _hipass = v;
    _transient.set_detector_hipass (v);
  }

  static constexpr auto get_parameter (detector_hipass_tag)
  {
    return frequency_parameter (10.0, 1000.0, 75.0);
  }
  //----------------------------------------------------------------------------
  struct detector_recovery_tag {};

  void set (detector_recovery_tag, float v)
  {
    v *= 0.01f;
    if (v == _recovery) {
      return;
    }
    _recovery = v;
    _transient.set_detector_recovery (v);
  }

  static constexpr auto get_parameter (detector_recovery_tag)
  {
    return float_param ("%", 0.0, 100.0, 35., 0.1);
  }
  //----------------------------------------------------------------------------
  struct detector_channels_tag {};

  void set (detector_channels_tag, uint v)
  {
    _transient.set_detector_channels (v);
  }

  static constexpr auto get_parameter (detector_channels_tag)
  {
    return choice_param (0, make_cstr_array ("Stereo", "Mid", "L", "R"), 16);
  }
  //----------------------------------------------------------------------------
  struct detector_shape_tag {};

  void set (detector_shape_tag, uint v) { _transient.set_detector_shape (v); }

  static constexpr auto get_parameter (detector_shape_tag)
  {
    return choice_param (
      0,
      make_cstr_array ("Linear", "Squared", "Cubic", "Quadratic", "Quintic"),
      16);
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};

  void set (decay_tag, float v)
  {
    v *= 0.01f;
    if (v == _decay) {
      return;
    }
    _decay = v;
    _transient.set_curve_decay (v);
  }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("%", 0.0, 100.0, 25., 0.1);
  }
  //----------------------------------------------------------------------------
  struct decay_shape_tag {};

  void set (decay_shape_tag, uint v) { _transient.set_curve_decay_shape (v); }

  static constexpr auto get_parameter (decay_shape_tag)
  {
    return get_parameter (detector_shape_tag {});
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    output_tag,
    detector_hipass_tag,
    detector_recovery_tag,
    detector_channels_tag,
    detector_shape_tag,
    decay_tag,
    decay_shape_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _recovery = -1.;
    _decay    = -1.;
    _hipass   = -1.;
    _is_tonal = false;
    _transient.reset (pc.get_sample_rate());
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    static constexpr uint in  = 0;
    static constexpr uint out = 1;

    std::array<std::array<double_x2, 2>, 32> io;

    size_t done = 0;
    while (done < samples) {
      uint blocksize = std::min (samples - done, io.size());
      // interleaving
      for (uint i = 0; i < blocksize; ++i) {
        io[i][in][0] = ins[0][done + i];
        io[i][in][1] = ins[1][done + i];
      }
      // process
      for (uint i = 0; i < blocksize; ++i) {
        io[i][out] = _transient.tick (io[i][in]);
      }
      // deinterleaving
      if (!_is_tonal) {
        for (uint i = 0; i < blocksize; ++i) {
          outs[0][done + i] = io[i][out][0];
          outs[1][done + i] = io[i][out][1];
        }
      }
      else {
        for (uint i = 0; i < blocksize; ++i) {
          auto diff         = io[i][in] - io[i][out];
          outs[0][done + i] = diff[0];
          outs[1][done + i] = diff[1];
        }
      }
      done += blocksize;
    };
  }
  //----------------------------------------------------------------------------
private:
  transient_gate _transient;
  float          _hipass;
  float          _decay;
  float          _recovery;
  bool           _is_tonal;
};
//------------------------------------------------------------------------------
} // namespace artv
