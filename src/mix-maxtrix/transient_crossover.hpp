#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/transient_gate.hpp"

namespace artv {

//------------------------------------------------------------------------------
struct mixmaxtrix_transient_crossover {
public:
  using V = vec<double, 1>;
  using T = double;

  static constexpr dsp_types dsp_type  = dsp_types::dynamics;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 2;
  //----------------------------------------------------------------------------
  struct transient_output {};

  void set (transient_output, int)
  {
    // dummy, to be used by processor
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_parameter (transient_output)
  {
    constexpr auto str_array = make_cstr_array (
      "bus2", "bus3", "bus4", "bus5", "bus6", "bus7", "bus8");
    constexpr uint n_future_choices = 16;

    return choice_param (1, str_array, n_future_choices);
  }
  //----------------------------------------------------------------------------
  struct detector_threshold_tag {};

  void set (detector_threshold_tag, float v)
  {
    _transient.set_detector_threshold (v);
  }

  static constexpr auto get_parameter (detector_threshold_tag)
  {
    return float_param ("dB", 7.0, 30.0, 12.0, 0.01);
  }
  //----------------------------------------------------------------------------
  struct detector_curve_tag {};

  void set (detector_curve_tag, uint v) { _transient.set_detector_curve (v); }

  static constexpr auto get_parameter (detector_curve_tag)
  {
    return choice_param (
      2,
      make_cstr_array (
        "Inv Cubic", "Inv Squared", "Linear", "Squared", "Cubic"),
      16);
  }
  //----------------------------------------------------------------------------
  struct detector_channel_mode_tag {};

  void set (detector_channel_mode_tag, uint v)
  {
    _transient.set_channel_mode (v);
  }

  static constexpr auto get_parameter (detector_channel_mode_tag)
  {
    return choice_param (0, make_cstr_array ("Stereo", "Mid", "L", "R"), 16);
  }
  //----------------------------------------------------------------------------
  struct detector_hysteresis_tag {};

  void set (detector_hysteresis_tag, float v)
  {
    _transient.set_detector_hysteresis (v);
  }

  static constexpr auto get_parameter (detector_hysteresis_tag)
  {
    return float_param ("", 0.0, 23.0, 7.0, 0.01);
  }
  //----------------------------------------------------------------------------
  struct detector_speed_tag {};

  void set (detector_speed_tag, float v)
  {
    if (v == _speed) {
      return;
    }
    _speed = v;
    _transient.set_detector_speed (v);
  }

  static constexpr auto get_parameter (detector_speed_tag)
  {
    return float_param ("", 0.0, 1.0, 0.35, 0.001);
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};

  void set (decay_tag, float v)
  {
    if (v == _decay) {
      return;
    }
    _decay = v;
    _transient.set_decay (v);
  }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("", 0.0, 1.0, 0.5, 0.001);
  }
  //----------------------------------------------------------------------------
  struct decay_curve_tag {};

  void set (decay_curve_tag, uint v) { _transient.set_decay_curve (v); }

  static constexpr auto get_parameter (decay_curve_tag)
  {
    return get_parameter (detector_curve_tag {});
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    transient_output,
    detector_threshold_tag,
    detector_hysteresis_tag,
    detector_curve_tag,
    detector_channel_mode_tag,
    detector_speed_tag,
    decay_tag,
    decay_curve_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _speed = -1.;
    _decay = -1.;
    _transient.reset (pc.get_sample_rate());
    mp11::mp_for_each<parameters> ([&] (auto type) {
      set (type, get_parameter (type).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    std::array<double_x2, 32>                in;
    std::array<std::array<double_x2, 2>, 32> out;

    size_t done = 0;
    while (done < samples) {
      uint blocksize = std::min (samples - done, in.size());
      // interleaving
      for (uint i = 0; i < blocksize; ++i) {
        in[i][0] = ins[0][done + i];
        in[i][1] = ins[1][done + i];
      }
      // process
      for (uint i = 0; i < blocksize; ++i) {
        auto spl       = in[i];
        auto transient = _transient.tick (spl);
        out[i][0]      = spl - transient;
        out[i][1]      = transient;
      }
      // deinterleaving
      for (uint i = 0; i < blocksize; ++i) {
        outs[0][done + i] = out[i][0][0];
        outs[1][done + i] = out[i][0][1];
        outs[2][done + i] = out[i][1][0];
        outs[3][done + i] = out[i][1][1];
      }
      done += blocksize;
    };
  }
  //----------------------------------------------------------------------------
private:
  transient_gate _transient;
  float          _speed;
  float          _decay;
};
//------------------------------------------------------------------------------
} // namespace artv
