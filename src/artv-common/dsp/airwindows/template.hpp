#pragma once

// Airwindows "awtemplate" Port.
//
// Original code left mostly intact, just run through clang-format and adapted
// to the interface of this project.
//
// https://www.patreon.com/airwindows
// https://github.com/airwindows/airwindows
// from commit 7bb7fde13a9c94af242835823b813e6bcd0f20c8

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/airwindows/common.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv { namespace airwindows {

class awtemplate {
public:
  // DSP------------------------------------------------------------------------
  awtemplate() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { sample_rate = pc.get_sample_rate(); }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

#if AIRWINDOWS_FP_DITHER_ENABLE
#endif
#if AIRWINDOWS_DENORMALIZATION_ENABLE 0
#endif
#if AIRWINDOWS_SUBSONIC_FILTER_ENABLE 0
#endif
  }
  // Parameters (call once per block) ------------------------------------------
  struct type {};
  void                  set (type, int v) {}
  static constexpr auto get_type_parameter()
  {
    return choice_param (0, make_cstr_array ("v1", "v2 "));
  }
  //----------------------------------------------------------------------------
  struct drive {};
  void set (drive, float v)
  {
    // try to keep more or less the same perceived loudness
  }
  static constexpr auto get_drive_parameter()
  {
    return float_param ("", 0.f, 200.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<type, drive>;
  //----------------------------------------------------------------------------
private:
  // Addon
  float sample_rate;
};

}} // namespace artv::airwindows
