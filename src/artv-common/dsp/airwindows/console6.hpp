#pragma once

// Airwindows Console6 Port.
//
// Original code left mostly intact, just run through clang-format and adapted
// to the interface of this project.
//
// https://www.patreon.com/airwindows
// https://github.com/airwindows/airwindows
// from commit 7bb7fde13a9c94af242835823b813e6bcd0f20c8

#include "artv-common/dsp/airwindows/common.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/stereo_summing_processor.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv { namespace airwindows {

class console6bus {
public:
  // DSP------------------------------------------------------------------------
  console6bus() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { A = 1.0; }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

    double gain = A;

    while (--samples < ((uint) -1ull)) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0].clamp (inputSampleL, T {});
      inputSampleR = dither[1].clamp (inputSampleR, T {});
#endif
      if (gain != 1.0) {
        inputSampleL *= gain;
        inputSampleR *= gain;
      }

      // TODO faster "pow" for -1 to 1 ranges?

      // encode/decode courtesy of torridgristle under the MIT license
      // Inverse Square 1-(1-x)^2 and 1-(1-x)^0.5

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      else if (inputSampleL > 0.0)
        inputSampleL = 1.0 - pow (1.0 - inputSampleL, 0.5);

      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      else if (inputSampleL < 0.0)
        inputSampleL = -1.0 + pow (1.0 + inputSampleL, 0.5);

      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      else if (inputSampleR > 0.0)
        inputSampleR = 1.0 - pow (1.0 - inputSampleR, 0.5);

      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      else if (inputSampleR < 0.0)
        inputSampleR = -1.0 + pow (1.0 + inputSampleR, 0.5);

#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0](inputSampleL, T {});
      inputSampleR = dither[1](inputSampleR, T {});
#endif
      *out1 = inputSampleL;
      *out2 = inputSampleR;

      in1++;
      in2++;
      out1++;
      out2++;
    }
  }
  // Parameters (call once per block) ------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float value) { A = db_to_gain (-value); }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag>;
  //----------------------------------------------------------------------------
private:
#if AIRWINDOWS_FP_DITHER_ENABLE
  get_dither_for<console6bus>::type dither[2];
#endif
  float A; // gain
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
class console6channel : public stereo_summing_processor {
public:
  // DSP------------------------------------------------------------------------
  console6channel() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) { A = 1.0; }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_adding (
    crange<T*>       outs,
    crange<T const*> ins,
    uint             samples,
    bool             dst_sum)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

    double gain = A;

    while (--samples < ((uint) -1ull)) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0].clamp (inputSampleL, T {});
      inputSampleR = dither[1].clamp (inputSampleR, T {});
#endif
      // changed: this optimization is a pesimization
      // if (gain != 1.0) {
      inputSampleL *= gain;
      inputSampleR *= gain;
      //}

      // TODO faster "pow" for -1 to 1 ranges?

      // encode/decode courtesy of torridgristle under the MIT license
      // Inverse Square 1-(1-x)^2 and 1-(1-x)^0.5

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      else if (inputSampleL > 0.0)
        inputSampleL = 1.0 - pow (1.0 - inputSampleL, 2.0);

      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      else if (inputSampleL < 0.0)
        inputSampleL = -1.0 + pow (1.0 + inputSampleL, 2.0);

      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      else if (inputSampleR > 0.0)
        inputSampleR = 1.0 - pow (1.0 - inputSampleR, 2.0);

      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      else if (inputSampleR < 0.0)
        inputSampleR = -1.0 + pow (1.0 + inputSampleR, 2.0);

#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0](inputSampleL, T {});
      inputSampleR = dither[1](inputSampleR, T {});
#endif
      // changed while porting
      if (dst_sum) {
        *out1 += inputSampleL;
        *out2 += inputSampleR;
      }
      else {
        *out1 = inputSampleL;
        *out2 = inputSampleR;
      }

      in1++;
      in2++;
      out1++;
      out2++;
    }
  }
  //----------------------------------------------------------------------------
  void process (
    crange<float*>       outs,
    crange<float const*> ins,
    uint                 samples,
    bool                 dst_sum) override
  {
    process_block_adding (outs, ins, samples, dst_sum);
  }
  //----------------------------------------------------------------------------
  void process (
    crange<double*>       outs,
    crange<double const*> ins,
    uint                  samples,
    bool                  dst_sum) override
  {
    process_block_adding (outs, ins, samples, dst_sum);
  }
  // Parameters (call once per block) ------------------------------------------
  void set (console6bus::drive_tag, float value) { A = db_to_gain (value); }
  //----------------------------------------------------------------------------
  using parameters = mp_list<>;
  //----------------------------------------------------------------------------
private:
#if AIRWINDOWS_FP_DITHER_ENABLE
  get_dither_for<console6channel>::type dither[2];
#endif
  float A; // gain
};
}} // namespace artv::airwindows
