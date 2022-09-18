#pragma once

// Airwindows console7 Port.
//
// Original code left mostly intact, just run through clang-format and adapted
// to the interface of this project.
//
// https://www.patreon.com/airwindows
// https://github.com/airwindows/airwindows
// from commit 7bb7fde13a9c94af242835823b813e6bcd0f20c8

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/stereo_summing_processor.hpp"
#include "artv-common/dsp/third_party/airwindows/common.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace airwindows {

class console7bus {
public:
  // DSP------------------------------------------------------------------------
  console7bus() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    A           = 1.0;
    gainchase   = -1.0;
    chasespeed  = 64.0;
    for (int x = 0; x < 15; x++) {
      biquadA[x] = 0.0;
      biquadB[x] = 0.0;
    }
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpdL = 1.0;
    while (fpdL < 16386)
      fpdL = rand() * UINT32_MAX;
    fpdR = 1.0;
    while (fpdR < 16386)
      fpdR = rand() * UINT32_MAX;
#endif
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

    long double inputgain = A * 1.03;

    if (gainchase != inputgain)
      chasespeed *= 2.0;
    if (chasespeed > samples)
      chasespeed = samples;
    if (gainchase < 0.0)
      gainchase = inputgain;

    biquadB[0] = biquadA[0] = 20000.0 / sample_rate;
    biquadA[1]              = 0.618033988749894848204586;
    biquadB[1]              = 0.5;

    double K    = tan (M_PI * biquadA[0]); // lowpass
    double norm = 1.0 / (1.0 + K / biquadA[1] + K * K);
    biquadA[2]  = K * K * norm;
    biquadA[3]  = 2.0 * biquadA[2];
    biquadA[4]  = biquadA[2];
    biquadA[5]  = 2.0 * (K * K - 1.0) * norm;
    biquadA[6]  = (1.0 - K / biquadA[1] + K * K) * norm;

    K          = tan (M_PI * biquadA[0]);
    norm       = 1.0 / (1.0 + K / biquadB[1] + K * K);
    biquadB[2] = K * K * norm;
    biquadB[3] = 2.0 * biquadB[2];
    biquadB[4] = biquadB[2];
    biquadB[5] = 2.0 * (K * K - 1.0) * norm;
    biquadB[6] = (1.0 - K / biquadB[1] + K * K) * norm;

    while (--samples < ((uint) -1ull)) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER
      if (abs (inputSampleL) < 1.18e-43)
        inputSampleL = fpdL * 1.18e-43;
      if (abs (inputSampleR) < 1.18e-43)
        inputSampleR = fpdR * 1.18e-43;
#endif
      long double outSampleL = biquadA[2] * inputSampleL
        + biquadA[3] * biquadA[7] + biquadA[4] * biquadA[8]
        - biquadA[5] * biquadA[9] - biquadA[6] * biquadA[10];
      biquadA[8]   = biquadA[7];
      biquadA[7]   = inputSampleL;
      inputSampleL = outSampleL;
      biquadA[10]  = biquadA[9];
      biquadA[9]   = inputSampleL; // DF1 left

      long double outSampleR = biquadA[2] * inputSampleR
        + biquadA[3] * biquadA[11] + biquadA[4] * biquadA[12]
        - biquadA[5] * biquadA[13] - biquadA[6] * biquadA[14];
      biquadA[12]  = biquadA[11];
      biquadA[11]  = inputSampleR;
      inputSampleR = outSampleR;
      biquadA[14]  = biquadA[13];
      biquadA[13]  = inputSampleR; // DF1 right

      chasespeed *= 0.9999;
      chasespeed -= 0.01;
      if (chasespeed < 64.0)
        chasespeed = 64.0;
      // we have our chase speed compensated for recent fader activity
      gainchase = (((gainchase * chasespeed) + inputgain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor
      if (1.0 != gainchase) {
        inputSampleL *= sqrt (gainchase);
        inputSampleR *= sqrt (gainchase);
      }
      // done with trim control

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      inputSampleL = ((asin (inputSampleL * abs (inputSampleL))
                       / ((abs (inputSampleL) == 0.0) ? 1 : abs (inputSampleL)))
                      * 0.618033988749894848204586)
        + (asin (inputSampleL) * 0.381966011250105);
      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      inputSampleR = ((asin (inputSampleR * abs (inputSampleR))
                       / ((abs (inputSampleR) == 0.0) ? 1 : abs (inputSampleR)))
                      * 0.618033988749894848204586)
        + (asin (inputSampleR) * 0.381966011250105);
      // this is an asin version of Spiral blended with regular asin
      // ConsoleBuss. It's blending between two different harmonics in the
      // overtones of the algorithm

      outSampleL = biquadB[2] * inputSampleL + biquadB[3] * biquadB[7]
        + biquadB[4] * biquadB[8] - biquadB[5] * biquadB[9]
        - biquadB[6] * biquadB[10];
      biquadB[8]   = biquadB[7];
      biquadB[7]   = inputSampleL;
      inputSampleL = outSampleL;
      biquadB[10]  = biquadB[9];
      biquadB[9]   = inputSampleL; // DF1 left

      outSampleR = biquadB[2] * inputSampleR + biquadB[3] * biquadB[11]
        + biquadB[4] * biquadB[12] - biquadB[5] * biquadB[13]
        - biquadB[6] * biquadB[14];
      biquadB[12]  = biquadB[11];
      biquadB[11]  = inputSampleR;
      inputSampleR = outSampleR;
      biquadB[14]  = biquadB[13];
      biquadB[13]  = inputSampleR; // DF1 right

      if (1.0 != gainchase) {
        inputSampleL *= sqrt (gainchase);
        inputSampleR *= sqrt (gainchase);
      }
      // we re-amplify after the distortion relative to how much we cut back
      // previously.
#if AIRWINDOWS_FP_DITHER_ENABLE
      // begin 64 bit stereo floating point dither
      int expon;
      frexp ((double) inputSampleL, &expon);
      fpdL ^= fpdL << 13;
      fpdL ^= fpdL >> 17;
      fpdL ^= fpdL << 5;
      inputSampleL
        += ((double (fpdL) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      frexp ((double) inputSampleR, &expon);
      fpdR ^= fpdR << 13;
      fpdR ^= fpdR >> 17;
      fpdR ^= fpdR << 5;
      inputSampleR
        += ((double (fpdR) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      // end 64 bit stereo floating point dither
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
  void set (drive_tag, float value)
  {
    A = db_to_gain (-value);
  }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag>;
  //----------------------------------------------------------------------------
private:
  double gainchase;
  double chasespeed;

  long double biquadA[15];
  long double biquadB[15];

#if AIRWINDOWS_FP_DITHER_ENABLE
  uint32_t fpdL;
  uint32_t fpdR;
#endif
  float A; // gain
  float sample_rate;
};
//------------------------------------------------------------------------------
class console7channel : public stereo_summing_processor {
public:
  // DSP------------------------------------------------------------------------
  console7channel() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    A           = 1.0;
    gainchase   = -1.0;
    chasespeed  = 64.0;
    for (int x = 0; x < 15; x++) {
      biquadA[x] = 0.0;
    }
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpdL = 1.0;
    while (fpdL < 16386)
      fpdL = rand() * UINT32_MAX;
    fpdR = 1.0;
    while (fpdR < 16386)
      fpdR = rand() * UINT32_MAX;
#endif
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_adding (
    xspan<T*>       outs,
    xspan<T const*> ins,
    uint            samples,
    bool            dst_sum)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

    long double inputgain = A * 1.272019649514069;
    // which is, in fact, the square root of 1.618033988749894848204586...
    // this happens to give us a boost factor where the track continues to get
    // louder even as it saturates and loses a bit of peak energy.
    // Console7Channel channels go to 12! (.272,etc) Neutral gain through the
    // whole system with a full scale sine ia 0.772 on the gain knob
    if (gainchase != inputgain)
      chasespeed *= 2.0;
    if (chasespeed > samples)
      chasespeed = samples;
    if (gainchase < 0.0)
      gainchase = inputgain;

    biquadA[0] = 20000.0 / sample_rate;
    biquadA[1] = 1.618033988749894848204586;

    double K    = tan (M_PI * biquadA[0]); // lowpass
    double norm = 1.0 / (1.0 + K / biquadA[1] + K * K);
    biquadA[2]  = K * K * norm;
    biquadA[3]  = 2.0 * biquadA[2];
    biquadA[4]  = biquadA[2];
    biquadA[5]  = 2.0 * (K * K - 1.0) * norm;
    biquadA[6]  = (1.0 - K / biquadA[1] + K * K) * norm;

    while (--samples < ((uint) -1ull)) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER
      if (abs (inputSampleL) < 1.18e-43)
        inputSampleL = fpdL * 1.18e-43;
      if (abs (inputSampleR) < 1.18e-43)
        inputSampleR = fpdR * 1.18e-43;
#endif
      long double outSampleL = biquadA[2] * inputSampleL
        + biquadA[3] * biquadA[7] + biquadA[4] * biquadA[8]
        - biquadA[5] * biquadA[9] - biquadA[6] * biquadA[10];
      biquadA[8]   = biquadA[7];
      biquadA[7]   = inputSampleL;
      inputSampleL = outSampleL;
      biquadA[10]  = biquadA[9];
      biquadA[9]   = inputSampleL; // DF1 left

      long double outSampleR = biquadA[2] * inputSampleR
        + biquadA[3] * biquadA[11] + biquadA[4] * biquadA[12]
        - biquadA[5] * biquadA[13] - biquadA[6] * biquadA[14];
      biquadA[12]  = biquadA[11];
      biquadA[11]  = inputSampleR;
      inputSampleR = outSampleR;
      biquadA[14]  = biquadA[13];
      biquadA[13]  = inputSampleR; // DF1 right

      chasespeed *= 0.9999;
      chasespeed -= 0.01;
      if (chasespeed < 64.0)
        chasespeed = 64.0;
      // we have our chase speed compensated for recent fader activity
      gainchase = (((gainchase * chasespeed) + inputgain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor
      if (1.0 != gainchase) {
        inputSampleL *= pow (gainchase, 3);
        inputSampleR *= pow (gainchase, 3);
      }
      // this trim control cuts back extra hard because we will amplify after
      // the distortion that will shift the distortion/antidistortion curve, in
      // order to make faded settings slightly 'expanded' and fall back in the
      // soundstage, subtly

      if (inputSampleL > 1.097)
        inputSampleL = 1.097;
      if (inputSampleL < -1.097)
        inputSampleL = -1.097;
      inputSampleL = ((sin (inputSampleL * abs (inputSampleL))
                       / ((abs (inputSampleL) == 0.0) ? 1 : abs (inputSampleL)))
                      * 0.8)
        + (sin (inputSampleL) * 0.2);
      if (inputSampleR > 1.097)
        inputSampleR = 1.097;
      if (inputSampleR < -1.097)
        inputSampleR = -1.097;
      inputSampleR = ((sin (inputSampleR * abs (inputSampleR))
                       / ((abs (inputSampleR) == 0.0) ? 1 : abs (inputSampleR)))
                      * 0.8)
        + (sin (inputSampleR) * 0.2);
      // this is a version of Spiral blended 80/20 with regular Density
      // ConsoleChannel. It's blending between two different harmonics in the
      // overtones of the algorithm

      if (1.0 != gainchase && 0.0 != gainchase) {
        inputSampleL /= gainchase;
        inputSampleR /= gainchase;
      }
      // we re-amplify after the distortion relative to how much we cut back
      // previously.

#if AIRWINDOWS_FP_DITHER_ENABLE
      // begin 64 bit stereo floating point dither
      int expon;
      frexp ((double) inputSampleL, &expon);
      fpdL ^= fpdL << 13;
      fpdL ^= fpdL >> 17;
      fpdL ^= fpdL << 5;
      inputSampleL
        += ((double (fpdL) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      frexp ((double) inputSampleR, &expon);
      fpdR ^= fpdR << 13;
      fpdR ^= fpdR >> 17;
      fpdR ^= fpdR << 5;
      inputSampleR
        += ((double (fpdR) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      // end 64 bit stereo floating point dither
#endif
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
    xspan<float*>       outs,
    xspan<float const*> ins,
    uint                samples,
    bool                dst_sum) override
  {
    process_block_adding (outs, ins, samples, dst_sum);
  }
  //----------------------------------------------------------------------------
  void process (
    xspan<double*>       outs,
    xspan<double const*> ins,
    uint                 samples,
    bool                 dst_sum) override
  {
    process_block_adding (outs, ins, samples, dst_sum);
  }
  // Parameters (call once per block) ------------------------------------------
  void set (console7bus::drive_tag, float value)
  {
    A = db_to_gain (value) * 0.772;
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<>;
  //----------------------------------------------------------------------------
private:
  double gainchase;
  double chasespeed;

  long double biquadA[15];

#if AIRWINDOWS_FP_DITHER_ENABLE
  uint32_t fpdL;
  uint32_t fpdR;
#endif

  float A; // gain
  float sample_rate;
};
}} // namespace artv::airwindows
