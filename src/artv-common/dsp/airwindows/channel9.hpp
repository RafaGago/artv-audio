#pragma once

// Airwindows Channel9 Port.
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

class channel9 {
public:
  // DSP-------------------------------------------------------------------------
  channel9() {}
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::exciter;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    B           = 0.0;
    C           = 1.0;
    for (int x = 0; x < 15; x++) {
      biquadA[x] = 0.0;
      biquadB[x] = 0.0;
    }
    iirSampleLA  = 0.0;
    iirSampleRA  = 0.0;
    iirSampleLB  = 0.0;
    iirSampleRB  = 0.0;
    lastSampleAL = lastSampleBL = lastSampleCL = 0.0;
    lastSampleAR = lastSampleBR = lastSampleCR = 0.0;
    flip                                       = false;
    iirAmount                                  = 0.005832;
    threshold = 0.33362176; // instantiating with Neve values
    cutoff    = 28811.0;
  }
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

    double overallscale = 1.0;
    overallscale /= 44100.0;
    overallscale *= sample_rate;
    double localiirAmount = iirAmount / overallscale;
    double localthreshold = threshold; // we've learned not to try and adjust
                                       // threshold for sample rate
    double density  = B * 2.0; // 0-2
    double phattity = density - 1.0;
    if (density > 1.0)
      density = 1.0; // max out at full wet for Spiral aspect
    if (phattity < 0.0)
      phattity = 0.0; //
    double nonLin = 5.0
      - density; // number is smaller for more intense, larger for more subtle
    biquadB[0] = biquadA[0] = cutoff / sample_rate;
    biquadA[1]              = 1.618033988749894848204586;
    biquadB[1]              = 0.618033988749894848204586;

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

#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0].clamp (inputSampleL, T {});
      inputSampleR = dither[1].clamp (inputSampleR, T {});
#endif

      long double tempSample;

      if (biquadA[0] < 0.49999) {
        tempSample = biquadA[2] * inputSampleL + biquadA[3] * biquadA[7]
          + biquadA[4] * biquadA[8] - biquadA[5] * biquadA[9]
          - biquadA[6] * biquadA[10];
        biquadA[8]   = biquadA[7];
        biquadA[7]   = inputSampleL;
        inputSampleL = tempSample;
        biquadA[10]  = biquadA[9];
        biquadA[9]   = inputSampleL; // DF1 left
        tempSample   = biquadA[2] * inputSampleR + biquadA[3] * biquadA[11]
          + biquadA[4] * biquadA[12] - biquadA[5] * biquadA[13]
          - biquadA[6] * biquadA[14];
        biquadA[12]  = biquadA[11];
        biquadA[11]  = inputSampleR;
        inputSampleR = tempSample;
        biquadA[14]  = biquadA[13];
        biquadA[13]  = inputSampleR; // DF1 right
      }

      double dielectricScaleL = abs (2.0 - ((inputSampleL + nonLin) / nonLin));
      double dielectricScaleR = abs (2.0 - ((inputSampleR + nonLin) / nonLin));

      if (flip) {
        iirSampleLA
          = (iirSampleLA * (1.0 - (localiirAmount * dielectricScaleL)))
          + (inputSampleL * localiirAmount * dielectricScaleL);
        inputSampleL = inputSampleL - iirSampleLA;
        iirSampleRA
          = (iirSampleRA * (1.0 - (localiirAmount * dielectricScaleR)))
          + (inputSampleR * localiirAmount * dielectricScaleR);
        inputSampleR = inputSampleR - iirSampleRA;
      }
      else {
        iirSampleLB
          = (iirSampleLB * (1.0 - (localiirAmount * dielectricScaleL)))
          + (inputSampleL * localiirAmount * dielectricScaleL);
        inputSampleL = inputSampleL - iirSampleLB;
        iirSampleRB
          = (iirSampleRB * (1.0 - (localiirAmount * dielectricScaleR)))
          + (inputSampleR * localiirAmount * dielectricScaleR);
        inputSampleR = inputSampleR - iirSampleRB;
      }
      // highpass section
      long double drySampleL = inputSampleL;
      long double drySampleR = inputSampleR;

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      long double phatSampleL = sin (inputSampleL * 1.57079633);
      inputSampleL *= 1.2533141373155;
      // clip to 1.2533141373155 to reach maximum output, or 1.57079633 for pure
      // sine 'phat' version

      long double distSampleL = sin (inputSampleL * abs (inputSampleL))
        / ((abs (inputSampleL) == 0.0) ? 1 : abs (inputSampleL));

      inputSampleL = distSampleL; // purest form is full Spiral
      if (density < 1.0)
        inputSampleL = (drySampleL * (1 - density))
          + (distSampleL * density); // fade Spiral aspect
      if (phattity > 0.0)
        inputSampleL = (inputSampleL * (1 - phattity))
          + (phatSampleL * phattity); // apply original Density on top

      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      long double phatSampleR = sin (inputSampleR * 1.57079633);
      inputSampleR *= 1.2533141373155;
      // clip to 1.2533141373155 to reach maximum output, or 1.57079633 for pure
      // sine 'phat' version

      long double distSampleR = sin (inputSampleR * abs (inputSampleR))
        / ((abs (inputSampleR) == 0.0) ? 1 : abs (inputSampleR));

      inputSampleR = distSampleR; // purest form is full Spiral
      if (density < 1.0)
        inputSampleR = (drySampleR * (1 - density))
          + (distSampleR * density); // fade Spiral aspect
      if (phattity > 0.0)
        inputSampleR = (inputSampleR * (1 - phattity))
          + (phatSampleR * phattity); // apply original Density on top

      // begin L
      double clamp = (lastSampleBL - lastSampleCL) * 0.381966011250105;
      clamp -= (lastSampleAL - lastSampleBL) * 0.6180339887498948482045;
      clamp += inputSampleL - lastSampleAL; // regular slew clamping added

      lastSampleCL = lastSampleBL;
      lastSampleBL = lastSampleAL;
      lastSampleAL = inputSampleL; // now our output relates off lastSampleB

      if (clamp > localthreshold)
        inputSampleL = lastSampleBL + localthreshold;
      if (-clamp > localthreshold)
        inputSampleL = lastSampleBL - localthreshold;

      lastSampleAL = (lastSampleAL * 0.381966011250105)
        + (inputSampleL
           * 0.6180339887498948482045); // split the difference between raw and
                                        // smoothed for buffer
      // end L

      // begin R
      clamp = (lastSampleBR - lastSampleCR) * 0.381966011250105;
      clamp -= (lastSampleAR - lastSampleBR) * 0.6180339887498948482045;
      clamp += inputSampleR - lastSampleAR; // regular slew clamping added

      lastSampleCR = lastSampleBR;
      lastSampleBR = lastSampleAR;
      lastSampleAR = inputSampleR; // now our output relates off lastSampleB

      if (clamp > localthreshold)
        inputSampleR = lastSampleBR + localthreshold;
      if (-clamp > localthreshold)
        inputSampleR = lastSampleBR - localthreshold;

      lastSampleAR = (lastSampleAR * 0.381966011250105)
        + (inputSampleR
           * 0.6180339887498948482045); // split the difference between raw and
                                        // smoothed for buffer
      // end R

      flip = !flip;

      if (C < 1.0) {
        inputSampleL *= C;
        inputSampleR *= C;
      }

      if (biquadB[0] < 0.49999) {
        tempSample = biquadB[2] * inputSampleL + biquadB[3] * biquadB[7]
          + biquadB[4] * biquadB[8] - biquadB[5] * biquadB[9]
          - biquadB[6] * biquadB[10];
        biquadB[8]   = biquadB[7];
        biquadB[7]   = inputSampleL;
        inputSampleL = tempSample;
        biquadB[10]  = biquadB[9];
        biquadB[9]   = inputSampleL; // DF1 left
        tempSample   = biquadB[2] * inputSampleR + biquadB[3] * biquadB[11]
          + biquadB[4] * biquadB[12] - biquadB[5] * biquadB[13]
          - biquadB[6] * biquadB[14];
        biquadB[12]  = biquadB[11];
        biquadB[11]  = inputSampleR;
        inputSampleR = tempSample;
        biquadB[14]  = biquadB[13];
        biquadB[13]  = inputSampleR; // DF1 right
      }

#if AIRWINDOWS_FP_DITHER_ENABLE
      // changed while porting
      inputSampleL = dither[0](inputSampleL, T {});
      inputSampleR = dither[1](inputSampleR, T {});
#endif
      *out1 = (T) inputSampleL;
      *out2 = (T) inputSampleR;

      in1++;
      in2++;
      out1++;
      out2++;
    }
  }
  // Parameters (call once per block) ------------------------------------------
  struct type_tag {};
  void set (type_tag, int v)
  {
    switch (v) {
    case 0:
      iirAmount = 0.005832;
      threshold = 0.33362176;
      cutoff    = 28811.0;
      break; // Neve
    case 1:
      iirAmount = 0.004096;
      threshold = 0.59969536;
      cutoff    = 27216.0;
      break; // API
    case 2:
      iirAmount = 0.004913;
      threshold = 0.84934656;
      cutoff    = 23011.0;
      break; // SSL
    case 3:
      iirAmount = 0.009216;
      threshold = 0.149;
      cutoff    = 18544.0;
      break; // Teac
    case 4:
      iirAmount = 0.011449;
      threshold = 0.092;
      cutoff    = 19748.0;
      break; // Mackie
    default:
      break; // should not happen
    }
  }
  static constexpr auto get_parameter (type_tag)
  {
    return choice_param (
      0, make_cstr_array ("Neve", "API ", "SSL", "Teac ", "Mackie"));
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    B = v / 200;
    C = 1. - (0.4 * B); // try to keep more or less the same perceived loudness
  }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("", 0.f, 200.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = typelist<type_tag, drive_tag>;
  //----------------------------------------------------------------------------
private:
#if AIRWINDOWS_FP_DITHER_ENABLE
  get_dither_for<channel9>::type dither[2];
#endif

  double      iirSampleLA;
  double      iirSampleRA;
  double      iirSampleLB;
  double      iirSampleRB;
  double      lastSampleAL;
  double      lastSampleBL;
  double      lastSampleCL;
  double      lastSampleAR;
  double      lastSampleBR;
  double      lastSampleCR;
  long double biquadA[15];
  long double biquadB[15];
  double      iirAmount;
  double      threshold;
  double      cutoff;
  bool        flip;
  // float       A; // type
  float B; // drive
  float C; // output

  float sample_rate;
};

}} // namespace artv::airwindows
