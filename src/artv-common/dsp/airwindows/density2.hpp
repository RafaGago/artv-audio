#pragma once

// Airwindows "density2" Port.
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
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace airwindows {
class density2 {
public:
  // DSP------------------------------------------------------------------------
  density2() {}
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::distortion;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    A = 0.2;
    B = 0.0;
    D = 1.0;
    gain_compensate();

    ataAL = ataBL = ataCL = lastDiffSampleL = 0.0;
    iirSampleAL = iirSampleBL = last3SampleL = last2SampleL = last1SampleL
      = 0.0;
    ataAR = ataBR = ataCR = lastDiffSampleR = 0.0;
    iirSampleAR = iirSampleBR = last3SampleR = last2SampleR = last1SampleR
      = 0.0;

    sample_rate = pc.get_sample_rate();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, int samples)
  {
    T* in1  = chnls[0];
    T* in2  = chnls[1];
    T* out1 = chnls[0];
    T* out2 = chnls[1];

    double overallscale = 1.0;
    overallscale /= 44100.0;
    overallscale *= sample_rate;
    T density = A;
    T out     = abs (density);
    while (out > 1.0)
      out = out - 1.0;
    density     = density * abs (density);
    T iirAmount = pow (B, 3) / overallscale;
    T output    = C;
    T wet       = D;

    while (--samples >= 0) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER_ENABLE
      if (abs (inputSampleL) < 1.18e-43)
        inputSampleL = fpdL * 1.18e-43;
      if (abs (inputSampleR) < 1.18e-43)
        inputSampleR = fpdR * 1.18e-43;
#endif
      long double drySampleL = inputSampleL;
      long double drySampleR = inputSampleR;

      long double halfwaySampleL
        = (inputSampleL + last1SampleL
           + ((-last2SampleL + last3SampleL) * 0.0414213562373095048801688))
        / 2.0;
      long double halfDrySampleL = halfwaySampleL;

      long double halfwaySampleR
        = (inputSampleR + last1SampleR
           + ((-last2SampleR + last3SampleR) * 0.0414213562373095048801688))
        / 2.0;
      long double halfDrySampleR = halfwaySampleR;

      last3SampleL = last2SampleL;
      last2SampleL = last1SampleL;
      last1SampleL = inputSampleL;
      last3SampleR = last2SampleR;
      last2SampleR = last1SampleR;
      last1SampleR = inputSampleR;

      iirSampleBL
        = (iirSampleBL * (1.0 - iirAmount)) + (halfwaySampleL * iirAmount);
      halfwaySampleL -= iirSampleBL; // highpass section
      iirSampleBR
        = (iirSampleBR * (1.0 - iirAmount)) + (halfwaySampleR * iirAmount);
      halfwaySampleR -= iirSampleBR; // highpass section

      long double bridgerectifier;

      double count = density;
      while (count > 1.0) {
        bridgerectifier = abs (halfwaySampleL) * 1.57079633;
        if (bridgerectifier > 1.57079633)
          bridgerectifier = 1.57079633;
        bridgerectifier = sin (bridgerectifier);
        if (halfwaySampleL > 0.0)
          halfwaySampleL = bridgerectifier;
        else
          halfwaySampleL = -bridgerectifier;
        count = count - 1.0;
      }
      bridgerectifier = abs (halfwaySampleL) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.57079633;
      if (density > 0)
        bridgerectifier = sin (bridgerectifier);
      else
        bridgerectifier = 1
          - cos (bridgerectifier); // produce either boosted or starved version
      if (halfwaySampleL > 0)
        halfwaySampleL
          = (halfwaySampleL * (1.0 - out)) + (bridgerectifier * out);
      else
        halfwaySampleL = (halfwaySampleL * (1.0 - out))
          - (bridgerectifier * out); // blend according to density control

      count = density;
      while (count > 1.0) {
        bridgerectifier = abs (halfwaySampleR) * 1.57079633;
        if (bridgerectifier > 1.57079633)
          bridgerectifier = 1.57079633;
        bridgerectifier = sin (bridgerectifier);
        if (halfwaySampleR > 0.0)
          halfwaySampleR = bridgerectifier;
        else
          halfwaySampleR = -bridgerectifier;
        count = count - 1.0;
      }
      bridgerectifier = abs (halfwaySampleR) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.57079633;
      if (density > 0)
        bridgerectifier = sin (bridgerectifier);
      else
        bridgerectifier = 1
          - cos (bridgerectifier); // produce either boosted or starved version
      if (halfwaySampleR > 0)
        halfwaySampleR
          = (halfwaySampleR * (1.0 - out)) + (bridgerectifier * out);
      else
        halfwaySampleR = (halfwaySampleR * (1.0 - out))
          - (bridgerectifier * out); // blend according to density control

      ataCL = halfwaySampleL - halfDrySampleL;
      ataAL *= 0.915965594177219015;
      ataBL *= 0.915965594177219015;
      ataBL += ataCL;
      ataAL -= ataCL;
      ataCL                       = ataBL;
      long double halfDiffSampleL = ataCL * 0.915965594177219015;

      ataCR = halfwaySampleR - halfDrySampleR;
      ataAR *= 0.915965594177219015;
      ataBR *= 0.915965594177219015;
      ataBR += ataCR;
      ataAR -= ataCR;
      ataCR                       = ataBR;
      long double halfDiffSampleR = ataCR * 0.915965594177219015;

      iirSampleAL
        = (iirSampleAL * (1.0 - iirAmount)) + (inputSampleL * iirAmount);
      inputSampleL -= iirSampleAL; // highpass section
      iirSampleAR
        = (iirSampleAR * (1.0 - iirAmount)) + (inputSampleR * iirAmount);
      inputSampleR -= iirSampleAR; // highpass section

      count = density;
      while (count > 1.0) {
        bridgerectifier = abs (inputSampleL) * 1.57079633;
        if (bridgerectifier > 1.57079633)
          bridgerectifier = 1.57079633;
        bridgerectifier = sin (bridgerectifier);
        if (inputSampleL > 0.0)
          inputSampleL = bridgerectifier;
        else
          inputSampleL = -bridgerectifier;
        count = count - 1.0;
      }
      bridgerectifier = abs (inputSampleL) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.57079633;
      if (density > 0)
        bridgerectifier = sin (bridgerectifier);
      else
        bridgerectifier = 1
          - cos (bridgerectifier); // produce either boosted or starved version
      if (inputSampleL > 0)
        inputSampleL = (inputSampleL * (1 - out)) + (bridgerectifier * out);
      else
        inputSampleL = (inputSampleL * (1 - out))
          - (bridgerectifier * out); // blend according to density control

      count = density;
      while (count > 1.0) {
        bridgerectifier = abs (inputSampleR) * 1.57079633;
        if (bridgerectifier > 1.57079633)
          bridgerectifier = 1.57079633;
        bridgerectifier = sin (bridgerectifier);
        if (inputSampleR > 0.0)
          inputSampleR = bridgerectifier;
        else
          inputSampleR = -bridgerectifier;
        count = count - 1.0;
      }
      bridgerectifier = abs (inputSampleR) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.57079633;
      if (density > 0)
        bridgerectifier = sin (bridgerectifier);
      else
        bridgerectifier = 1
          - cos (bridgerectifier); // produce either boosted or starved version
      if (inputSampleR > 0)
        inputSampleR = (inputSampleR * (1 - out)) + (bridgerectifier * out);
      else
        inputSampleR = (inputSampleR * (1 - out))
          - (bridgerectifier * out); // blend according to density control

      ataCL = inputSampleL - drySampleL;
      ataAL *= 0.915965594177219015;
      ataBL *= 0.915965594177219015;
      ataAL += ataCL;
      ataBL -= ataCL;
      ataCL                   = ataAL;
      long double diffSampleL = ataCL * 0.915965594177219015;

      ataCR = inputSampleR - drySampleR;
      ataAR *= 0.915965594177219015;
      ataBR *= 0.915965594177219015;
      ataAR += ataCR;
      ataBR -= ataCR;
      ataCR                   = ataAR;
      long double diffSampleR = ataCR * 0.915965594177219015;

      inputSampleL = drySampleL
        + ((diffSampleL + halfDiffSampleL + lastDiffSampleL) / 1.187);
      lastDiffSampleL = diffSampleL / 2.0;
      inputSampleL *= output;
      inputSampleL = (drySampleL * (1.0 - wet)) + (inputSampleL * wet);

      inputSampleR = drySampleR
        + ((diffSampleR + halfDiffSampleR + lastDiffSampleR) / 1.187);
      lastDiffSampleR = diffSampleR / 2.0;
      inputSampleR *= output;
      inputSampleR = (drySampleR * (1.0 - wet)) + (inputSampleR * wet);

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
  // Parameters (call once per block)
  // ---------------------------------------------------------------------------
  void gain_compensate()
  {
    float density_gain = 1.;
    if (A > 0) {
      density_gain = db_to_gain (-2. * A);
    }
    else {
      density_gain = db_to_gain ((1. * A) * 6.);
    }
    C = 1. * db_to_gain (6. * B) * density_gain;
  }
  // ---------------------------------------------------------------------------
  struct density_tag {};
  void set (density_tag, float v)
  {
    A = v;
    gain_compensate();
  }
  static constexpr auto get_parameter (density_tag)
  {
    return float_param ("", -1.f, 4.f, 0.f, 0.01f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct highpass_tag {};
  void set (highpass_tag, float v)
  {
    B = v;
    gain_compensate();
  }
  static constexpr auto get_parameter (highpass_tag)
  {
    return float_param ("", 0.f, 1.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct dry_wet_tag {};
  void                  set (dry_wet_tag, int v) { D = v / 100.; }
  static constexpr auto get_parameter (dry_wet_tag)
  {
    return float_param ("%", 0.f, 100.f, 100.f, 0.01f, 0.5f);
  }
  using parameters = mp_list<density_tag, highpass_tag, dry_wet_tag>;
  //----------------------------------------------------------------------------
private:
#if AIRWINDOWS_FP_DITHER_ENABLE
  uint32_t fpdL;
  uint32_t fpdR;
#endif
  // default stuff
  long double last3SampleL;
  long double last2SampleL;
  long double last1SampleL;
  long double ataAL;
  long double ataBL;
  long double ataCL;
  long double lastDiffSampleL;
  long double iirSampleAL;
  long double iirSampleBL;

  long double last3SampleR;
  long double last2SampleR;
  long double last1SampleR;
  long double ataAR;
  long double ataBR;
  long double ataCR;
  long double lastDiffSampleR;
  long double iirSampleAR;
  long double iirSampleBR;

  float A;
  float B;
  float C;
  float D;
  // Addon
  float sample_rate;
};
//------------------------------------------------------------------------------
}} // namespace artv::airwindows
