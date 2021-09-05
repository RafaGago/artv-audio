#pragma once

// Airwindows Channel4 Port.
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

class channel4 {
public:
  // DSP------------------------------------------------------------------------
  channel4() {}
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::exciter;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    drive_p     = 0.0;
    gain        = 1.f;
    fpNShapeLA  = 0.0;
    fpNShapeLB  = 0.0;
    fpNShapeRA  = 0.0;
    fpNShapeRB  = 0.0;
    fpFlip      = true;
    iirSampleLA = 0.0;
    iirSampleRA = 0.0;
    iirSampleLB = 0.0;
    iirSampleRB = 0.0;
    lastSampleL = 0.0;
    lastSampleR = 0.0;
    iirAmount   = 0.005832;
    threshold   = 0.33362176; // instantiating with Neve values
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
    double fpTemp; // this is different from singlereplacing
    double fpOld = 0.618033988749894848204586; // golden ratio!
    double fpNew = 1.0 - fpOld;

    const double localiirAmount = iirAmount / overallscale;
    const double localthreshold = threshold / overallscale;
    const double density        = pow (
      drive_p, 2); // this doesn't relate to the plugins Density and Drive much
    double      clamp;
    long double bridgerectifier;

    long double inputSampleL;
    long double inputSampleR;

    while (--samples < ((uint) -1ull)) {
      inputSampleL = *in1;
      inputSampleR = *in2;
#if AIRWINDOWS_DENORMALIZATION_ENABLE
      if (inputSampleL < 1.2e-38 && -inputSampleL < 1.2e-38) {
        static int noisesource = 0;
        // this declares a variable before anything else is compiled. It won't
        // keep assigning it to 0 for every sample, it's as if the declaration
        // doesn't exist in this context, but it lets me add this
        // denormalization fix in a single place rather than updating it in
        // three different locations. The variable isn't thread-safe but this is
        // only a random seed and we can share it with whatever.
        noisesource = noisesource % 1700021;
        noisesource++;
        int residue = noisesource * noisesource;
        residue     = residue % 170003;
        residue *= residue;
        residue = residue % 17011;
        residue *= residue;
        residue = residue % 1709;
        residue *= residue;
        residue = residue % 173;
        residue *= residue;
        residue             = residue % 17;
        double applyresidue = residue;
        applyresidue *= 0.00000001;
        applyresidue *= 0.00000001;
        inputSampleL = applyresidue;
      }
      if (inputSampleR < 1.2e-38 && -inputSampleR < 1.2e-38) {
        static int noisesource = 0;
        noisesource            = noisesource % 1700021;
        noisesource++;
        int residue = noisesource * noisesource;
        residue     = residue % 170003;
        residue *= residue;
        residue = residue % 17011;
        residue *= residue;
        residue = residue % 1709;
        residue *= residue;
        residue = residue % 173;
        residue *= residue;
        residue             = residue % 17;
        double applyresidue = residue;
        applyresidue *= 0.00000001;
        applyresidue *= 0.00000001;
        inputSampleR = applyresidue;
        // this denormalization routine produces a white noise at -300 dB which
        // the noise shaping will interact with to produce a bipolar output, but
        // the noise is actually all positive. That should stop any variables
        // from going denormal, and the routine only kicks in if digital black
        // is input. As a final touch, if you save to 24-bit the silence will
        // return to being digital black again.
      }
#endif
      if (fpFlip) {
        iirSampleLA = (iirSampleLA * (1 - localiirAmount))
          + (inputSampleL * localiirAmount);
        inputSampleL = inputSampleL - iirSampleLA;
        iirSampleRA  = (iirSampleRA * (1 - localiirAmount))
          + (inputSampleR * localiirAmount);
        inputSampleR = inputSampleR - iirSampleRA;
      }
      else {
        iirSampleLB = (iirSampleLB * (1 - localiirAmount))
          + (inputSampleL * localiirAmount);
        inputSampleL = inputSampleL - iirSampleLB;
        iirSampleRB  = (iirSampleRB * (1 - localiirAmount))
          + (inputSampleR * localiirAmount);
        inputSampleR = inputSampleR - iirSampleRB;
      }
      // highpass section

      bridgerectifier = abs (inputSampleL) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.0;
      else
        bridgerectifier = sin (bridgerectifier);
      if (inputSampleL > 0)
        inputSampleL
          = (inputSampleL * (1 - density)) + (bridgerectifier * density);
      else
        inputSampleL
          = (inputSampleL * (1 - density)) - (bridgerectifier * density);

      bridgerectifier = abs (inputSampleR) * 1.57079633;
      if (bridgerectifier > 1.57079633)
        bridgerectifier = 1.0;
      else
        bridgerectifier = sin (bridgerectifier);
      if (inputSampleR > 0)
        inputSampleR
          = (inputSampleR * (1 - density)) + (bridgerectifier * density);
      else
        inputSampleR
          = (inputSampleR * (1 - density)) - (bridgerectifier * density);
      // drive section

      clamp = inputSampleL - lastSampleL;
      if (clamp > localthreshold)
        inputSampleL = lastSampleL + localthreshold;
      if (-clamp > localthreshold)
        inputSampleL = lastSampleL - localthreshold;
      lastSampleL = inputSampleL;

      clamp = inputSampleR - lastSampleR;
      if (clamp > localthreshold)
        inputSampleR = lastSampleR + localthreshold;
      if (-clamp > localthreshold)
        inputSampleR = lastSampleR - localthreshold;
      lastSampleR = inputSampleR;
      // slew section
      // noise shaping to 64-bit floating point
      if (fpFlip) {
        fpTemp     = inputSampleL;
        fpNShapeLA = (fpNShapeLA * fpOld) + ((inputSampleL - fpTemp) * fpNew);
        inputSampleL += fpNShapeLA;
        fpTemp     = inputSampleR;
        fpNShapeRA = (fpNShapeRA * fpOld) + ((inputSampleR - fpTemp) * fpNew);
        inputSampleR += fpNShapeRA;
      }
      else {
        fpTemp     = inputSampleL;
        fpNShapeLB = (fpNShapeLB * fpOld) + ((inputSampleL - fpTemp) * fpNew);
        inputSampleL += fpNShapeLB;
        fpTemp     = inputSampleR;
        fpNShapeRB = (fpNShapeRB * fpOld) + ((inputSampleR - fpTemp) * fpNew);
        inputSampleR += fpNShapeRB;
      }
      fpFlip = !fpFlip;
      // end noise shaping on 64 bit output
      *out1 = inputSampleL * gain;
      *out2 = inputSampleR * gain;

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
      break; // Neve
    case 1:
      iirAmount = 0.004096;
      threshold = 0.59969536;
      break; // API
    case 2:
      iirAmount = 0.004913;
      threshold = 0.84934656;
      break; // SSL
    default:
      break; // should not happen
    }
  }
  static constexpr auto get_parameter (type_tag)
  {
    return choice_param (0, make_cstr_array ("Neve", "API ", "SSL"));
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    drive_p = v;
    gain    = 1. - (0.33 * v); // trying to compensate for loudness increases
  }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("", 0.f, 1.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = typelist<type_tag, drive_tag>;
  //----------------------------------------------------------------------------
private:
  float sample_rate = 0.;

  long double fpNShapeLA;
  long double fpNShapeLB;
  long double fpNShapeRA;
  long double fpNShapeRB;
  bool        fpFlip;
  // default stuff
  double iirSampleLA;
  double iirSampleRA;
  double iirSampleLB;
  double iirSampleRB;
  double lastSampleL;
  double lastSampleR;
  double iirAmount;
  double threshold;

  float drive_p; // parameters. Always 0-1, and we scale/alter them elsewhere.
  float gain;
};

}} // namespace artv::airwindows
