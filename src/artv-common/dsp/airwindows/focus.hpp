#pragma once

// Airwindows "focus" Port.
//
// Original code left mostly intact, just run through clang-format and adapted
// to the interface of this project.
//
// https://www.patreon.com/airwindows
// https://github.com/airwindows/airwindows
// from commit 7bb7fde13a9c94af242835823b813e6bcd0f20c8

#include "artv-common/dsp/airwindows/common.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace airwindows {

class focus {
public:
  // DSP------------------------------------------------------------------------
  focus() {}
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::exciter;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();

    A = 0.0;
    B = 0.5;
    C = 0.5;
    D = 1.0;
    E = 1.0;
    for (int x = 0; x < 9; x++) {
      figureL[x] = 0.0;
      figureR[x] = 0.0;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, int samples)
  {
    float* in1  = chnls[0];
    float* in2  = chnls[1];
    float* out1 = chnls[0];
    float* out2 = chnls[1];

    //[0] is frequency: 0.000001 to 0.499999 is near-zero to near-Nyquist
    //[1] is resonance, 0.7071 is Butterworth. Also can't be zero
    double boost = A;
    figureL[0]   = figureR[0]
      = 3515.775 / sample_rate; // fixed frequency, 3.515775k
    figureL[1] = figureR[1] = pow (pow (B, 3) * 2, 2) + 0.0001; // resonance
    int    mode             = C;
    double output           = D;
    double wet              = E;

    double K    = tan (M_PI * figureR[0]);
    double norm = 1.0 / (1.0 + K / figureR[1] + K * K);
    figureL[2] = figureR[2] = K / figureR[1] * norm;
    figureL[4] = figureR[4] = -figureR[2];
    figureL[5] = figureR[5] = 2.0 * (K * K - 1.0) * norm;
    figureL[6] = figureR[6] = (1.0 - K / figureR[1] + K * K) * norm;

    while (--samples >= 0) {
      long double inputSampleL = *in1;
      long double inputSampleR = *in2;
#if AIRWINDOWS_FP_DITHER_ENABLE
      if (abs (inputSampleL) < 1.18e-43)
        inputSampleL = fpd * 1.18e-43;
      if (abs (inputSampleR) < 1.18e-43)
        inputSampleR = fpd * 1.18e-43;
#endif
      long double drySampleL = inputSampleL;
      long double drySampleR = inputSampleR;

      inputSampleL = sin (inputSampleL);
      inputSampleR = sin (inputSampleR);
      // encode Console5: good cleanness

      long double tempSample = (inputSampleL * figureL[2]) + figureL[7];
      figureL[7]             = -(tempSample * figureL[5]) + figureL[8];
      figureL[8]   = (inputSampleL * figureL[4]) - (tempSample * figureL[6]);
      inputSampleL = tempSample;

      tempSample   = (inputSampleR * figureR[2]) + figureR[7];
      figureR[7]   = -(tempSample * figureR[5]) + figureR[8];
      figureR[8]   = (inputSampleR * figureR[4]) - (tempSample * figureR[6]);
      inputSampleR = tempSample;

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      // without this, you can get a NaN condition where it spits out DC offset
      // at full blast!
      inputSampleL = asin (inputSampleL);
      inputSampleR = asin (inputSampleR);
      // decode Console5

      long double groundSampleL = drySampleL - inputSampleL; // set up UnBox
      long double groundSampleR = drySampleR - inputSampleR; // set up UnBox
      inputSampleL *= boost; // now, focussed area gets cranked before distort
      inputSampleR *= boost; // now, focussed area gets cranked before distort

      switch (mode) {
      case 0: // Density
        if (inputSampleL > 1.570796326794897)
          inputSampleL = 1.570796326794897;
        if (inputSampleL < -1.570796326794897)
          inputSampleL = -1.570796326794897;
        if (inputSampleR > 1.570796326794897)
          inputSampleR = 1.570796326794897;
        if (inputSampleR < -1.570796326794897)
          inputSampleR = -1.570796326794897;
        // clip to 1.570796326794897 to reach maximum output
        inputSampleL = sin (inputSampleL);
        inputSampleR = sin (inputSampleR);
        break;
      case 1: // Drive
        if (inputSampleL > 1.0)
          inputSampleL = 1.0;
        if (inputSampleL < -1.0)
          inputSampleL = -1.0;
        if (inputSampleR > 1.0)
          inputSampleR = 1.0;
        if (inputSampleR < -1.0)
          inputSampleR = -1.0;
        inputSampleL
          -= (inputSampleL * (abs (inputSampleL) * 0.6) * (abs (inputSampleL) * 0.6));
        inputSampleR
          -= (inputSampleR * (abs (inputSampleR) * 0.6) * (abs (inputSampleR) * 0.6));
        inputSampleL *= 1.6;
        inputSampleR *= 1.6;
        break;
      case 2: // Spiral
        if (inputSampleL > 1.2533141373155)
          inputSampleL = 1.2533141373155;
        if (inputSampleL < -1.2533141373155)
          inputSampleL = -1.2533141373155;
        if (inputSampleR > 1.2533141373155)
          inputSampleR = 1.2533141373155;
        if (inputSampleR < -1.2533141373155)
          inputSampleR = -1.2533141373155;
        // clip to 1.2533141373155 to reach maximum output
        inputSampleL = sin (inputSampleL * abs (inputSampleL))
          / ((abs (inputSampleL) == 0.0) ? 1 : abs (inputSampleL));
        inputSampleR = sin (inputSampleR * abs (inputSampleR))
          / ((abs (inputSampleR) == 0.0) ? 1 : abs (inputSampleR));
        break;
      case 3: // Mojo
        long double mojo;
        mojo = pow (abs (inputSampleL), 0.25);
        if (mojo > 0.0)
          inputSampleL
            = (sin (inputSampleL * mojo * M_PI * 0.5) / mojo) * 0.987654321;
        mojo = pow (abs (inputSampleR), 0.25);
        if (mojo > 0.0)
          inputSampleR
            = (sin (inputSampleR * mojo * M_PI * 0.5) / mojo) * 0.987654321;
        // mojo is the one that flattens WAAAAY out very softly before
        // wavefolding
        break;
      case 4: // Dyno
        long double dyno;
        dyno = pow (abs (inputSampleL), 4);
        if (dyno > 0.0)
          inputSampleL = (sin (inputSampleL * dyno) / dyno) * 1.1654321;
        dyno = pow (abs (inputSampleR), 4);
        if (dyno > 0.0)
          inputSampleR = (sin (inputSampleR * dyno) / dyno) * 1.1654321;
        // dyno is the one that tries to raise peak energy
        break;
      }

      inputSampleL *= output;
      inputSampleR *= output;

      inputSampleL += groundSampleL; // effectively UnBox
      inputSampleR += groundSampleR; // effectively UnBox

      if (wet != 1.0) {
        inputSampleL = (inputSampleL * wet) + (drySampleL * (1.0 - wet));
        inputSampleR = (inputSampleR * wet) + (drySampleR * (1.0 - wet));
      }
#if AIRWINDOWS_FP_DITHER_ENABLE
      // begin 64 bit stereo floating point dither
      int expon;
      frexp ((double) inputSampleL, &expon);
      fpd ^= fpd << 13;
      fpd ^= fpd >> 17;
      fpd ^= fpd << 5;
      inputSampleL
        += ((double (fpd) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      frexp ((double) inputSampleR, &expon);
      fpd ^= fpd << 13;
      fpd ^= fpd >> 17;
      fpd ^= fpd << 5;
      inputSampleR
        += ((double (fpd) - uint32_t (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
      // end 64 bit stereo floating point dither
#endif
      *out1 = inputSampleL;
      *out2 = inputSampleR;

      *in1++;
      *in2++;
      *out1++;
      *out2++;
    }
  }
  // Parameters (call once per block) ------------------------------------------
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    // no need for separate input and output gain at this level
    A = db_to_gain (v);
    D = 1. / A;
  }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct focus_tag {};
  void                  set (focus_tag, float v) { B = v; }
  static constexpr auto get_parameter (focus_tag)
  {
    return float_param ("", 0.f, 1.f, 0.5f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void                  set (mode_tag, int v) { C = v; }
  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (
      0, make_cstr_array ("Density", "Drive", "Spiral", "Mojo", "Dyno"));
  }
  //----------------------------------------------------------------------------
  struct dry_wet_tag {};
  void                  set (dry_wet_tag, int v) { E = v / 100.; }
  static constexpr auto get_parameter (dry_wet_tag)
  {
    return float_param ("%", 0.f, 100.f, 100.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag, mode_tag, focus_tag, dry_wet_tag>;
  //----------------------------------------------------------------------------
private:
  long double figureL[9];
  long double figureR[9];
#if AIRWINDOWS_FP_DITHER_ENABLE
  uint32_t fpd;
#endif
  // default stuff

  float A; // Drive
  float B; // Focus
  int   C; // Mode
  float D; // inv drive
  float E; // dry-wet

  // Addon
  float sample_rate;
};

}} // namespace artv::airwindows
