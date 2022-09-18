#pragma once

// Airwindows console5 Port.
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
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace airwindows {

class console5bus {
public:
  // DSP------------------------------------------------------------------------
  console5bus() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate     = pc.get_sample_rate();
    A               = 1.0;
    lastSampleBussL = 0.0;
    lastSampleBussR = 0.0;
    lastFXBussL     = 0.0;
    lastFXBussR     = 0.0;
    iirCorrectL     = 0.0;
    iirCorrectR     = 0.0;
    gainchase       = -90.0;
    settingchase    = -90.0;
    chasespeed      = 350.0;
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpNShapeL = 0.0;
    fpNShapeR = 0.0;
#endif
  }
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    T const* in1  = ins[0];
    T const* in2  = ins[1];
    T*       out1 = outs[0];
    T*       out2 = outs[1];

    double overallscale = 1.0;
    overallscale /= 44100.0;
    overallscale *= sample_rate;

    double      inputgain = A;
    double      differenceL;
    double      differenceR;
    double      nearZeroL;
    double      nearZeroR;
    double      servoTrim = 0.0000001 / overallscale;
    double      bassTrim  = 0.005 / overallscale;
    long double inputSampleL;
    long double inputSampleR;

    if (settingchase != inputgain) {
      chasespeed *= 2.0;
      settingchase = inputgain;
    }
    if (chasespeed > 2500.0)
      chasespeed = 2500.0;
    if (gainchase < 0.0)
      gainchase = inputgain;

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
      chasespeed *= 0.9999;
      chasespeed -= 0.01;
      if (chasespeed < 350.0)
        chasespeed = 350.0;
      // we have our chase speed compensated for recent fader activity

      gainchase = (((gainchase * chasespeed) + inputgain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor

      inputSampleL *= gainchase;
      inputSampleR *= gainchase;
      // done with trim control

      if (inputSampleL > 1.0)
        inputSampleL = 1.0;
      if (inputSampleL < -1.0)
        inputSampleL = -1.0;
      inputSampleL = asin (inputSampleL);
      // amplitude aspect

      if (inputSampleR > 1.0)
        inputSampleR = 1.0;
      if (inputSampleR < -1.0)
        inputSampleR = -1.0;
      inputSampleR = asin (inputSampleR);
      // amplitude aspect

      differenceL     = lastSampleBussL - inputSampleL;
      differenceR     = lastSampleBussR - inputSampleR;
      lastSampleBussL = inputSampleL;
      lastSampleBussR = inputSampleR;
      // derive slew part off direct sample measurement + from last time

      if (differenceL > 1.57079633)
        differenceL = 1.57079633;
      if (differenceL < -1.57079633)
        differenceL = -1.57079633;
      if (differenceR > 1.57079633)
        differenceR = 1.57079633;
      if (differenceR < -1.57079633)
        differenceR = -1.57079633;

      differenceL = lastFXBussL + sin (differenceL);
      differenceR = lastFXBussR + sin (differenceR);
      // we're about to use this twice and then not use difference again, so
      // we'll reuse it enhance slew is arcsin(): cutting it back is sin()

      iirCorrectL += inputSampleL - differenceL;
      iirCorrectR += inputSampleR - differenceR;
      inputSampleL = differenceL;
      inputSampleR = differenceR;
      // apply the slew to stored value: can develop DC offsets.
      // store the change we made so we can dial it back

      lastFXBussL = inputSampleL;
      lastFXBussR = inputSampleR;
      if (lastFXBussL > 1.0)
        lastFXBussL = 1.0;
      if (lastFXBussL < -1.0)
        lastFXBussL = -1.0;
      if (lastFXBussR > 1.0)
        lastFXBussR = 1.0;
      if (lastFXBussR < -1.0)
        lastFXBussR = -1.0;
      // build new signal off what was present in output last time

      nearZeroL = pow (abs (abs (lastFXBussL) - 1.0), 2);
      nearZeroR = pow (abs (abs (lastFXBussR) - 1.0), 2);
      // if the sample is very near zero this number is higher.
      if (iirCorrectL > 0)
        iirCorrectL -= servoTrim;
      if (iirCorrectL < 0)
        iirCorrectL += servoTrim;
      if (iirCorrectR > 0)
        iirCorrectR -= servoTrim;
      if (iirCorrectR < 0)
        iirCorrectR += servoTrim;
      // cut back the servo by which we're pulling back the DC
      lastFXBussL += (iirCorrectL * 0.0000005);
      lastFXBussR += (iirCorrectR * 0.0000005);
      // apply the servo to the stored value, pulling back the DC
      lastFXBussL *= (1.0 - (nearZeroL * bassTrim));
      lastFXBussR *= (1.0 - (nearZeroR * bassTrim));
      // this cuts back the DC offset directly, relative to how near zero we are

#if AIRWINDOWS_FP_DITHER_ENABLE
      // stereo 64 bit dither, made small and tidy.
      int expon;
      frexp ((double) inputSampleL, &expon);
      long double dither
        = (rand() / (RAND_MAX * 7.737125245533627e+25)) * pow (2, expon + 62);
      dither /= 536870912.0; // needs this to scale to 64 bit zone
      inputSampleL += (dither - fpNShapeL);
      fpNShapeL = dither;
      frexp ((double) inputSampleR, &expon);
      dither
        = (rand() / (RAND_MAX * 7.737125245533627e+25)) * pow (2, expon + 62);
      dither /= 536870912.0; // needs this to scale to 64 bit zone
      inputSampleR += (dither - fpNShapeR);
      fpNShapeR = dither;
      // end 64 bit dither
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
  double lastSampleBussL;
  double lastFXBussL;
  double lastSampleBussR;
  double lastFXBussR;
  double iirCorrectL;
  double iirCorrectR;
  double gainchase;
  double settingchase;
  double chasespeed;
#if AIRWINDOWS_FP_DITHER_ENABLE
  long double fpNShapeL;
  long double fpNShapeR;
#endif

  float A; // gain
  float sample_rate;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
class console5channel : public stereo_summing_processor {
public:
  // DSP------------------------------------------------------------------------
  console5channel() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate        = pc.get_sample_rate();
    A                  = 1.0;
    lastSampleChannelL = 0.0;
    lastSampleChannelR = 0.0;
    lastFXChannelL     = 0.0;
    lastFXChannelR     = 0.0;
    iirCorrectL        = 0.0;
    iirCorrectR        = 0.0;
    gainchase          = -90.0;
    settingchase       = -90.0;
    chasespeed         = 350.0;
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpNShapeL = 0.0;
    fpNShapeR = 0.0;
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

    double overallscale = 1.0;
    overallscale /= 44100.0;
    overallscale *= sample_rate;

    double      inputgain = A;
    double      differenceL;
    double      differenceR;
    double      nearZeroL;
    double      nearZeroR;
    double      servoTrim = 0.0000001 / overallscale;
    double      bassTrim  = 0.005 / overallscale;
    long double inputSampleL;
    long double inputSampleR;

    if (settingchase != inputgain) {
      chasespeed *= 2.0;
      settingchase = inputgain;
    }
    if (chasespeed > 2500.0)
      chasespeed = 2500.0;
    if (gainchase < 0.0)
      gainchase = inputgain;

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
      chasespeed *= 0.9999;
      chasespeed -= 0.01;
      if (chasespeed < 350.0)
        chasespeed = 350.0;
      // we have our chase speed compensated for recent fader activity

      gainchase = (((gainchase * chasespeed) + inputgain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor

      if (1.0 != gainchase) {
        inputSampleL *= gainchase;
        inputSampleR *= gainchase;
      }
      // done with trim control

      differenceL        = lastSampleChannelL - inputSampleL;
      lastSampleChannelL = inputSampleL;
      differenceR        = lastSampleChannelR - inputSampleR;
      lastSampleChannelR = inputSampleR;
      // derive slew part off direct sample measurement + from last time

      if (differenceL > 1.0)
        differenceL = 1.0;
      if (differenceL < -1.0)
        differenceL = -1.0;
      if (differenceR > 1.0)
        differenceR = 1.0;
      if (differenceR < -1.0)
        differenceR = -1.0;
      // clamp the slew correction to prevent invalid math results

      differenceL = lastFXChannelL + asin (differenceL);
      differenceR = lastFXChannelR + asin (differenceR);
      // we're about to use this twice and then not use difference again, so
      // we'll reuse it enhance slew is arcsin(): cutting it back is sin()

      iirCorrectL += inputSampleL - differenceL;
      inputSampleL = differenceL;
      iirCorrectR += inputSampleR - differenceR;
      inputSampleR = differenceR;
      // apply the slew to stored value: can develop DC offsets.
      // store the change we made so we can dial it back

      lastFXChannelL = inputSampleL;
      lastFXChannelR = inputSampleR;
      if (lastFXChannelL > 1.0)
        lastFXChannelL = 1.0;
      if (lastFXChannelL < -1.0)
        lastFXChannelL = -1.0;
      if (lastFXChannelR > 1.0)
        lastFXChannelR = 1.0;
      if (lastFXChannelR < -1.0)
        lastFXChannelR = -1.0;
      // store current sample as new base for next offset

      nearZeroL = pow (abs (abs (lastFXChannelL) - 1.0), 2);
      nearZeroR = pow (abs (abs (lastFXChannelR) - 1.0), 2);
      // if the sample is very near zero this number is higher.
      if (iirCorrectL > 0)
        iirCorrectL -= servoTrim;
      if (iirCorrectL < 0)
        iirCorrectL += servoTrim;
      if (iirCorrectR > 0)
        iirCorrectR -= servoTrim;
      if (iirCorrectR < 0)
        iirCorrectR += servoTrim;
      // cut back the servo by which we're pulling back the DC
      lastFXChannelL += (iirCorrectL * 0.0000005);
      lastFXChannelR += (iirCorrectR * 0.0000005);
      // apply the servo to the stored value, pulling back the DC
      lastFXChannelL *= (1.0 - (nearZeroL * bassTrim));
      lastFXChannelR *= (1.0 - (nearZeroR * bassTrim));
      // this cuts back the DC offset directly, relative to how near zero we are

      if (inputSampleL > 1.57079633)
        inputSampleL = 1.57079633;
      if (inputSampleL < -1.57079633)
        inputSampleL = -1.57079633;
      inputSampleL = sin (inputSampleL);
      // amplitude aspect

      if (inputSampleR > 1.57079633)
        inputSampleR = 1.57079633;
      if (inputSampleR < -1.57079633)
        inputSampleR = -1.57079633;
      inputSampleR = sin (inputSampleR);
      // amplitude aspect

#if AIRWINDOWS_FP_DITHER_ENABLE
      // stereo 64 bit dither, made small and tidy.
      int expon;
      frexp ((double) inputSampleL, &expon);
      long double dither
        = (rand() / (RAND_MAX * 7.737125245533627e+25)) * pow (2, expon + 62);
      dither /= 536870912.0; // needs this to scale to 64 bit zone
      inputSampleL += (dither - fpNShapeL);
      fpNShapeL = dither;
      frexp ((double) inputSampleR, &expon);
      dither
        = (rand() / (RAND_MAX * 7.737125245533627e+25)) * pow (2, expon + 62);
      dither /= 536870912.0; // needs this to scale to 64 bit zone
      inputSampleR += (dither - fpNShapeR);
      fpNShapeR = dither;
// end 64 bit dither
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
  void set (console5bus::drive_tag, float value)
  {
    A = db_to_gain (value);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<>;
  //----------------------------------------------------------------------------
private:
  double lastSampleChannelL;
  double lastSampleChannelR;
  double lastFXChannelL;
  double lastFXChannelR;
  double iirCorrectL;
  double iirCorrectR;
  double gainchase;
  double settingchase;
  double chasespeed;
#if AIRWINDOWS_FP_DITHER_ENABLE
  long double fpNShapeL;
  long double fpNShapeR;
#endif
  float A; // gain
  float sample_rate;
};
}} // namespace artv::airwindows
