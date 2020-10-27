#pragma once

// Airwindows console4 Port.
//
// Original code left mostly intact, just run through clang-format and adapted
// to the interface of this project.
//
// https://www.patreon.com/airwindows
// https://github.com/airwindows/airwindows
// from commit 7bb7fde13a9c94af242835823b813e6bcd0f20c8

#include "artv-common/dsp/airwindows/common.hpp"
#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/own/stereo_summing_processor.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace airwindows {

class console4bus {
public:
  // DSP------------------------------------------------------------------------
  console4bus() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate  = pc.get_sample_rate();
    gain         = 1.0;
    lastSampleL  = 0.0;
    lastSampleR  = 0.0;
    gainchase    = -90.0;
    settingchase = -90.0;
    chasespeed   = 350.0;
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpNShapeL = 0.0;
    fpNShapeR = 0.0;
#endif
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, int samples)
  {
    T const* in1  = chnls[0];
    T const* in2  = chnls[1];
    T*       out1 = chnls[0];
    T*       out2 = chnls[1];

    double overallscale = 1.0;
    overallscale /= 44100.0;
    overallscale *= sample_rate;

    long double inputSampleL;
    long double inputSampleR;
    long double half;
    long double falf;
    long double slewcompensation;
    if (settingchase != gain) {
      chasespeed *= 2.0;
      settingchase = gain;
    }
    if (chasespeed > 2500.0)
      chasespeed = 2500.0;
    if (gainchase < 0.0)
      gainchase = gain;

    while (--samples >= 0) {
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

      gainchase = (((gainchase * chasespeed) + gain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor

      if (1.0 != gainchase) {
        inputSampleL *= gainchase;
        inputSampleR *= gainchase;
      }
      // done with trim control

      half = inputSampleL;
      falf = abs (half);
      half *= falf;
      half *= falf;
      slewcompensation = abs (inputSampleL - lastSampleL) * overallscale;
      // magnify effect at high sample rate so it will still register when
      // inter-sample changes are very small at high rates.
      if (slewcompensation > 1.0)
        slewcompensation = 1.0;
      // let's not invert the effect: maximum application is to cancel out half
      // entirely
      half *= (1.0 - slewcompensation);
      // apply it
      lastSampleL = inputSampleL;
      inputSampleL += half;
      // this is the inverse processing for Console: boosts but not so much if
      // there's slew. is this too subtle an effect?

      half = inputSampleR;
      falf = abs (half);
      half *= falf;
      half *= falf;
      slewcompensation = abs (inputSampleR - lastSampleR) * overallscale;
      // magnify effect at high sample rate so it will still register when
      // inter-sample changes are very small at high rates.
      if (slewcompensation > 1.0)
        slewcompensation = 1.0;
      // let's not invert the effect: maximum application is to cancel out half
      // entirely
      half *= (1.0 - slewcompensation);
      // apply it
      lastSampleR = inputSampleR;
      inputSampleR += half;
      // this is the inverse processing for Console: boosts but not so much if
      // there's slew. is this too subtle an effect?
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

      *in1++;
      *in2++;
      *out1++;
      *out2++;
    }
  }
  // Parameters (call once per block) ------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float value) { gain = db_to_gain (-value); }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag>;
  //----------------------------------------------------------------------------
private:
#if AIRWINDOWS_FP_DITHER_ENABLE
  double fpNShapeL;
  double fpNShapeR;
#endif
  // default stuff
  double lastSampleL;
  double lastSampleR;
  double gainchase;
  double settingchase;
  double chasespeed;

  float gain;

  float sample_rate;
};
//------------------------------------------------------------------------------
class console4channel : public stereo_summing_processor {
public:
  // DSP------------------------------------------------------------------------
  console4channel() {}
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    sample_rate  = pc.get_sample_rate();
    gain         = 1.0;
    gainchase    = -90.0;
    settingchase = -90.0;
    chasespeed   = 350.0;
#if AIRWINDOWS_FP_DITHER_ENABLE
    fpNShapeL = 0.0;
    fpNShapeR = 0.0;
#endif
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_adding (
    std::array<T*, 2>       dst,
    std::array<T const*, 2> src,
    int                     samples,
    bool                    dst_sum)
  {
    T const* in1  = src[0];
    T const* in2  = src[1];
    T*       out1 = dst[0];
    T*       out2 = dst[1];

    double overallscale = 1.0;
    overallscale /= 44100.0;
    long double inputSampleR;
    long double inputSampleL;
    long double half;
    long double falf;
    // replace inputgain with gain, serves same purpose. Stereo inputsample.
    if (settingchase != gain) {
      chasespeed *= 2.0;
      settingchase = gain;
    }
    if (chasespeed > 2500.0)
      chasespeed = 2500.0;
    if (gainchase < 0.0)
      gainchase = gain;
    // settings section from AU

    while (--samples >= 0) {
      inputSampleL = *in1;
      inputSampleR = *in2;
#if AIRWINDOWS_DENORMALIZATION_PREVENTION
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

      gainchase = (((gainchase * chasespeed) + gain) / (chasespeed + 1.0));
      // gainchase is chasing the target, as a simple multiply gain factor

      if (1.0 != gainchase) {
        inputSampleL *= gainchase;
        inputSampleR *= gainchase;
      }
      // done with trim control

      half = inputSampleL;
      falf = abs (half);
      half *= falf;
      half *= falf;
      inputSampleL -= half;

      half = inputSampleR;
      falf = abs (half);
      half *= falf;
      half *= falf;
      inputSampleR -= half;
// entire audio code. kthxbai!
// this is part of the Purest line: stuff that is on every track
// needs to be DAMN LOW ON MATH srsly guys

// stereo 64 bit dither, made small and tidy.
#if AIRWINDOWS_FP_DITHER_ENABLE
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
#endif
      if (dst_sum) {
        *out1 += inputSampleL;
        *out2 += inputSampleR;
      }
      else {
        *out1 = inputSampleL;
        *out2 = inputSampleR;
      }

      *in1++;
      *in2++;
      *out1++;
      *out2++;
    }
  }
  //----------------------------------------------------------------------------
  void process (
    std::array<float*, 2>       dst,
    std::array<float const*, 2> src,
    uint                        samples,
    bool                        dst_sum) override
  {
    process_block_adding (dst, src, samples, dst_sum);
  }
  //----------------------------------------------------------------------------
  void process (
    std::array<double*, 2>       dst,
    std::array<double const*, 2> src,
    uint                         samples,
    bool                         dst_sum) override
  {
    process_block_adding (dst, src, samples, dst_sum);
  }
  // Parameters (call once per block) ------------------------------------------
  void set (console4bus::drive_tag, float value) { gain = db_to_gain (value); }
  //----------------------------------------------------------------------------
  using parameters = mp_list<>;
  //----------------------------------------------------------------------------
private:
  double gainchase;
  double settingchase;
  double chasespeed;

#if AIRWINDOWS_FP_DITHER_ENABLE
  double fpNShapeL;
  double fpNShapeR;
#endif
  float gain;
  // default stuff
  float sample_rate;
};

}} // namespace artv::airwindows
