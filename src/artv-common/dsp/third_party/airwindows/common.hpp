#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/juce/math.hpp"
#include "artv-common/misc/short_ints.hpp"

#define AIRWINDOWS_FP_DITHER_ENABLE 0
#define AIRWINDOWS_DENORMALIZATION_ENABLE 0
#define AIRWINDOWS_SUBSONIC_FILTER_ENABLE 0

namespace artv { namespace airwindows {

#if AIRWINDOWS_FP_DITHER_ENABLE

// incomplete!!! If this is enabled it needs work

// Fx declarations
struct channel9;
struct console7bus;
struct console7channel;
//------------------------------------------------------------------------------
struct dither1 {
  dither1()
  {
    fpd = 1.0;
    while (fpd < 16386) {
      fpd = rand() * UINT32_MAX;
    }
  }
  //----------------------------------------------------------------------------
  long double clamp (long double v, float)
  {
    if (abs (v) < 1.18e-37) {
      v = fpd * 1.18e-37;
    }
    return v;
  }
  //----------------------------------------------------------------------------
  long double operator() (long double v, float)
  {
    int expon;
    frexpf ((float) v, &expon);
    fpd ^= fpd << 13;
    fpd ^= fpd >> 17;
    fpd ^= fpd << 5;
    v += ((double (fpd) - u32 {0x7fffffff}) * 5.5e-36l * pow (2, expon + 62));
    return v;
  }
  //----------------------------------------------------------------------------
  long double clamp (long double v, double)
  {
    if (abs (v) < 1.18e-43) {
      v = fpd * 1.18e-43;
    }
    return v;
  }
  //----------------------------------------------------------------------------
  long double operator() (long double v, double)
  {
    int expon;
    frexp ((double) v, &expon);
    fpd ^= fpd << 13;
    fpd ^= fpd >> 17;
    fpd ^= fpd << 5;
    v += ((double (fpd) - u32 (0x7fffffff)) * 1.1e-44l * pow (2, expon + 62));
    return v;
  }
  //----------------------------------------------------------------------------
private:
  u32 fpd;
};
//------------------------------------------------------------------------------
template <class T>
struct get_dither_for;

template <>
struct get_dither_for<channel9> {
  using type = dither1;
};

template <>
struct get_dither_for<console7bus> {
  using type = dither1;
};

template <>
struct get_dither_for<console7channel> {
  using type = dither1;
};
#endif // #AIRWINDOWS_FP_DITHER_ENABLE

}} // namespace artv::airwindows
