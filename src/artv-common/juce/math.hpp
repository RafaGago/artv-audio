#pragma once

#include <cmath>
#include <cstdint>
#include <type_traits>

#define JUCE_OPT_MATH 0

#if JUCE_OPT_MATH
#include <juce_audio_processors/juce_audio_processors.h>
using optmath = juce::dsp::FastMathApproximations;
#else
namespace optmath = std;
#endif
