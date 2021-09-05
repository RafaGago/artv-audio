// This file has to be compiled as C++11/14 to
#include <cstring>

#include "artv-common/dsp/dragonfly/compilation_firewalls.hpp"

using uint = unsigned int;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"

// clang-format off
#include "dragonfly-hall-reverb/DistrhoPluginInfo.h"
#include "dragonfly-hall-reverb/DSP.hpp"
// clang-format on
#pragma GCC diagnostic pop

namespace artv { namespace dragonfly {

static ::dragonfly::hall::DragonflyReverbDSP* cast (void* p)
{
  return reinterpret_cast<::dragonfly::hall::DragonflyReverbDSP*> (p);
}

hall_compile_firewall::hall_compile_firewall()
{
  _dsp = new ::dragonfly::hall::DragonflyReverbDSP {44100};
}

hall_compile_firewall::~hall_compile_firewall()
{
  if (_dsp) {
    delete cast (_dsp);
    _dsp = nullptr;
  }
}

void hall_compile_firewall::set_parameter (unsigned id, float v)
{
  cast (_dsp)->setParameterValue (id, v);
};

void hall_compile_firewall::reset (unsigned samplerate)
{
  cast (_dsp)->sampleRateChanged (samplerate);
};

void hall_compile_firewall::process (
  float**       out,
  const float** in,
  unsigned      samples)
{
  cast (_dsp)->run (in, out, samples);
};
}} // namespace artv::dragonfly
