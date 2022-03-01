// This file has to be compiled as C++11/14 to
#include <cstring>

#include "artv-common/dsp/third_party/dragonfly/compilation_firewalls.hpp"

using uint = unsigned int;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"

// clang-format off
#include "dragonfly-room-reverb/DistrhoPluginInfo.h"
#include "dragonfly-room-reverb/DSP.hpp"
// clang-format on
#pragma GCC diagnostic pop

namespace artv { namespace dragonfly {

static ::dragonfly::room::DragonflyReverbDSP* cast (void* p)
{
  return reinterpret_cast<::dragonfly::room::DragonflyReverbDSP*> (p);
}

room_compile_firewall::room_compile_firewall()
{
  _dsp = new ::dragonfly::room::DragonflyReverbDSP {44100};
}

room_compile_firewall::~room_compile_firewall()
{
  if (_dsp) {
    delete cast (_dsp);
    _dsp = nullptr;
  }
}

void room_compile_firewall::set_parameter (unsigned id, float v)
{
  cast (_dsp)->setParameterValue (id, v);
};

void room_compile_firewall::reset (unsigned samplerate)
{
  cast (_dsp)->sampleRateChanged (samplerate);
};

void room_compile_firewall::process (
  float**       out,
  const float** in,
  unsigned      samples)
{
  cast (_dsp)->run (in, out, samples);
};
}} // namespace artv::dragonfly
