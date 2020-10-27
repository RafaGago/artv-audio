// This file has to be compiled as C++11/14 to

#include "artv-common/dsp/dragonfly/compilation_firewalls.hpp"

using uint = unsigned int;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"

// clang-format off
#include "dragonfly-early-reflections/DistrhoPluginInfo.h"
#include "dragonfly-early-reflections/DSP.hpp"
// clang-format on
#pragma GCC diagnostic pop

namespace artv { namespace dragonfly {

static ::dragonfly::er::DragonflyReverbDSP* cast (void* p)
{
  return reinterpret_cast<::dragonfly::er::DragonflyReverbDSP*> (p);
}

er_compile_firewall::er_compile_firewall()
{
  _dsp = new ::dragonfly::er::DragonflyReverbDSP {44100};
}

er_compile_firewall::~er_compile_firewall()
{
  if (_dsp) {
    delete cast (_dsp);
    _dsp = nullptr;
  }
}

void er_compile_firewall::set_parameter (unsigned id, float v)
{
  cast (_dsp)->setParameterValue (id, v);
};

void er_compile_firewall::reset (unsigned samplerate)
{
  cast (_dsp)->sampleRateChanged (samplerate);
};

void er_compile_firewall::process (float** io, unsigned samples)
{
  cast (_dsp)->run ((const float**) io, io, samples);
};
}} // namespace artv::dragonfly
