// This file has to be compiled as C++11/14 to

#include "artv-common/dsp/dragonfly/compilation_firewalls.hpp"

using uint = unsigned int;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"

// clang-format off
#include "dragonfly-plate-reverb/DistrhoPluginInfo.h"
#include "dragonfly-plate-reverb/DSP.hpp"
// clang-format on
#pragma GCC diagnostic pop

namespace artv { namespace dragonfly {

static ::dragonfly::plate::DragonflyReverbDSP* cast (void* p)
{
  return reinterpret_cast<::dragonfly::plate::DragonflyReverbDSP*> (p);
}

plate_compile_firewall::plate_compile_firewall()
{
  _dsp = new ::dragonfly::plate::DragonflyReverbDSP {44100};
}

plate_compile_firewall::~plate_compile_firewall()
{
  if (_dsp) {
    auto p = cast (_dsp);
    delete p;
    _dsp = nullptr;
  }
}

void plate_compile_firewall::set_parameter (unsigned id, float v)
{
  cast (_dsp)->setParameterValue (id, v);
};

void plate_compile_firewall::reset (unsigned samplerate)
{
  cast (_dsp)->sampleRateChanged (samplerate);
};

void plate_compile_firewall::process (float** io, unsigned samples)
{
  cast (_dsp)->run ((const float**) io, io, samples);
};
}} // namespace artv::dragonfly
