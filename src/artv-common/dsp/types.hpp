#pragma once
namespace artv {

enum class dsp_types {
  console,
  crossover,
  mixer,
  delay,
  distortion,
  saturation,
  dynamics,
  filter,
  eq,
  exciter,
  modulation,
  reverb,
  stereo,
  waveshaper,
  upsampler,
  downsampler,
  pitch,
  other
};

enum class bus_types {
  dummy,
  mono,
  stereo,
};

} // namespace artv
