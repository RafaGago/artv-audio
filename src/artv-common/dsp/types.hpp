#pragma once
namespace artv {

enum class dsp_types {
  console,
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
  other
};

enum class bus_types {
  dummy,
  mono,
  stereo,
};

} // namespace artv
