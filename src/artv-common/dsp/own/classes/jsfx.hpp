#pragma once

#include "artv-common/dsp/own/classes/fft.hpp"

namespace artv { namespace jsfx {

using fft = initialized_ffts<
  float,
  true,
  16,
  32,
  64,
  128,
  256,
  512,
  1024,
  2048,
  4096,
  8192,
  16384,
  32768>;

}} // namespace artv::jsfx
