#pragma once

#include "artv-common/dsp/own/mufft.hpp"

namespace artv { namespace jsfx {

using fft = mufft::initialized_ffts<
  float,
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
