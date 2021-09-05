#pragma once

#include <array>

namespace artv {

struct stereo_summing_processor {
  virtual ~stereo_summing_processor() {}

  virtual void process (
    crange<float*>       outs,
    crange<float const*> ins,
    unsigned             samples,
    bool                 sum)
    = 0;

  virtual void process (
    crange<double*>       outs,
    crange<double const*> ins,
    unsigned              samples,
    bool                  sum)
    = 0;
};

} // namespace artv
