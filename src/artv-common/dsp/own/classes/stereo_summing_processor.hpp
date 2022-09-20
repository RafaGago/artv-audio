#pragma once

#include "artv-common/misc/xspan.hpp"
#include <array>

namespace artv {

struct stereo_summing_processor {
  virtual ~stereo_summing_processor() {}

  virtual void process (
    xspan<float*>       outs,
    xspan<float const*> ins,
    unsigned            samples,
    bool                sum)
    = 0;

  virtual void process (
    xspan<double*>       outs,
    xspan<double const*> ins,
    unsigned             samples,
    bool                 sum)
    = 0;
};

} // namespace artv
