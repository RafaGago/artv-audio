#pragma once

#include <array>

namespace artv {

struct stereo_summing_processor {
  virtual ~stereo_summing_processor() {}

  virtual void process (
    std::array<float*, 2>       dst,
    std::array<float const*, 2> src,
    unsigned                    samples,
    bool                        sum)
    = 0;

  virtual void process (
    std::array<double*, 2>       dst,
    std::array<double const*, 2> src,
    unsigned                     samples,
    bool                         sum)
    = 0;
};

} // namespace artv
