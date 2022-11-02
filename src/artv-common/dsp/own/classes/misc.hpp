#pragma once

#include "artv-common/misc/short_ints.hpp"
#include <array>

namespace artv {

//------------------------------------------------------------------------------
// Using a high quality hash function with a counter.
// https://www.kvraudio.com/forum/viewtopic.php?p=8096257
class white_noise_generator {
public:
  //----------------------------------------------------------------------------
  std::array<float, 2> operator()() // -1 to 1
  {
    constexpr float scale = 1.0f / 0xffff;

    auto v = raw();
    return {scale * (s16) v, scale * (s16) (v >> 16)};
  }
  //----------------------------------------------------------------------------
  float get() // wastes resolutin
  {
    constexpr float scale = 1.0f / (float) 0xffffffff;
    return scale * (s32) raw();
  }
  //----------------------------------------------------------------------------
  u32 raw() { return lowbias32 (++_c); }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // from https://github.com/skeeto/hash-prospector
  u32 lowbias32 (u32 x)
  {
    x ^= x >> 15;
    x *= 0xd168aaadu;
    x ^= x >> 15;
    x *= 0xaf723597u;
    x ^= x >> 15;
    return x;
  }
  //----------------------------------------------------------------------------
  u32 _c {};
};
//------------------------------------------------------------------------------

} // namespace artv
