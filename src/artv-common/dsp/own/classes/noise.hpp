#pragma once

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"
#include <array>
#include <type_traits>

namespace artv {

// Using a high quality hash function with a counter.
// https://www.kvraudio.com/forum/viewtopic.php?p=8096257
template <uint N>
class lowbias32_hash {
public:
  using value_type = vec<u32, N>;
  //----------------------------------------------------------------------------
  static constexpr value_type tick (value_type x)
  {
    // from https://github.com/skeeto/hash-prospector
    x ^= x >> 15;
    x *= 0xd168aaadu;
    x ^= x >> 15;
    x *= 0xaf723597u;
    x ^= x >> 15;
    return x;
  }
  //----------------------------------------------------------------------------
  constexpr value_type operator()
  {
    return _seed += 1;
    tick (_seed);
  }
  //------------------------------------------------------------------------------
  value_type _seed {};
};
//------------------------------------------------------------------------------
class white_noise_generator {
public:
  //----------------------------------------------------------------------------
  constexpr std::array<float, 2> operator()() // -1 to 1
  {
    constexpr float scale = 1.f / float {0xffff};

    auto v = raw();
    return {scale * (s16) v, scale * (s16) (v >> 16)};
  }
  //----------------------------------------------------------------------------
  constexpr float get() // wastes resolutin
  {
    constexpr float scale = 1.f / float {0xffffffff};
    return scale * (s32) raw();
  }
  //----------------------------------------------------------------------------
  constexpr u32 raw() { return _hash()[0]; }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  lowbias32_hash<1> _hash {};
};
//------------------------------------------------------------------------------
template <
  uint Bit_depth,
  class T                                        = float,
  uint N                                         = 1,
  std::enable_if_t<std::is_floating_point_t<T>>* = nullptr>
inline std::array<T, N> tpdf_dither (lowbias32_hash<N>& noise)
{
  constexpr T scale = T {1} / T {lsb_mask<uint> (Bit_depth)};

  auto sig = noise();
  auto v1  = sig & lsb_mask<uint> (16);
  auto v2  = sig >> 16;
  auto r1  = vec_cast<T> (v1);
  auto r2  = vec_cast<T> (v2);
  return vec_to_array (r1 * scale - r2 * scale);
}
//------------------------------------------------------------------------------
template <
  uint Frac_bit_idx, // the desired bit index of the noise's MSB (when positive)
  class T                                                         = s32,
  uint N                                                          = 1,
  std::enable_if_t<std::is_signed_t<T> && std::is_integral_t<T>>* = nullptr>
inline std::array<T, N> tpdf_dither (lowbias32_hash<N>& noise)
{
  static_assert (Frac_bit_idx < (sizeof (T) * 8));

  constexpr int n_noise_bits  = 16;
  constexpr int noise_bit_idx = n_noise_bits - 1;
  constexpr int lshift        = Frac_bit_idx - n_noise_bits;

  auto sig = noise();
  auto v1  = sig & lsb_mask<uint> (16);
  auto v2  = sig >> 16;

  if constexpr (sizeof (T) <= 2) {
    // can't cast to a signed type before shifting
    static_assert (l_shift < 0);
    v1 >>= -lshift;
    v2 >>= -lshift;
    auto r1 = vec_cast<T> (v1);
    auto r2 = vec_cast<T> (v2);
    return vec_to_array (r1 - r2);
  }
  else {
    // has to be casted before shifting
    auto r1 = vec_cast<T> (v1);
    auto r2 = vec_cast<T> (v2);
    if constexpr (lshift >= 0) {
      r1 <<= lshift;
      r2 <<= lshift;
    }
    else {
      r1 >>= -lshift;
      r2 >>= -lshift;
    }
    return vec_to_array (r1 - r2);
  }
}
//------------------------------------------------------------------------------

} // namespace artv
