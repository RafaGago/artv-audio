#pragma once

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_util.hpp"
#include <array>
#include <type_traits>

namespace artv {

// Using a high quality hash function with a counter.
// https://www.kvraudio.com/forum/viewtopic.php?p=8096257
template <uint N>
class lowbias32_hash {
public:
  //----------------------------------------------------------------------------
  using value_type = vec<u32, N>;
  lowbias32_hash()
  {
    for (uint i = 0; i < N; ++i) {
      _seed[i] = i;
    }
  }
  //----------------------------------------------------------------------------
  static value_type tick (value_type x)
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
  value_type operator()()
  {
    _seed += N;
    return tick (_seed);
  }
  //------------------------------------------------------------------------------
  value_type _seed;
};
//------------------------------------------------------------------------------
class white_noise_generator {
public:
  //----------------------------------------------------------------------------
  std::array<float, 2> operator()() // -1 to 1
  {
    constexpr float scale = 1.f / float {0xffff};

    auto v = raw();
    return {scale * (s16) v, scale * (s16) (v >> 16)};
  }
  //----------------------------------------------------------------------------
  float get() // wastes resolutin
  {
    constexpr float scale = 1.f / float {0xffffffff};
    return scale * (s32) raw();
  }
  //----------------------------------------------------------------------------
  u32 raw() { return _hash()[0]; }
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
  std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
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
  std::enable_if_t<std::is_signed_v<T> && std::is_integral_v<T>>* = nullptr>
inline std::array<T, N> tpdf_dither (lowbias32_hash<N>& noise)
{
  using T_u = std::make_unsigned_t<T>;
  static_assert (Frac_bit_idx < (sizeof (T) * 8));

  constexpr int n_noise_bits  = 16;
  constexpr int noise_bit_idx = n_noise_bits - 1;
  // the - 1 is to compensate the bit gained on the final subtraction
  constexpr int lshift = Frac_bit_idx - noise_bit_idx - 1;

  static_assert (((Frac_bit_idx + lshift) > 0), "under representable range");
  static_assert (
    ((Frac_bit_idx + lshift) < (sizeof (T) * 8)), "undesired truncation?");

  auto sig = noise();
  auto v1  = vec_cast<s16> (sig & lsb_mask<uint> (16));
  auto v2  = vec_cast<s16> (sig >> 16);

  if constexpr (sizeof (T) <= 2) {
    // don't cast before shifting for simplicity
    static_assert (lshift < 0);
    v1      = ashr<-lshift> (v1);
    v2      = ashr<-lshift> (v2);
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
      r1 = ashr<-lshift> (r1);
      r2 = ashr<-lshift> (r2);
    }
    return vec_to_array (r1 - r2);
  }
}
//------------------------------------------------------------------------------

} // namespace artv
