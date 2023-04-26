#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>

#include "artv-common/misc/float.hpp"

namespace artv {
//------------------------------------------------------------------------------
TEST (float_packing, float_dtz_clamp)
{
  using float16
    = f16pack<5, -1, f16pack_dftz | f16pack_clamp | f16pack_nearest>;
  float max = float16::max();
  // correct maximum
  EXPECT_NEAR (0.99999f, max, 0.001f);
  float got = float16::decode (float16::encode (max));
  // regular value
  EXPECT_NEAR (max, got, max * 0.001f);
  got = float16::decode (float16::encode (float16::denormal_min()));
  // denornal
  EXPECT_EQ (0.f, got);
  // clamping
  got = float16::decode (float16::encode (max * 4));
  EXPECT_NEAR (max, got, max * 0.001f);
}
//------------------------------------------------------------------------------
TEST (float_packing, float_clamp)
{
  using float16 = f16pack<5, -1, f16pack_clamp | f16pack_nearest>;
  float max     = float16::max();
  // correct maximum
  EXPECT_NEAR (0.99999f, max, 0.001f);
  float got = float16::decode (float16::encode (max));
  // regular value
  EXPECT_NEAR (max, got, max * 0.001f);
  // denornal
  float f16_denormal = float16::denormal_min();
  got                = float16::decode (float16::encode (f16_denormal));
  EXPECT_EQ (f16_denormal, got);
  got = float16::decode (float16::encode (max * 4));
  EXPECT_NEAR (max, got, max * 0.001f);
}
//------------------------------------------------------------------------------
TEST (float_packing, vec_float_dtfz_clamp)
{
  using float16
    = f16pack<5, -1, f16pack_dftz | f16pack_clamp | f16pack_nearest>;
  f32_x4 expected {float16::max(), -0.51f, -0.251f, -0.1251f};
  f32_x4 got = float16::decode (float16::encode (expected));
  // regular value
  for (uint i = 0; i < 4; ++i) {
    EXPECT_NEAR (expected[i], got[i], abs (expected[i] * 0.001f));
  }
  got = float16::decode (
    float16::encode (vec_set<f32_x4> (float16::denormal_min())));
  for (uint i = 0; i < 4; ++i) {
    EXPECT_EQ (0, got[i]);
  }
  got = float16::decode (
    float16::encode (vec_set<f32_x4> (float16::max() * 4.f)));
  for (uint i = 0; i < 4; ++i) {
    EXPECT_EQ (float16::max(), got[i]);
  }
}
//------------------------------------------------------------------------------
TEST (float_packing, vec_float_clamp)
{
  using float16 = f16pack<5, -1, f16pack_clamp | f16pack_nearest>;
  f32_x4 expected {float16::max(), -0.51f, -0.251f, -0.1251f};
  f32_x4 got = float16::decode (float16::encode (expected));
  // regular value
  for (uint i = 0; i < 4; ++i) {
    EXPECT_NEAR (expected[i], got[i], abs (expected[i] * 0.001f));
  }
  expected = f32_x4 {
    float16::denormal_min(),
    float16::denormal_min() * 2,
    float16::denormal_min() * 3,
    float16::denormal_min() * 4};

  got = float16::decode (float16::encode (expected));
  for (uint i = 0; i < 4; ++i) {
    EXPECT_NEAR (expected[i], got[i], abs (expected[i] * 0.001f));
  }
  got = float16::decode (
    float16::encode (vec_set<f32_x4> (float16::max() * 4.f)));
  for (uint i = 0; i < 4; ++i) {
    EXPECT_EQ (float16::max(), got[i]);
  }
}
//------------------------------------------------------------------------------
} // namespace artv
