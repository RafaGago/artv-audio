#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gtest/gtest.h>

#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/vec_approx_math.hpp"
#include "artv-common/misc/vec_math.hpp"

namespace artv {

template <class F1, class F2>
static void simd_check (
  F1     cmath_func,
  F2     vec_func,
  double range_start = 0.,
  double range_end   = 10.,
  double range_step  = 0.01,
  double epsilon     = 0.01) // 1%
{
  double min = 1. - epsilon;
  double max = 1. + epsilon;

  for (double i = range_start; i < range_end; i += range_step) {
    double cm  = cmath_func (i);
    double vec = vec_func (vec_set<f64_x2> (i))[0];
    double fac = cm / vec;

    if ((fac < min || fac > max) && !(cm == 0. && vec == 0.)) {
      printf ("input: %f, std: %f, vec: %f\n", i, cm, vec);
      ASSERT_EQ (cm, vec); // error!
    }
  }
}

//------------------------------------------------------------------------------
TEST (simd, sanity_check)
{
  vec<int, 4> a {0, 1, 2, 3};
  vec<int, 4> b {3, 2, 1, 0};
  auto        r = a + b;
  ASSERT_EQ (r[0], 3);
  ASSERT_EQ (r[1], 3);
  ASSERT_EQ (r[2], 3);
  ASSERT_EQ (r[3], 3);
}
//------------------------------------------------------------------------------
TEST (simd, broadcast)
{
  using T = vec<int, 4>;
  auto r  = vec_set<T> (3);
  ASSERT_EQ (r[0], 3);
  ASSERT_EQ (r[1], 3);
  ASSERT_EQ (r[2], 3);
  ASSERT_EQ (r[3], 3);
}
//------------------------------------------------------------------------------
TEST (simd, unaligned_load_store)
{
  using array_type = std::array<u64, 2>;
  using vec_type   = vec<u64, 2>;

  // enforce unalignment of the src and dst memory blocks.
  u8   pool[(sizeof (array_type) * 2) + sizeof (array_type::value_type) - 1];
  uint memstart = 0;

  constexpr uint mask = sizeof (array_type) - 1;
  while ((reinterpret_cast<uintptr_t> (&pool[memstart]) & mask)
         != sizeof (array_type::value_type)) {
    ++memstart;
  }
  array_type* src = new (&pool[memstart]) array_type {};
  array_type* dst = new (&pool[memstart + sizeof *src]) array_type {};

  for (auto& v : *src) {
    v = (array_type::value_type) rand();
  }
  auto a = vec_load_unaligned<vec_type> (*src);
  vec_store_unaligned (*dst, a);

  ASSERT_EQ (*src, *dst);
}
//------------------------------------------------------------------------------
TEST (simd, aligned_load_store)
{
  using array_type = std::array<u64, 2>;
  using vec_type   = vec<u64, 2>;

  alignas (vec_type) array_type src, dst;
  for (auto& v : src) {
    v = (array_type::value_type) rand();
  }
  auto a = vec_load<vec_type> (src);
  vec_store (dst, a);

  ASSERT_EQ (src, dst);
}
//------------------------------------------------------------------------------
TEST (simd, cast)
{
  vec<float, 4> src = {1.4, 2.4, 3.4, 4.4};
  auto          dst = vec_cast<int> (src);
  ASSERT_EQ (dst[0], 1);
  ASSERT_EQ (dst[1], 2);
  ASSERT_EQ (dst[2], 3);
  ASSERT_EQ (dst[3], 4);
}
//------------------------------------------------------------------------------
TEST (simd, shuffle)
{
  vec<int, 4> src  = {1, 2, 3, 4};
  vec<int, 4> src2 = {5, 6, 7, 8};

  auto dst = vec_shuffle (src, src2, 0, 1, 4, 5);

  ASSERT_EQ (dst[0], 1);
  ASSERT_EQ (dst[1], 2);
  ASSERT_EQ (dst[2], 5);
  ASSERT_EQ (dst[3], 6);
}
//------------------------------------------------------------------------------
TEST (simd, exp)
{
  simd_check (
    [] (auto v) { return exp (v); }, [] (auto v) { return vec_exp (v); });
}
//------------------------------------------------------------------------------
TEST (simd, log)
{
  simd_check (
    [] (auto v) { return log (v); }, [] (auto v) { return vec_log (v); });
}
//------------------------------------------------------------------------------
TEST (simd, log2)
{
  simd_check (
    [] (auto v) { return log2 (v); }, [] (auto v) { return vec_log2 (v); });
}
//------------------------------------------------------------------------------
TEST (simd, exp2)
{
  simd_check (
    [] (auto v) { return exp2 (v); }, [] (auto v) { return vec_exp2 (v); });
}
//------------------------------------------------------------------------------
TEST (simd, pow)
{
  simd_check (
    [] (auto v) { return pow (v, 2.); },
    [] (auto v) { return vec_pow (v, 2.); });
}
//------------------------------------------------------------------------------
TEST (simd, sin)
{
  simd_check (
    [] (auto v) { return sin (v); }, [] (auto v) { return vec_sin (v); });
}
//------------------------------------------------------------------------------
TEST (simd, cos)
{
  simd_check (
    [] (auto v) { return cos (v); }, [] (auto v) { return vec_cos (v); });
}
//------------------------------------------------------------------------------
TEST (simd, tan)
{
  simd_check (
    [] (auto v) { return tan (v); }, [] (auto v) { return vec_tan (v); });
}
//------------------------------------------------------------------------------
TEST (simd, sinh)
{
  simd_check (
    [] (auto v) { return sinh (v); }, [] (auto v) { return vec_sinh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, cosh)
{
  simd_check (
    [] (auto v) { return cosh (v); }, [] (auto v) { return vec_cosh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, tanh)
{
  simd_check (
    [] (auto v) { return tanh (v); }, [] (auto v) { return vec_tanh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, asin)
{
  simd_check (
    [] (auto v) { return asin (v); }, [] (auto v) { return vec_asin (v); });
}
//------------------------------------------------------------------------------
TEST (simd, acos)
{
  simd_check (
    [] (auto v) { return acos (v); }, [] (auto v) { return vec_acos (v); });
}
//------------------------------------------------------------------------------
TEST (simd, atan)
{
  simd_check (
    [] (auto v) { return atan (v); }, [] (auto v) { return vec_atan (v); });
}
//------------------------------------------------------------------------------
TEST (simd, atan2)
{
  simd_check (
    [] (auto v) { return atan2 (v, v); },
    [] (auto v) { return vec_atan2 (v, v); });
}
//------------------------------------------------------------------------------
TEST (simd, asinh)
{
  simd_check (
    [] (auto v) { return sinh (v); }, [] (auto v) { return vec_sinh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, acosh)
{
  simd_check (
    [] (auto v) { return acosh (v); }, [] (auto v) { return vec_acosh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, atanh)
{
  simd_check (
    [] (auto v) { return atanh (v); }, [] (auto v) { return vec_atanh (v); });
}
//------------------------------------------------------------------------------
TEST (simd, sqrt)
{
  simd_check (
    [] (auto v) { return sqrt (v); }, [] (auto v) { return vec_sqrt (v); });
}
//------------------------------------------------------------------------------
TEST (simd, all_ones)
{
  u64_x2 v {1, 3};
  ASSERT_FALSE (vec_is_all_ones (v > 2));
  ASSERT_FALSE (vec_is_all_ones (v > 3));
  ASSERT_TRUE (vec_is_all_ones (v > 0));
}
//------------------------------------------------------------------------------
TEST (simd, all_zeros)
{
  u64_x2 v {1, 3};
  ASSERT_FALSE (vec_is_all_zeros (v > 2));
  ASSERT_FALSE (vec_is_all_zeros (v > 0));
  ASSERT_TRUE (vec_is_all_zeros (v > 3));
}
//------------------------------------------------------------------------------
TEST (simd, approx_tan_prewarp)
{
  for (uint i = 1; i < 49; ++i) {
    float  x  = 0.01 * i;
    f32_x4 r1 = vec_tan_pix_prewarp_2dat (vec_set<f32_x4> (x));
    float  r2 = std::tan (x * M_PI);
    EXPECT_NEAR (r1[0], r2, r2 / 1000000.);
  }
}
//------------------------------------------------------------------------------
TEST (simd, approx_exp)
{
  for (int i = -100; i < 110; ++i) {
    float  x  = 0.1f * i;
    f32_x4 r1 = vec_exp_2dat (vec_set<f32_x4> (x));
    float  r2 = std::exp (x);
    EXPECT_NEAR (r1[0], r2, r2 / 1000000.);
  }
}
//------------------------------------------------------------------------------
TEST (simd, approx_log)
{
  for (int i = 1; i < 11; ++i) {
    float  x  = i;
    f32_x4 r1 = vec_log_2dat (vec_set<f32_x4> (x));
    float  r2 = std::log (x);
    EXPECT_NEAR (r1[0], r2, r2 / 1000000.);
  }
}
//------------------------------------------------------------------------------
TEST (simd, approx_pow)
{
  float y = 0.29f;
  for (int i = 1; i < 11; ++i) {
    float  x  = i;
    f32_x4 r1 = vec_pow_2dat (vec_set<f32_x4> (x), vec_set<f32_x4> (y));
    float  r2 = std::pow (x, y);
    EXPECT_NEAR (r1[0], r2, r2 / 1000000.);
  }
  float x = 2;
  for (int i = 1; i < 11; ++i) {
    y         = i;
    f32_x4 r1 = vec_pow_2dat (vec_set<f32_x4> (x), vec_set<f32_x4> (y));
    float  r2 = std::pow (x, y);
    EXPECT_NEAR (r1[0], r2, r2 / 1000000.);
  }
}
//------------------------------------------------------------------------------
TEST (simd, approx_tanh)
{
  for (int i = -100; i < 110; ++i) {
    float  x  = 0.1f * i;
    f32_x4 r1 = vec_tanh_ln2_div2_2dat (vec_set<f32_x4> (x));
    float  r2 = std::tanh (x * M_LN2 * 0.5);
    EXPECT_NEAR (r1[0], r2, abs (r2 / 1000000.));
  }
}
//------------------------------------------------------------------------------
} // namespace artv
