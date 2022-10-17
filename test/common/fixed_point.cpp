#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gtest/gtest.h>

#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
TEST (fixed_point, float_set_and_cast)
{
  auto a = fp_int<1, 0, 15>::from_float (0.6);
  EXPECT_NEAR (a.as_float(), 0.6, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_equal_frac_bits)
{
  auto a = fp_int<1, 0, 15>::from_float (0.6);
  auto b = a.from_float (-0.7f);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 1, 15>>);
  EXPECT_NEAR (c.as_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_less_frac_bits)
{
  auto a = fp_int<1, 0, 15>::from_float (0.6);
  auto b = fp_int<1, 3, 12>::from_float (-0.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 15>>);
  EXPECT_NEAR (c.as_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_more_frac_bits)
{
  auto a = fp_int<1, 3, 12>::from_float (0.6);
  auto b = fp_int<1, 0, 15>::from_float (-0.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 15>>);
  EXPECT_NEAR (c.as_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_equal_frac_bits)
{
  auto a = fp_int<1, 0, 15>::from_float (0.6);
  auto b = fp_int<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 1, 15>>);
  EXPECT_NEAR (c.as_float(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_less_frac_bits)
{
  auto a = fp_int<1, 0, 15>::from_float (0.6);
  auto b = fp_int<1, 3, 12>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 15>>);
  EXPECT_NEAR (c.as_float(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_more_frac_bits)
{
  auto a = fp_int<1, 3, 12>::from_float (0.6);
  auto b = fp_int<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 15>>);
  EXPECT_NEAR (c.as_float(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_equal_frac_bits)
{
  auto a = fp_int<1, 1, 14>::from_float (1.6);
  auto b = fp_int<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 2, 28>>);
  EXPECT_NEAR (c.as_float(), 1.6 * -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_less_frac_bits)
{
  auto a = fp_int<1, 1, 14>::from_float (1.6);
  auto b = fp_int<1, 3, 12>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 26>>);
  EXPECT_NEAR (c.as_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_more_frac_bits)
{
  auto a = fp_int<1, 3, 12>::from_float (1.6);
  auto b = fp_int<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 4, 26>>);
  EXPECT_NEAR (c.as_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_equal_frac_bits)
{
  auto a = fp_int<1, 1, 14>::from_float (1.6);
  auto b = fp_int<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 15, 15>>);
  EXPECT_NEAR (c.as_float(), 1.6 / -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_less_frac_bits)
{
  auto a = fp_int<1, 1, 14>::from_float (1.6);
  auto b = fp_int<1, 3, 12>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 13, 17>>);
  EXPECT_NEAR (c.as_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_more_frac_bits)
{
  auto a = fp_int<1, 3, 12>::from_float (1.6);
  auto b = fp_int<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 17, 13>>);
  EXPECT_NEAR (c.as_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_reciprocal)
{
  auto a = fp_int<1, 1, 0>::from_int (1);
  auto b = fp_int<1, 0, 15>::from_float (-0.9999999999);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fp_int<1, 16, 0>>);
  EXPECT_NEAR (b.as_float(), 1. / -0.9999999999, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, normalize)
{
  using fp = fp_int<1, 0, 4>;

  s8   v = 0x7f;
  auto a = fp::from (v);
  a.normalize();
  EXPECT_EQ (a.value(), (s8) 0x0f); // 15

  v = 0x80;
  a = fp::from (v);
  a.normalize();
  EXPECT_EQ (a.value(), (s8) 0xf0); // -16
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_cast)
{
  auto a  = fp_int<1, 15, 16>::from_float (23.1123232456);
  auto r1 = a.as_float();
  auto r  = a.cast<fp_int<1, 5, 10>>().as_float();
  EXPECT_NEAR (r, 23.1123232456, 0.001);
}
//------------------------------------------------------------------------------

} // namespace artv
