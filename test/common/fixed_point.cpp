#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gtest/gtest.h>

#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/vec_math.hpp"

namespace artv {

//------------------------------------------------------------------------------
TEST (fixed_point, float_set_and_cast)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  EXPECT_NEAR (a.to_float(), 0.6, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, limits)
{
  using T = fixpt<1, 0, 15>;
  EXPECT_EQ (T::min_raw(), std::numeric_limits<s16>::min());
  EXPECT_EQ (T::max_raw(), std::numeric_limits<s16>::max());
  EXPECT_EQ (T::min_int(), -1);
  EXPECT_EQ (T::max_int(), 0);
  EXPECT_NEAR (T::min_float(), -1., 0.0001);
  EXPECT_NEAR (T::max_float(), 1., 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_equal_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = a.from_float (-0.7f);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 1, 15>>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_less_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 3, 12>::from_float (-0.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_lossy)
{
  auto a = fixpt<1, 3, 12, 0>::from_float (0.6);
  auto b = fixpt<1, 0, 15, 0>::from_float (-0.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<1, 0, 15, flags>::from_float (-0.7);
  auto           c     = a + b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_equal_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 1, 15>>);
  EXPECT_NEAR (c.to_float(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_less_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 3, 12>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_float(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_float(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_lossy)
{
  auto a = fixpt<1, 2, 13, 0>::from_float (1.6);
  auto b = fixpt<1, 5, 10, 0>::from_float (-1.7);
  auto c = a + b;
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_float(), -0.1, 0.0005);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<1, 0, 15, flags>::from_float (-0.7);
  auto           c     = a - b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_float(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_equal_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 2, 28>>);
  EXPECT_NEAR (c.to_float(), 1.6 * -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_less_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 3, 12>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26>>);
  EXPECT_NEAR (c.to_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26>>);
  EXPECT_NEAR (c.to_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_lossy)
{
  auto a = fixpt<1, 3, 12, 0>::from_float (1.6);
  auto b = fixpt<1, 5, 10, 0>::from_float (-1.7);
  auto c = a * b;
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a * b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26, flags>>);
  EXPECT_NEAR (c.to_float(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_equal_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 15, 15>>);
  EXPECT_NEAR (c.to_float(), 1.6 / -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_less_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 3, 12>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 13, 17>>);
  EXPECT_NEAR (c.to_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13>>);
  EXPECT_NEAR (c.to_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_reciprocal)
{
  auto a = fixpt<1, 1, 0>::from_int (1);
  auto b = fixpt<1, 0, 15>::from_float (-0.9999999999);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 16, 0>>);
  EXPECT_NEAR (b.to_float(), 1. / -0.9999999999, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_lossy)
{
  auto a = fixpt<1, 3, 12, false>::from_float (1.6);
  auto b = fixpt<1, 1, 14, false>::from_float (-1.7);
  auto c = a / b;
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13, flags>>);
  EXPECT_NEAR (c.to_float(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, normalize)
{
  using fp = fixpt<1, 0, 4>;

  s8   v = 0x7f;
  auto a = fp::from (v);
  a.normalize();
  EXPECT_EQ (a.value(), (s8) 0x1f); // 31

  v = 0x80;
  a = fp::from (v);
  a.normalize();
  EXPECT_EQ (a.value(), (s8) 0xe0); // -32
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_cast)
{
  auto a  = fixpt<1, 15, 16>::from_float (23.1123232456);
  auto r1 = a.to_float();
  auto r  = a.cast<fixpt<1, 5, 10>>().to_float();
  EXPECT_NEAR (r, 23.1123232456, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_cast_w_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 0, 7, flags>::from (0x05); // 0000 0101
  auto           b     = a.resize<0, -1>(); // drop 1 fractional bit
  EXPECT_EQ (b.value(), 0x03); //                               _000 0011
}
//------------------------------------------------------------------------------

} // namespace artv
