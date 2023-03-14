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
  EXPECT_NEAR (a.to_floatp(), 0.6, 0.0001);
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
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 1, 15>>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_less_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 3, 12>::from_float (-0.7);
  auto c = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_lossy)
{
  auto a = fixpt<1, 3, 12, 0>::from_float (0.6);
  auto b = fixpt<1, 0, 15, 0>::from_float (-0.7);
  auto c = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_rounding)
{
  constexpr auto flags = fixpt_rounding;
  auto           a     = fixpt<1, 4, 11, flags>::from_float (0.6);
  auto           b     = fixpt<1, 4, 11, flags>::from_float (-2.52);
  fixpt<1, 15, 0, flags | fixpt_relaxed_frac_assign> c {};
  c = a + b;
  EXPECT_NEAR (c.to_floatp(), -2., 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_unsigned)
{
  constexpr auto flags = fixpt_dynamic;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<0, 0, 15, flags>::from_float (0.7);
  auto           c     = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), 0.6 + 0.7, 0.00011);
}
//------------------------------------------------------------------------------
TEST (fixed_point, add_gaining_sign)
{
  constexpr auto flags = fixpt_dynamic;
  auto           a     = fixpt<0, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<1, 0, 15, flags>::from_float (-0.7);
  auto           c     = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_equal_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 1, 15>>);
  EXPECT_NEAR (c.to_floatp(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_less_frac_bits)
{
  auto a = fixpt<1, 0, 15>::from_float (0.6);
  auto b = fixpt<1, 3, 12>::from_float (-0.7);
  auto c = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_floatp(), 1.3, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (0.6);
  auto b = fixpt<1, 0, 15>::from_float (-0.7);
  auto c = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15>>);
  EXPECT_NEAR (c.to_floatp(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_lossy)
{
  auto a = fixpt<1, 2, 13, 0>::from_float (1.6);
  auto b = fixpt<1, 5, 10, 0>::from_float (-1.7);
  auto c = a + b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0005);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<1, 0, 15, flags>::from_float (-0.7);
  auto           c     = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_unsigned)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<0, 0, 15, flags>::from_float (0.7);
  auto           c     = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), -0.1, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_gaining_sign)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<0, 3, 12, flags>::from_float (0.6);
  auto           b     = fixpt<1, 0, 15, flags>::from_float (-0.7);
  auto           c     = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.3, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, sub_unsigned_gains_no_bits)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<0, 3, 12, flags>::from_float (0.7);
  auto           b     = fixpt<0, 0, 15, flags>::from_float (0.6);
  auto           c     = a - b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<0, 3, 15, flags>>);
  EXPECT_NEAR (c.to_floatp(), 0.1, 0.0002);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_equal_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 2, 28>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_less_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 3, 12>::from_float (-1.7);
  auto c = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_lossy)
{
  auto a = fixpt<1, 3, 12, 0>::from_float (1.6);
  auto b = fixpt<1, 5, 10, 0>::from_float (-1.7);
  auto c = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_unsigned)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<0, 1, 14, flags>::from_float (1.7);
  auto           c     = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * 1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, mul_gaining_sign)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<0, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a * b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 4, 26, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 * -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_equal_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 15, 15>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_less_frac_bits)
{
  auto a = fixpt<1, 1, 14>::from_float (1.6);
  auto b = fixpt<1, 3, 12>::from_float (-1.7);
  auto c = a / b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 13, 17>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rhs_more_frac_bits)
{
  auto a = fixpt<1, 3, 12>::from_float (1.6);
  auto b = fixpt<1, 1, 14>::from_float (-1.7);
  auto c = a / b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_reciprocal)
{
  auto a = fixpt<1, 1, 0>::from_int (1);
  auto b = fixpt<1, 0, 15>::from_float (-0.9999999999);
  auto c = a / b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), fixpt<1, 16, 0>>);
  EXPECT_NEAR (b.to_floatp(), 1. / -0.9999999999, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_lossy)
{
  auto a = fixpt<1, 3, 12, false>::from_float (1.6);
  auto b = fixpt<1, 1, 14, false>::from_float (-1.7);
  auto c = a / b;
  c.normalize();
  static_assert (std::is_same_v<decltype (c), decltype (a)>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_rounding)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_unsigned)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<1, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<0, 1, 14, flags>::from_float (1.7);
  auto           c     = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / 1.7, 0.001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, div_gaining_sign)
{
  constexpr auto flags = fixpt_dynamic | fixpt_rounding;
  auto           a     = fixpt<0, 3, 12, flags>::from_float (1.6);
  auto           b     = fixpt<1, 1, 14, flags>::from_float (-1.7);
  auto           c     = a / b;
  static_assert (std::is_same_v<decltype (c), fixpt<1, 17, 13, flags>>);
  EXPECT_NEAR (c.to_floatp(), 1.6 / -1.7, 0.001);
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
  auto r1 = a.to_floatp();
  auto r  = a.cast<fixpt<1, 5, 10>>().to_floatp();
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
TEST (fixed_point, fixed_point_pow2_ratio)
{
  auto fp = fixpt_from_ratio<1, 2>();
  static_assert (fp.n_int == 0);
  static_assert (fp.n_frac == 1);
  EXPECT_EQ (static_cast<double> (fp), 0.5);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_decimals_ratio)
{
  auto fp = fixpt_from_ratio<1, 3>();
  static_assert (fp.n_int == 0);
  static_assert (fp.n_frac == 63);
  EXPECT_NEAR (static_cast<double> (fp), .333333333333333333, .00000000001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_pow2_ratio_add)
{
  auto a = fixpt<1, 15, 16>::from_float (23.);
  auto b = a + 0.5_r;
  b.normalize();
  static_assert (b.n_int == 16);
  static_assert (b.n_frac == 16);
  EXPECT_NEAR (b.to_floatp(), 23.5, .00000000001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_infinite_decimals_ratio_add)
{
  auto a = fixpt<1, 15, 16, fixpt_mixed>::from_float (23.);
  auto b = a + (1_r / 3_r);
  b.normalize();
  static_assert (b.n_int == 16);
  static_assert (b.n_frac == 16);
  EXPECT_NEAR (b.to_floatp(), 23.3333, .0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_pow2_ratio_sub)
{
  auto a = fixpt<1, 15, 16>::from_float (23.);
  auto b = a - 0.5_r;
  b.normalize();
  static_assert (b.n_int == 16);
  static_assert (b.n_frac == 16);
  EXPECT_NEAR (b.to_floatp(), 22.5, .00000000001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_infinite_decimals_ratio_sub)
{
  auto a = fixpt<1, 15, 16, fixpt_mixed>::from_float (23.);
  auto b = a - (1_r / 3_r);
  b.normalize();
  static_assert (b.n_int == 16);
  static_assert (b.n_frac == 16);
  EXPECT_NEAR (b.to_floatp(), 22.6666, .0001);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_pow2_ratio_mul)
{
  auto a = fixpt<1, 2, 0>::from_float (2.);
  auto b = a * 0.5_r;
  b.normalize();
  static_assert (b.n_int == 1);
  static_assert (b.n_frac == 1);
  EXPECT_EQ (b.to_floatp(), 1.);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_infinite_decimals_ratio_mul)
{
  using T = fixpt<1, 3, 0, fixpt_mixed>;
  auto a  = T::max();
  auto b  = a * ((1_r / 3_r) & 15_max_ratio_fracb);
  b.normalize();
  static_assert (b.n_int == 2);
  static_assert (b.n_frac == 15);
  EXPECT_NEAR (b.to_floatp(), a.to_floatp() * (1. / 3.), .00015);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_integer_ratio_mul)
{
  using T = fixpt<1, 3, 1, fixpt_mixed>;
  auto a  = T::max();
  // 6 = 110. shifts 3 bits, reduces the fraction 1 bit (trailing 0)
  auto b = a * 6_r;
  b.normalize();
  static_assert (b.n_int == 6);
  static_assert (b.n_frac == 0);
  EXPECT_NEAR (b.to_floatp(), a.to_floatp() * 6., .00015);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_integer_non_pure_ratio_mul)
{
  using T = fixpt<1, 3, 1, fixpt_mixed>;
  auto a  = T::max();
  auto b  = a * ((14_r / 10_r) & 15_max_ratio_fracb);
  b.normalize();
  static_assert (b.n_int == 4);
  static_assert (b.n_frac == 16);
  EXPECT_NEAR (b.to_floatp(), a.to_floatp() * 1.4, .00015);
}
//------------------------------------------------------------------------------
// divisions are implemented as inverse multiplications, not a lot to test...
TEST (fixed_point, fixed_point_pow2_ratio_div)
{
  auto a = fixpt<1, 15, 16>::from_float (2.);
  auto b = a / 0.5_r;
  b.normalize();
  static_assert (b.n_int == 16);
  static_assert (b.n_frac == 15);
  EXPECT_EQ (b.to_floatp(), 4.);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_infinite_decimals_ratio_div)
{
  auto a = fixpt<1, 3, 16, fixpt_mixed>::from_float (3.);
  auto b = a / 3_r;
  b.normalize();
  static_assert (b.n_int == 2);
  static_assert (b.n_frac == 32);
  EXPECT_NEAR (b.to_floatp(), 1., .00015);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_mixed_add_saturation)
{
  constexpr uint intb  = 62;
  constexpr uint fracb = 1;

  auto a = fixpt<1, intb, fracb, fixpt_mixed>::from_float (3.);
  auto b = a + a;
  static_assert (b.n_int == 63);
  static_assert (b.n_frac == 0);
  b.normalize();
  EXPECT_NEAR (b.to_floatp(), 3 + 3, .00015);
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_mixed_sub_saturation)
{
  constexpr uint intb  = 62;
  constexpr uint fracb = 1;

  auto a = fixpt<1, intb, fracb, fixpt_mixed>::from_float (3.);
  auto b = a - fixpt<1, intb, fracb, fixpt_mixed>::from_float (-3.);
  static_assert (b.n_int == 63);
  static_assert (b.n_frac == 0);
  b.normalize();
  EXPECT_NEAR (b.to_floatp(), 3 - -3, .00015);
  // auto c = b + a; // Triggers assert!
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_mixed_mul_saturation)
{
  constexpr uint intb  = 15;
  constexpr uint fracb = 16;

  auto a = fixpt<1, intb, fracb, fixpt_mixed>::from_float (3.);
  auto b = a * a;
  static_assert (b.n_int == (intb * 2));
  static_assert (b.n_frac == (fracb * 2));
  b.normalize();
  EXPECT_NEAR (b.to_floatp(), 3 * 3, .00015);
  auto c = b * a;
  static_assert (c.n_int == (intb * 3));
  static_assert (c.n_frac == (64 - (intb * 3) - 1));
  c.normalize();
  EXPECT_NEAR (c.to_floatp(), 3 * 3 * 3, .00015);
  auto d = c * a;
  static_assert (d.n_int == (intb * 4));
  static_assert (d.n_frac == (64 - (intb * 4) - 1));
  d.normalize();
  EXPECT_NEAR (d.to_floatp(), 3 * 3 * 3 * 3, .00015);
  // auto e = d * a; // Triggers assert!
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_mixed_div_saturation)
{
  constexpr uint intb  = 62;
  constexpr uint fracb = 1;

  s64    aval = (1ull << 40) - 1;
  double bval = 7.1234;

  auto a = fixpt<1, 40, 23, fixpt_mixed>::from_int (aval);
  auto b = fixpt<1, 40, 20, fixpt_mixed>::from_float (bval);
  auto c = a / b;
  static_assert (c.n_int == 60);
  static_assert (c.n_frac == 3);
  c.normalize();
  // errors in the denominator magnify a lot
  auto err = 1. - (b.to_floatp() / bval);
  EXPECT_NEAR (c.to_floatp(), aval / bval, (aval / bval) * err);
  // auto c = b + a; // Triggers assert!
}
//------------------------------------------------------------------------------
TEST (fixed_point, fixed_point_implicit_enabled)
{
  // implicit assignment and arithmetic with float and int types.
  constexpr uint intb  = 62;
  constexpr uint fracb = 1;

  auto a = fixpt<1, 61, 2, fixpt_mixed | fixpt_implicit>::from_float (3.);
  auto b = a + 3;
  static_assert (b.n_int == 62);
  static_assert (b.n_frac == 1);
  b.normalize();
  EXPECT_NEAR (b.to_floatp(), 3 + 3, .00015);
  auto c = b + 3.;
  static_assert (c.n_int == 63);
  static_assert (c.n_frac == 0);
  b.normalize();
  EXPECT_NEAR (c.to_floatp(), 3 + 3 + 3, .00015);
  a = 7.f;
  EXPECT_NEAR (a.to_floatp(), 7, .00015);
  a = 4;
  EXPECT_NEAR (a.to_floatp(), 4, .00015);
}
//------------------------------------------------------------------------------
TEST (fixed_point, rounding_and_truncation)
{
  using base       = fixpt_d<1, 5, 10>;
  using truncating = fixpt_d<1, 5, 0>;
  using rounding   = truncating::rounding_twin;

  static_assert (rounding::rounds_nearest);

  auto v     = base::from_float (1.9);
  auto initv = v.to_floatp();

  EXPECT_EQ (((truncating) v).to_floatp(), 1.);
  EXPECT_EQ (((rounding) v).to_floatp(), 2.);

  v = base::from_float (-1.9);
  EXPECT_EQ (((truncating) v).to_floatp(), -1.);
  EXPECT_EQ (((rounding) v).to_floatp(), -2.);
}
//------------------------------------------------------------------------------
TEST (fixed_point, integer_operators)
{
  auto v = fixpt_d<1, 7, 8>::from_float (1.5f);

  EXPECT_TRUE (v != v.max());
  EXPECT_TRUE (v < v.max());
  EXPECT_TRUE (v > v.min());
  EXPECT_TRUE (v == v);
  EXPECT_TRUE (v <= v);
  EXPECT_TRUE (v >= v);
  EXPECT_TRUE (v >= v);
  auto w = v.resize<1, -1>();
  EXPECT_TRUE (v == w);
  EXPECT_TRUE (v <= w);
  EXPECT_TRUE (v >= w);
  EXPECT_TRUE (v >= w);
}

//------------------------------------------------------------------------------
} // namespace artv
