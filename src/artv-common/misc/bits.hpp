#pragma once

#include <cassert>
#include <limits>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <
  class T,
  std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>>* = nullptr>
constexpr T lsb_mask (uint bits)
{
  constexpr T ones      = (T) ~0ull;
  constexpr T type_bits = sizeof (T) * 8;

  assert (bits <= type_bits);
  return ones >> (type_bits - bits);
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
constexpr T bit (uint bit, T value = 1)
{
  assert (bit < sizeof (T) * 8);
  return value << bit;
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
constexpr void set_bit (T& v, uint bit_idx, bool on = true)
{
  assert (bit_idx < sizeof (T) * 8);

  auto clear_mask = bit<T> (bit_idx);
  auto set_mask   = bit<T> (bit_idx, on);
  v &= ~clear_mask;
  v |= set_mask;
}
//------------------------------------------------------------------------------
namespace detail {
template <class T, class F, class FL, class FLL>
T unary_int_func_dispatch (T v, F&& f, FL&& fl, FLL&& fll)
{
  static_assert (std::is_integral_v<T>);
  if constexpr (sizeof v <= sizeof (int)) {
    return f ((unsigned int) v);
  }
  else if constexpr (sizeof v <= sizeof (long)) {
    return fl ((unsigned long) v);
  }
  else if constexpr (sizeof v <= sizeof (long long)) {
    return fll ((unsigned long long) v);
  }
  else {
    static_assert (sizeof v == 0, "unsupported type");
    return T {};
  }
}
} // namespace detail
//------------------------------------------------------------------------------
// 1-based indexes
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline T first_bit_set (T v)
{
  return detail::unary_int_func_dispatch (
    v,
    [] (auto v) { return __builtin_ffs (v); },
    [] (auto v) { return __builtin_ffsl (v); },
    [] (auto v) { return __builtin_ffsll (v); });
}
//------------------------------------------------------------------------------
// 1-based indexes
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline T last_bit_set (T v)
{
  return detail::unary_int_func_dispatch (
    v,
    [] (auto v) { return v ? (sizeof v * 8) - __builtin_clz (v) : 0; },
    [] (auto v) { return v ? (sizeof v * 8) - __builtin_clzl (v) : 0; },
    [] (auto v) { return v ? (sizeof v * 8) - __builtin_clzll (v) : 0; });
}
//------------------------------------------------------------------------------
template <class T, class F, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
void iterate_set_bits (T v, F f)
{
  unsigned pos = 0;
  while (unsigned curr_pos = first_bit_set (v)) {
    pos += curr_pos;
    v >>= curr_pos;
    f (pos - 1);
  }
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline auto popcount (T v)
{
  return detail::unary_int_func_dispatch (
    v,
    [] (auto v) { return __builtin_popcount (v); },
    [] (auto v) { return __builtin_popcountl (v); },
    [] (auto v) { return __builtin_popcountll (v); });
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr bool is_pow2 (T v)
{
  return (v != 0) && ((v & (v - 1)) == 0);
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline T pow2_round_floor (T v)
{
  uint lsb = last_bit_set (v);
  return lsb ? bit<T> (lsb - 1) : v;
}
//------------------------------------------------------------------------------
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline T pow2_round_ceil (T v)
{
  using overload_type = std::conditional_t<sizeof (T) == 8, u64, u32>;
  uint lsb            = last_bit_set ((overload_type) v);
  return (is_pow2 (v) || !v) ? v : bit<T> (lsb);
}
//------------------------------------------------------------------------------
// just to avoid undefined and implementation defined behavior on signed
// integer shiftings at the expense of probably (slightly) worse code. Notice
// that it is assumed:
//
// - Obvious: the optimizers can see multiplications and divisions by powers
//   of two as shifts.
// - For unsigned types the optimizer is able to deduce to a plain shift, so
//   not bothering writing versions for both signedness types.
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift Right
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ashr (T v, uint n)
{

  // might be a vector type (hence the weird write to write the shifts
  // below)
  assert (n < (sizeof (T) * 8));
  return v / ((T {} + 1) << (T {} + n));
}

template <uint N, class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ashr (T v)
{
  static_assert (N < (sizeof (T) * 8)); // just making the assertion comptime
  return ashr (v, N);
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift Left
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ashl (T v, uint n)
{
  assert (n < (sizeof (T) * 8));
  return v * ((T {} + 1) << (T {} + n));
}

template <uint N, class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ashl (T v)
{
  static_assert (N < (sizeof (T) * 8));
  return ashl (v, N);
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift. Positive values increase the
// result (shift left). Negative ones decrease it (shift right).
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ash (T v, int n)
{
  return n >= 0 ? ashl (v, n) : ashr (v, -n);
}

template <int N, class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T ash (T v)
{
  if constexpr (N >= 0) {
    return ashl<N> (v);
  }
  else {
    return ashr<-N> (v);
  }
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) LSB mask
template <class T, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline constexpr T alsb_mask (T v, uint bits)
{
  using T_u = std::make_unsigned_t<T>;

  auto const mask = lsb_mask<T_u> (bits);
  if constexpr (std::is_unsigned_v<T>) {
    return v & mask;
  }
  else {
    if (v >= 0) {
      // masking a positive is filling the MSBs with 0's
      return (T) ((T_u) v & mask);
    }
    else {
      // masking a negative is filling the MSBs with 1's
      return (T) ((T_u) v | ~mask);
    }
  }
}
//------------------------------------------------------------------------------

} // namespace artv
