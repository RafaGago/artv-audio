#pragma once

#include <cassert>
#include <limits>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <class T>
constexpr T lsb_mask (uint bits)
{
  static_assert (std::is_integral_v<T> && std::is_unsigned_v<T>, "");
  constexpr T ones      = (T) ~0ull;
  constexpr T type_bits = sizeof (T) * 8;

  assert (bits <= type_bits);
  return ones >> (type_bits - bits);
}
//------------------------------------------------------------------------------
template <class T>
constexpr T bit (uint bit, T value = 1)
{
  static_assert (std::is_integral_v<T>, "");
  assert (bit < sizeof (T) * 8);
  return value << bit;
}
//------------------------------------------------------------------------------
template <class T>
constexpr void set_bit (T& v, uint bit_idx, bool on = true)
{
  static_assert (std::is_integral_v<T>, "");
  assert (bit_idx < sizeof (T) * 8);

  auto clear_mask = bit<T> (bit_idx);
  auto set_mask   = bit<T> (bit_idx, on);
  v &= ~clear_mask;
  v |= set_mask;
}
//------------------------------------------------------------------------------
#if defined(__GNUC__)

#if 0 // quick and dirty. don't feel like doing proper integer type
      // (clusterfuck) type wrapping now.
inline unsigned first_bit_set (int v)
{
    return __builtin_ffs (v);
}

inline unsigned first_bit_set (long v)
{
    return __builtin_ffsl (v);
}

inline unsigned last_bit_set (int v)
{
  return v ? (sizeof v * 8) - __builtin_clz (v) : 0;
}

inline unsigned last_bit_set (long v)
{
    return v ? (sizeof v * 8) - __builtin_clzl (v) : 0;
}

#endif
//------------------------------------------------------------------------------
// returns 1 based indexes, 0 for no bits set
inline unsigned first_bit_set (long long v)
{
  return __builtin_ffsll (v);
}
//------------------------------------------------------------------------------
// returns 1 based indexes, 0 for no bits set
inline unsigned last_bit_set (long long v)
{
  return v ? (sizeof v * 8) - __builtin_clzll (v) : 0;
}
//------------------------------------------------------------------------------
#elif defined(_MSC_VER)

} // namespace artv

#include <intrin.h>

namespace artv {

// returns 1 based indexes, 0 for no bits set
inline unsigned first_bit_set (u32 v)
{
  unsigned long index;
  char          r = _BitScanForward (&index, v);
  return r ? index + 1 : 0;
}

inline unsigned first_bit_set (u64 v)
{
  unsigned long index;
  char          r = _BitScanForward (&index, v);
  return r ? index + 1 : 0;
}

// returns 1 based indexes, 0 for no bits set
inline unsigned last_bit_set (u32 v)
{
  unsigned long index;
  char          r = _BitScanReverse (&index, v);
  return r ? index + 1 : 0;
}

inline unsigned last_bit_set (u64 v)
{
  unsigned long index;
  char          r = _BitScanReverse (&index, v);
  return r ? index + 1 : 0;
}

#else
static_assert (false, "TBI");
#endif
//------------------------------------------------------------------------------
template <class T, class F>
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
#if 0 // unused, would need a MSVC/Clang-cl impl
// TODO: set -march
#if defined(__GNUC__)
inline int popcount (uint x)
{
  return __builtin_popcount (x);
}
#else
static_assert (false, "To be implemented");
#endif
#endif
//------------------------------------------------------------------------------
template <class T>
inline constexpr bool is_pow2 (T v)
{
  static_assert (std::is_integral_v<T>, "");
  return (v != 0) && ((v & (v - 1)) == 0);
}
//------------------------------------------------------------------------------
template <class T>
inline T pow2_round_floor (T v)
{
  static_assert (std::is_integral_v<T>, "");
  uint lsb = last_bit_set (v);
  return lsb ? bit<T> (lsb - 1) : v;
}
//------------------------------------------------------------------------------
template <class T>
inline T pow2_round_ceil (T v)
{
  using overload_type = std::conditional_t<sizeof (T) == 8, u64, u32>;
  static_assert (std::is_integral_v<T>, "");
  uint lsb = last_bit_set ((overload_type) v);
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

// Arithmetic (complement 2 respecting) Shift Right
template <class T>
inline constexpr T ashr (T v, uint n)
{
  if constexpr (std::is_integral_v<T>) {
    // might be a vector type (hence the weird write to write the shifts below)
    assert (n < (sizeof (T) * 8));
  }
  return v / ((T {} + 1) << (T {} + n));
}

template <uint N, class T>
inline constexpr T ashr (T v)
{
  if constexpr (std::is_integral_v<T>) {
    // might be a vector type
    static_assert (N < (sizeof (T) * 8)); // just making the assertion comptime
  }
  return ashr (v, N);
}
// Arithmetic (complement 2 respecting) Shift Left
template <class T>
inline constexpr T ashl (T v, uint n)
{
  if constexpr (std::is_integral_v<T>) {
    // might be a vector type (hence the weird write to write the shifts below)
    assert (n < (sizeof (T) * 8));
  }
  return v * ((T {} + 1) << (T {} + n));
}

template <uint N, class T>
inline constexpr T ashl (T v)
{
  if constexpr (std::is_integral_v<T>) {
    // might be a vector type
    static_assert (N < (sizeof (T) * 8)); // just making the assertion comptime
  }
  return ashl (v, N);
}
// Arithmetic (complement 2 respecting) Shift. Positive values increase the
// result (shift left). Negative ones decrease it (shift right).
template <class T>
inline constexpr T ash (T v, int n)
{
  return v >= 0 ? ashl (v, n) : ashr (v, -n);
}

template <int N, class T>
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
template <class T>
inline constexpr T alsb_mask (T v, uint bits)
{
  static_assert (std::is_integral_v<T>);
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

template <class T, uint size_of>
inline constexpr T alsb_mask (
  T __attribute__ ((vector_size (size_of))) v,
  uint bits)
{
  static_assert (std::is_integral_v<T>);
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
