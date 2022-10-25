#pragma once

#include <cassert>
#include <limits>
#include <type_traits>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_util.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <class V, enable_if_uint_vec_t<V>* = nullptr>
constexpr inline V lsb_mask (V bits)
{
  using T               = vec_value_type_t<V>;
  constexpr T ones      = (T) ~0ull;
  constexpr T type_bits = sizeof (T) * 8;

  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    assert (bits[i] <= type_bits);
  }
  return ones >> (type_bits - bits);
}

template <class V, enable_if_uint_vec_t<V>* = nullptr>
constexpr inline V lsb_mask (uint bits)
{
  return lsb_mask<V> (vec_set<V> (bits));
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
constexpr inline V bit (V bit, V value = vec_set<V> (1))
{
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    assert (bit[i] < (sizeof (value[0]) * 8));
  }
  return value << bit;
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
constexpr inline V bit (V bitidx, uint value = 1)
{
  return bit<V> (bitidx, vec_set<V> (value));
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
constexpr inline void set_bit (V& v, V bit_idx, V on = vec_set<V> (1))
{
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    assert (bit_idx[i] < sizeof (v[0]) * 8);
  }
  auto clear_mask = bit<V> (bit_idx);
  auto set_mask   = bit<V> (bit_idx, on);
  v &= ~clear_mask;
  v |= set_mask;
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
constexpr inline void set_bit (V& v, V bit_idx, bool on = true)
{
  return set_bit<V> (v, bit_idx, vec_set<V> (on));
}
//------------------------------------------------------------------------------
// TODO : better implementation
// 1-based indexes
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline V first_bit_set (V v)
{
  V r;
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    r[0] = first_bit_set (v[0]);
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline V last_bit_set (V v)
{
  V r;
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    r[0] = last_bit_set (v[0]);
  }
  return r;
}
//------------------------------------------------------------------------------
// "iterate_set_bits" TODO
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline V popcount (V v)
{
  V r;
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    r[0] = popcount (v[0]);
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
constexpr inline bool is_pow2 (V v)
{
  return (v != 0) && ((v & (v - 1)) == 0);
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline V pow2_round_floor (V v)
{
  V lsb = last_bit_set (v);
  return lsb ? bit<V> (lsb - 1) : v;
}
//------------------------------------------------------------------------------
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline V pow2_round_ceil (V v)
{
  V lsb = last_bit_set (v);
  return (is_pow2 (v) || !v) ? v : bit<V> (lsb);
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift Right
template <
  class V,
  class U,
  enable_if_any_int_vec_t<V>* = nullptr,
  enable_if_any_int_vec_t<U>* = nullptr>
inline constexpr V ashr (V v, U n)
{
  static_assert (vec_traits_t<U>::size == vec_traits_t<V>::size);
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    // might be a vector type (hence the weird write to write the shifts
    // below)
    assert (n[i] < (sizeof (v[0]) * 8));
    assert (n[i] >= 0);
  }
  return v / (vec_set<V> (1) << n);
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ashr (V v, uint n)
{
  return ashr<V> (v, vec_set<V> (n));
}

template <uint N, class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ashr (V v)
{
  static_assert (N < (sizeof (v[0]) * 8));
  return ashr (v, N);
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift left
template <
  class V,
  class U,
  enable_if_any_int_vec_t<V>* = nullptr,
  enable_if_any_int_vec_t<U>* = nullptr>
inline constexpr V ashl (V v, U n)
{
  static_assert (vec_traits_t<U>::size == vec_traits_t<V>::size);
  for (uint i = 0; i < vec_traits_t<V>::size; ++i) {
    // might be a vector type (hence the weird write to write the shifts
    // below)
    assert (n[i] < (sizeof (v[0]) * 8));
    assert (n[i] >= 0);
  }
  return v * (vec_set<V> (1) << n);
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ashl (V v, uint n)
{
  return ashl<V> (v, vec_set<V> (n));
}

template <uint N, class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ashl (V v)
{
  static_assert (N < (sizeof (v[0]) * 8));
  return ashr (v, N);
}
//------------------------------------------------------------------------------
// Arithmetic (complement 2 respecting) Shift. Positive values increase the
// result (shift left). Negative ones decrease it (shift right).
template <
  class V,
  class U,
  enable_if_any_int_vec_t<V>* = nullptr,
  enable_if_any_int_vec_t<U>* = nullptr>
inline constexpr V ash (V v, U n)
{
  static_assert (std::is_signed_v<vec_value_type_t<U>>);
  static_assert (vec_traits_t<U>::size == vec_traits_t<V>::size);

  return n >= 0 ? ashl (v, n) : ashr (v, -n);
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ash (V v, int n)
{
  return v >= 0 ? ashl (v, n) : ashr (v, -n);
}

template <int N, class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V ash (V v)
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
template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V alsb_mask (V v, V bits)
{
  using T  = vec_traits_t<V>;
  using Tu = std::make_unsigned_t<T>;
  using Vu = vec<Tu, vec_traits_t<V>::size>;

  auto const mask = lsb_mask<Vu> (vec_cast<Tu> (bits));
  if constexpr (std::is_unsigned_v<T>) {
    return v & mask;
  }
  else {
    if (v >= 0) {
      // masking a positive is filling the MSBs with 0's
      return vec_cast<T> ((vec_cast<Tu> (v) & mask));
    }
    else {
      // masking a negative is filling the MSBs with 1's
      return vec_cast<T> ((vec_cast<Tu> (v) | ~mask));
    }
  }
}

template <class V, enable_if_any_int_vec_t<V>* = nullptr>
inline constexpr V alsb_mask (V v, uint bits)
{
  return alsb_mask (v, vec_set<V> (bits));
}
//------------------------------------------------------------------------------

} // namespace artv
