#pragma once

#include <cstdint>
#include <type_traits>

namespace artv {

using s8  = int8_t;
using u8  = uint8_t;
using s16 = int16_t;
using u16 = uint16_t;
using s32 = int32_t;
using u32 = uint32_t;
using s64 = int64_t;
using u64 = uint64_t;
// clang and GCC only!
using s128 = __int128;
using u128 = unsigned __int128;

using uint   = unsigned int;
using ulong  = unsigned long;
using u2long = unsigned long long;
using sint   = signed int;
using slong  = signed long;
using s2long = signed long long;

namespace detail {

template <uint size>
struct int_for_size;

template <uint size>
struct uint_for_size;

// clang-format off
template <> struct int_for_size<1> { using type = s8; };
template <> struct int_for_size<2> { using type = s16; };
template <> struct int_for_size<4> { using type = s32; };
template <> struct int_for_size<8> { using type = s64; };
template <> struct int_for_size<16> { using type = s128; };

template <> struct uint_for_size<1> { using type = u8; };
template <> struct uint_for_size<2> { using type = u16; };
template <> struct uint_for_size<4> { using type = u32; };
template <> struct uint_for_size<8> { using type = u64; };
template <> struct uint_for_size<16> { using type = u128; };
// clang-format on

} // namespace detail

template <uint size>
using int_for_size = typename detail::int_for_size<size>::type;

template <uint size>
using uint_for_size = typename detail::uint_for_size<size>::type;

template <class T>
using same_size_int = int_for_size<sizeof (T)>;

template <class T>
using same_size_uint = uint_for_size<sizeof (T)>;

} // namespace artv
