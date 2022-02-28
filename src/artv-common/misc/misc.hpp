#pragma once

// miscellanious things that probably don't deserve its own header.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include <utility>
#include <vector>

#include <gcem.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

// TODO: Split this.

namespace artv {

#ifdef __GNUC__
#define likely(x) __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)
#else
#define likely(x) x
#define unlikely(x) x
#endif

#define VERSION_GET(major, minor, rev) (major * 1000000 + minor * 1000 + rev)

#define PP_EXPAND(x) x
#define PP_TO_STR(x) #x

// usage error_with_template_type<whatever> v; The compiler will print the type.
template <class... Ts>
struct error_with_template_type;
//------------------------------------------------------------------------------
#define array_elems(x) sizeof (x) / sizeof (x[0])

template <class T>
struct is_std_array : public std::false_type {};

template <class T, size_t N>
struct is_std_array<std::array<T, N>> : public std::true_type {};

template <class T>
static constexpr bool is_std_array_v = is_std_array<T>::value;
//------------------------------------------------------------------------------
// mostly to be used on decltype statements.
template <class T, class... args>
static constexpr auto parampack_get_first (T v, args&&... vargs)
{
  return v;
}
//------------------------------------------------------------------------------
template <template <class...> class Dst, class Src>
struct pass_template_parameters;

template <
  template <class...>
  class Dst,
  template <class...>
  class Src,
  class... Ts>
struct pass_template_parameters<Dst, Src<Ts...>> {
  using type = Dst<Ts...>;
};
//------------------------------------------------------------------------------
template <class T>
static constexpr T* declptr()
{
  return (T*) 0;
}
//------------------------------------------------------------------------------
template <class T>
struct type_wrapper {
  using type = T;
};
//------------------------------------------------------------------------------
// unfortunately doesn't work with non type/class template parameters.
template <template <class...> class T, class U>
struct is_same_template : public std::false_type {};

template <template <class...> class T, class... Ts>
struct is_same_template<T, T<Ts...>> : public std::true_type {};

template <template <class...> class T, class... Ts>
struct is_same_template<T, const T<Ts...>> : public std::true_type {};

template <template <class...> class T, class U>
using is_same_template_t = typename is_same_template<T, U>::type;

template <template <class...> class T, class U>
static constexpr bool is_same_template_v = is_same_template<T, U>::value;

//------------------------------------------------------------------------------
template <class T = void, class... Ts>
static constexpr auto make_array (Ts&&... args)
{
  using U
    = std::conditional_t<std::is_same_v<T, void>, std::common_type_t<Ts...>, T>;
  static_assert (
    !std::disjunction_v<std::is_same<
      std::reference_wrapper<U>,
      std::remove_volatile_t<std::remove_const_t<Ts>>>...>,
    "reference_wrappers are not allowed");
  return std::array<U, sizeof...(Ts)> {std::forward<Ts> (args)...};
}

template <class... Ts>
static constexpr auto make_cstr_array (Ts&&... args)
{
  return make_array<char const*> (std::forward<Ts> (args)...);
}
//------------------------------------------------------------------------------
namespace detail {

template <uint Offset, class T, size_t N_dst, size_t... N>
static constexpr void array_join (std::array<T, N_dst>& dst)
{}

template <uint Offset, class T, size_t N_dst, size_t N_src, size_t... N>
static constexpr void array_join (
  std::array<T, N_dst>&  dst,
  std::array<T, N_src>&& src,
  std::array<T, N>&&... arrs)
{
  for (uint i = 0; i < N_src; ++i) {
    dst[Offset + i] = src[i];
  }
  array_join<Offset + N_src> (dst, std::forward<std::array<T, N>> (arrs)...);
}

}; // namespace detail

template <class T, size_t... N>
static constexpr auto array_join (std::array<T, N>&&... arrs)
{
  constexpr size_t size = (N + ...);

  std::array<T, size> ret;
  detail::array_join<0> (ret, std::forward<std::array<T, N>> (arrs)...);
  return ret;
}

template <uint Start, uint Size, class T, size_t ArrSize>
static constexpr auto array_slice (std::array<T, ArrSize> v)
{
  static_assert (Start + Size <= ArrSize);

  std::array<T, Size> ret;
  for (uint i = 0; i < Size; ++i) {
    ret[i] = v[Start + i];
  }
  return ret;
}
//------------------------------------------------------------------------------
// Double pointers don't const convert, this is annoying on crange.
// https://stackoverflow.com/questions/5055655/double-pointer-const-correctness-warnings-in-c
template <class U, class T, size_t N>
static constexpr auto array_static_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = static_cast<U> (in[i]);
  }
  return ret;
}

template <class U, class T, size_t N>
static constexpr auto array_const_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = const_cast<U> (in[i]);
  }
  return ret;
}

template <class U, class T, size_t N>
static constexpr auto array_reinterpret_cast (std::array<T, N>& in)
{
  std::array<U, N> ret;
  for (uint i = 0; i < N; ++i) {
    ret[i] = reinterpret_cast<U> (in[i]);
  }
  return ret;
}
//------------------------------------------------------------------------------
// Apply a function to each parameter on a variadic argument pack
template <class funct>
static void apply (funct const& f)
{}

template <class funct, typename T, typename... args>
static void apply (funct const& f, T& v, args&&... vargs)
{
  f (v);
  apply (f, std::forward<args> (vargs)...);
}

template <class funct, typename T, typename... args>
static void apply (funct const& f, contiguous_range<T> r, args&&... vargs)
{
  for (int i = 0; i < r.size(); ++i) {
    f (r[i]);
  }
  apply (f, std::forward<args> (vargs)...);
}
//------------------------------------------------------------------------------
// divide memory on blocks. mainly to be used to wrap SIMD loop
template <class T, size_t N, class FunctorB, class FunctorT>
static void block_divide (
  uint              blocksize, // bytes
  std::array<T*, N> elems,
  uint              elem_count,
  FunctorB          f_blocks,
  FunctorT          f_tail)
{
  for (auto v : elems) {
    assert ((((uintptr_t) v) & (uintptr_t {blocksize} - 1)) == 0);
  }

  uint blocks = elem_count / (blocksize / sizeof (T));
  if (blocks) {
    f_blocks (elems, blocks);
  }
  uint offset = blocks * (blocksize / sizeof (T));
  uint tail   = elem_count - offset;
  if (tail) {
    for (auto& v : elems) {
      v += offset;
    }
    f_tail (elems, tail);
  }
}
//------------------------------------------------------------------------------
template <class T, size_t N, size_t Align>
class simd_mem {
public:
  using value_type                  = T;
  static constexpr size_t n_elems   = N;
  static constexpr size_t alignment = Align;

  T& operator[] (uint i)
  {
    assert (i < N);
    return mem[i];
  }
  T const& operator[] (uint i) const
  {
    assert (i < N);
    return mem[i];
  }
  T*               data() { return mem.data(); }
  T const*         data() const { return mem.data(); }
  constexpr size_t size() { return mem.size(); }

private:
  alignas (Align) std::array<T, N> mem;
};
//------------------------------------------------------------------------------
template <class T>
static constexpr T div_ceil (T num, T div)
{
  static_assert (std::is_integral<T>::value, "");
  return (num + div - 1) / div;
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T round_ceil (T num, T round)
{
  static_assert (std::is_integral<T>::value, "");
  return div_ceil (num, round) * round;
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? gcem::pow (T {10.0}, db * T {0.05}) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static T db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? std::exp (db * T {0.05} * T {M_LN10}) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_gain_to_db (T gain, T m_inf_db = T {-130.})
{
  constexpr auto inv_ln10 = 1. / M_LN10;
  return gain > T {0.} ? gcem::log (gain) * *T {20. * inv_ln10} : m_inf_db;
}
//------------------------------------------------------------------------------
template <class T>
static T gain_to_db (T gain, T m_inf_db = T {-130.})
{
  return gain > T {0.} ? std::log10 (gain) * T {20.} : m_inf_db;
}
//------------------------------------------------------------------------------
template <class T>
static T sgn_no_zero (T v, T neg = (T) -1., T pos_zero = (T) 1.)
{
  return (T) ((v < (T) 0) ? neg : pos_zero);
}
//------------------------------------------------------------------------------
template <class T>
static T sgn (T v)
{
  return (T) ((v > (T) 0) - (v < (T) 0));
}
//------------------------------------------------------------------------------
// 0 to 44100/48000. 1 to 88/96KHz, etc... Assumes multiples of 44100 or 48000.
static uint get_samplerate_order (uint sample_rate)
{
  uint sr_order = sample_rate;
  if ((sr_order % 44100) != 0) {
    // assuming multiple of 48Khz
    auto srate_f = (double) sr_order;
    srate_f *= 44100. / 48000.;
    sr_order = (uint) srate_f;
    assert (sr_order % 44100 == 0 && "precission issues");
  }
  sr_order /= 44100; // 1, 2, 4, 8, 16, 32 ...
  sr_order = last_bit_set (sr_order); // 0, 1, 2, 3, 4, 5 ...
  return sr_order;
}
//------------------------------------------------------------------------------
template <uint N>
using k_uint = std::integral_constant<uint, N>;

template <int N>
using k_int = std::integral_constant<int, N>;

template <bool N>
using k_bool = std::integral_constant<bool, N>;

template <s8 N>
using k_s8 = std::integral_constant<s8, N>;

template <u8 N>
using k_u8 = std::integral_constant<u8, N>;

template <s16 N>
using k_s16 = std::integral_constant<s16, N>;

template <u16 N>
using k_u16 = std::integral_constant<u16, N>;

template <s32 N>
using k_s32 = std::integral_constant<s32, N>;

template <u32 N>
using k_u32 = std::integral_constant<u32, N>;

template <s64 N>
using k_s64 = std::integral_constant<s64, N>;

template <u64 N>
using k_u64 = std::integral_constant<u64, N>;
//------------------------------------------------------------------------------
namespace detail {
template <class tuple_like, class Func, size_t... Idxs>
auto tuple_unpack (
  tuple_like&& t,
  Func&&       unpack_f,
  std::index_sequence<Idxs...>)
{
  using ret_t = decltype (unpack_f (std::get<Idxs> (t)...));

  if constexpr (std::is_same_v<ret_t, void>) {
    unpack_f (std::get<Idxs> (t)...);
    return nullptr;
  }
  else {
    return unpack_f (std::get<Idxs> (t)...);
  }
}
} // namespace detail

template <template <class...> class tuple_like, class Func, class... Ts>
auto tuple_unpack (tuple_like<Ts...>&& t, Func&& unpack_f)
{
  return detail::tuple_unpack (
    std::forward<tuple_like<Ts...>> (t),
    std::forward<Func> (unpack_f),
    std::index_sequence_for<Ts...> {});
}
//------------------------------------------------------------------------------
template <class T>
constexpr bool is_aligned_to (uint align, T* v)
{
  assert (is_pow2 (align));
  auto addr = reinterpret_cast<same_size_uint<T*>> (v);
  return (addr & (align - 1)) == 0;
}
//------------------------------------------------------------------------------
struct null_type {};
//------------------------------------------------------------------------------
// Add more stuff if/when required...
template <class T>
struct fraction {
  using value_type = T;

  constexpr void reset (value_type n, value_type d)
  {
    num = n;
    den = d;
  }

  value_type num;
  value_type den;
};
//------------------------------------------------------------------------------
} // namespace artv
