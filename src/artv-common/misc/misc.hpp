#pragma once

// miscellanious things that probably don't deserve its own header.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <tuple>
#include <type_traits>
#include <utility>

#include <gcem.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

// TODO: Split this.

namespace artv {

#ifdef __GNUC__
#define likely(x)   __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)
#else
#define likely(x)   x
#define unlikely(x) x
#endif

#define VERSION_GET(major, minor, rev) (major * 1000000 + minor * 1000 + rev)

#define PP_EXPAND(x) x
#define PP_TO_STR(x) #x

// usage print_type_as_error<whatever> v; The compiler will print the type.
template <class... Ts>
struct print_type_as_error;

template <class T>
static constexpr void print_type_as_error_f (T&& v)
{
  print_type_as_error<T> dummy;
}
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
struct is_same_template<T, T<Ts...> const> : public std::true_type {};

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

template <size_t N, class T>
static constexpr auto array_broadcast (T v)
{
  std::array<T, N> ret;
  ret.fill (v);
  return ret;
}

template <class U, class T, size_t N>
static constexpr auto array_cast (std::array<T, N> v)
{
  std::array<U, N> ret {};
  for (uint i = 0; i < N; ++i) {
    ret[i] = static_cast<U> (v[i]);
  }
  return ret;
}

template <size_t N1, size_t N2, class T>
static constexpr auto array_flatten (std::array<std::array<T, N2>, N1> v)
{
  std::array<T, N1 * N2> ret;
  for (uint i = 0; i < v.size(); ++i) {
    for (uint j = 0; j < v[i].size(); ++j) {
      ret[i * v.size() + j] = v[i][j];
    }
  }
  return ret;
}

//------------------------------------------------------------------------------
namespace detail {

template <uint Offset, class T, size_t N_dst, size_t... N>
static constexpr void array_cat (std::array<T, N_dst>& dst)
{}

template <uint Offset, class T, size_t N_dst, size_t N_src, size_t... N>
static constexpr void array_cat (
  std::array<T, N_dst>&       dst,
  std::array<T, N_src> const& src,
  std::array<T, N> const&... arrs)
{
  memcpy (&dst[Offset], src.data(), sizeof src);
  array_cat<Offset + N_src> (dst, arrs...);
}

}; // namespace detail

template <class T, size_t... N>
static constexpr auto array_cat (std::array<T, N> const&... arrs)
{
  constexpr size_t size = (N + ...);

  std::array<T, size> ret;
  detail::array_cat<0> (ret, arrs...);
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
// Double pointers don't const convert, this is annoying on xspan.
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
// Apply a function to each parameter on a variadic argument pack. iterating
// stuff (could be much smarter iterating containers, but it's not needed now)
template <class funct>
static void apply (funct const& f)
{}

template <class funct, typename T, typename... args>
static void apply (funct const& f, T&& v, args&&... vargs)
{
  f (v);
  apply (f, std::forward<args> (vargs)...);
}

template <class funct, typename T, typename... args>
static void apply (funct const& f, xspan<T> r, args&&... vargs)
{
  for (int i = 0; i < r.size(); ++i) {
    f (r[i]);
  }
  apply (f, std::forward<args> (vargs)...);
}

template <class funct, typename T, std::size_t N, typename... args>
static void apply (funct const& f, std::array<T, N>&& r, args&&... vargs)
{
  apply (f, xspan {r}, std::forward<args> (vargs)...);
}

template <class funct, typename T, std::size_t N, typename... args>
static void apply (funct const& f, std::array<T, N>& r, args&&... vargs)
{
  apply (f, xspan {r}, std::forward<args> (vargs)...);
}
//------------------------------------------------------------------------------
// divide memory on blocks. Invoke a function for each block and one for the
// remainder, if any.
template <class T, size_t N, class FunctorB, class FunctorT>
static void block_divide (
  uint              blocksize, // bytes
  std::array<T*, N> elems,
  uint              elem_count,
  FunctorB          f_blocks,
  FunctorT          f_tail)
{
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
static constexpr T round_floor (T num, T round)
{
  static_assert (std::is_integral<T>::value, "");
  return (num / round) * round;
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? gcem::exp (db * (T) (M_LN10 / 20.)) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static T db_to_gain (T db, T m_inf_db = T {-130.})
{
  return db > m_inf_db ? exp (db * (T) (M_LN10 / 20.)) : T {0.};
}
//------------------------------------------------------------------------------
template <class T>
static constexpr T constexpr_gain_to_db (T gain, T m_inf_db = T {-130.})
{
  auto absval = (T) gcem::abs (gain);
  return gain > T {0.} ? gcem::log (absval) * (T) (20. / M_LN10) : m_inf_db;
}
//------------------------------------------------------------------------------
template <class T>
static T gain_to_db (T gain, T m_inf_db = T {-130.})
{
  auto absval = (T) abs (gain);
  return gain > T {0.} ? log (absval) * (T) (20. / M_LN10) : m_inf_db;
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
// 0 to 44100/48000. 1 to 88/96KHz, etc... Assumes multiples of 44100 or
// 48000.
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
template <class T, size_t D1, size_t D2>
using array2d = std::array<std::array<T, D1>, D2>;
//------------------------------------------------------------------------------
constexpr void constexpr_memcpy (void* dst, void const* src, std::size_t size)
{
  __builtin_memcpy (dst, src, size);
}
//------------------------------------------------------------------------------
template <int I>
struct priority_tag : priority_tag<I - 1> {};
template <>
struct priority_tag<0> {};
//------------------------------------------------------------------------------
namespace detail {
// Inspired on this
// https://twitter.com/ericniebler/status/852192542653329408

template <class T>
char (&is_array_subscriptable (priority_tag<0>))[1]; // no

template <class T, class = decltype (std::declval<T>()[0])>
char (&is_array_subscriptable (priority_tag<1>))[2]; // yes
} // namespace detail

template <typename T>
struct is_array_subscriptable
  : k_bool<
      sizeof (detail::is_array_subscriptable<T> (priority_tag<1> {})) == 2> {};

template <typename T>
constexpr bool is_array_subscriptable_v = is_array_subscriptable<T>::value;
//------------------------------------------------------------------------------
// forward a value to a lambda, useful for e.g. slightly modifying constexpr
// variables.
template <class V, class F>
constexpr decltype (auto) lambda_forward (V&& value, F&& func)
{
  return func (std::forward<V> (value));
}
//------------------------------------------------------------------------------
template <uint Val, class T>
struct index_seq_add;

template <uint Val, uint... Idx>
struct index_seq_add<Val, std::index_sequence<Idx...>> {
  using type = std::index_sequence<(Val + Idx)...>;
};

template <uint Offset, class T>
using index_seq_add_t = typename index_seq_add<Offset, T>::type;
//------------------------------------------------------------------------------
template <uint Val, class T>
struct index_seq_mul;

template <uint Val, uint... Idx>
struct index_seq_mul<Val, std::index_sequence<Idx...>> {
  using type = std::index_sequence<(Idx * Val)...>;
};

template <uint Val, class T>
using index_seq_mul_t = typename index_seq_mul<Val, T>::type;
//------------------------------------------------------------------------------
// visit a variadic pack in a lambda
//------------------------------------------------------------------------------
template <std::size_t... Idxs, class... Ts, class Fn>
constexpr auto vpack_visit (
  std::index_sequence<Idxs...>,
  Fn&& visitor,
  Ts&&... args)
{
  auto tpl = std::forward_as_tuple (args...);
  return visitor (std::get<Idxs> (tpl)...);
}

template <uint Offset, uint N, class... Ts, class Fn>
constexpr auto vpack_visit (Fn&& visit, Ts&&... args)
{
  static_assert ((N + Offset) <= sizeof...(args));
  return vpack_visit (
    index_seq_add_t<Offset, std::make_index_sequence<N>> {},
    std::forward<Fn> (visit),
    std::forward<Ts> (args)...);
}

template <uint Offset, class... Ts, class Fn>
constexpr auto vpack_visit (Fn&& visit, Ts&&... args)
{
  static_assert (Offset <= sizeof...(args));
  return vpack_visit<Offset, sizeof...(Ts) - Offset> (
    std::forward<Fn> (visit), std::forward<Ts> (args)...);
}

template <class... Ts, class Fn>
constexpr auto vpack_visit (Fn&& visit, Ts&&... args)
{
  return vpack_visit<0> (std::forward<Fn> (visit), std::forward<Ts> (args)...);
}
//------------------------------------------------------------------------------
template <std::size_t... Idxs, class... Ts>
constexpr auto forward_range_as_tuple (
  std::index_sequence<Idxs...> s,
  Ts&&... args)
{
  return vpack_visit (
    s,
    [] (auto&&... targs) { return std::forward_as_tuple (targs...); },
    std::forward<Ts> (args)...);
}

template <uint Offset, uint N, class... Ts>
constexpr auto forward_range_as_tuple (Ts&&... args)
{
  return vpack_visit<Offset, N> (
    [] (auto&&... targs) { return std::forward_as_tuple (targs...); },
    std::forward<Ts> (args)...);
}

template <uint Offset, class... Ts>
constexpr auto forward_range_as_tuple (Ts&&... args)
{
  return vpack_visit<Offset> (
    [] (auto&&... targs) { return std::forward_as_tuple (targs...); },
    std::forward<Ts> (args)...);
}

//------------------------------------------------------------------------------
#if 0
//clang-format off
template <class P, class M>
size_t cpp_offsetof (M P::const* member)
{
  return (size_t) & (reinterpret_cast<P*> (0)->*member);
}

template <class P, class M>
P* container_of (M* ptr, M P::const* member)
{
  return (P*) ((char*) ptr - cpp_offsetof (member));
}
//clang-format on
#endif
} // namespace artv
