// Non-math (as in math.h) fuctions to help working with vectors

#pragma once

#include <array>
#include <immintrin.h>
#include <type_traits>

#include "artv-common/misc/hana.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/vec.hpp"

namespace artv {

//------------------------------------------------------------------------------
// just an array rounded to a simd width.
template <class T, size_t N, size_t simd_reg_bytes>
using simd_array
  = std::array<T, round_ceil<size_t> (N, (simd_reg_bytes / sizeof (T)))>;

// another array rounded to simd width, but this one contains SIMD vectors.
// "vec_sizeof" is totally equivalent to "simd_reg_bytes" above.
template <class T, size_t N, size_t vec_sizeof>
using simd_vec_array = std::array<
  vec<T, vec_sizeof / sizeof (T)>,
  round_ceil<size_t> (N * sizeof (T), vec_sizeof) / vec_sizeof>;

//------------------------------------------------------------------------------
template <class T, class = void>
struct make_vector;

template <class T>
struct make_vector<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
  using type = vec<T, 1>;
};

template <class T>
struct make_vector<T, enable_if_vec_t<T>> {
  using type = T;
};

template <class T>
using make_vector_t = typename make_vector<T>::type;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_to_intrin (V simdvec)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  if constexpr (std::is_integral_v<T> && traits.bytes == sse_bytes) {
    return *reinterpret_cast<__m128i*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, double> && traits.bytes == sse_bytes) {
    return *reinterpret_cast<__m128d*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, float> && traits.bytes == sse_bytes) {
    return *reinterpret_cast<__m128*> (&simdvec);
  }
#if 0
  //TODO: these types are not found when compiling for Windows. Keeping them
  // disabled until they are needed, as I wasn't able to fix it quickly.
  else if constexpr (std::is_integral_v<T> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256i*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, double> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256d*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, float> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256*> (&simdvec);
  }
#endif
  else {
    static_assert (!std::is_same_v<V, V>, "Unkown intrinsic conversion");
    return 0;
  }
}
//------------------------------------------------------------------------------
// slice N elements from vec as a new vector of the appropiate size
template <uint N, class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_slice (V src, uint offset)
{
  using T               = vec_value_type_t<V>;
  constexpr auto traits = vec_traits<V>();

  static_assert (traits.size >= N);
  assert ((offset + N) <= traits.size);

  vec<T, N> dst;
  // no memcpy, ordering inside the vector not guaranteed.
  for (uint i = 0; i < N; ++i) {
    dst[i] = src[offset + i];
  }
  return dst;
}
//------------------------------------------------------------------------------
// copy all elements from "src" into "dst" at a given offset.
template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_cp_slice (V1& dst, V2 src, uint offset)
{
  using T = vec_value_type_t<V1>;
  static_assert (std::is_same_v<T, vec_value_type_t<V2>>);

  constexpr auto traits1 = vec_traits<V1>();
  constexpr auto traits2 = vec_traits<V2>();

  static_assert (traits1.size >= traits2.size);
  assert ((offset + traits2.size) <= traits1.size);

  // no memcpy, ordering inside the vector not guaranteed.
  for (uint i = 0; i < traits2.size; ++i) {
    dst[offset + i] = src[i];
  }
  return dst;
}
//------------------------------------------------------------------------------
// Notice:
//
// - All these loads, stores and casts could be done with "memcpy" but I
//   lazily went with __may_alias__ and type punning to avoid verifying that
//   memcpy does the right thing.
// - No object orientation, as it would defeat the purpuse of using the
//   builtin compiler wrappers.
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V vec_load (vec_value_type_t<V> const* src)
{
  vec_traits<V>(); // check type validity only
  return *reinterpret_cast<V const*> (src);
}

template <uint N, class T>
constexpr inline vec<T, N> vec_load (T const* src)
{
  return vec_load<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_load (V& dst, vec_value_type_t<V> const* src)
{
  dst = vec_load<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V vec_load (xspan<vec_value_type_t<V> const> src)
{
  constexpr auto traits = vec_traits<V>();
  assert (src.size() >= traits.size);
  return vec_load<V> (src.data());
}

template <uint N, class T>
constexpr inline vec<T, N> vec_load (xspan<T const> src)
{
  return vec_load<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_load (V& dst, xspan<vec_value_type_t<V> const> src)
{
  dst = vec_load<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V vec_load_unaligned (vec_value_type_t<V> const* src)
{
  constexpr auto traits = vec_traits<V>();
  using vec_u_type      = typename decltype (traits)::type_u;

  return *reinterpret_cast<vec_u_type const*> (src);
}

template <uint N, class T>
constexpr inline vec<T, N> vec_load_unaligned (T const* src)
{
  return vec_load_unaligned<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_load_unaligned (
  V&                         dst,
  vec_value_type_t<V> const* src)
{
  dst = vec_load_unaligned<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V vec_load_unaligned (xspan<vec_value_type_t<V> const> src)
{
  constexpr auto traits = vec_traits<V>();
  assert (src.size() >= traits.size);
  return vec_load_unaligned<V> (src.data());
}

template <uint N, class T>
constexpr inline vec<T, N> vec_load_unaligned (xspan<T const> src)
{
  return vec_load_unaligned<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_load_unaligned (
  V&                               dst,
  xspan<vec_value_type_t<V> const> src)
{
  dst = vec_load_unaligned<V> (src);
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_store (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<V*> (dst) = src;
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_store (xspan<vec_value_type_t<V>> dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store (dst.data(), src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_store_unaligned (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using vec_u_type      = typename decltype (traits)::type_u;
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<vec_u_type*> (dst) = src;
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline void vec_store_unaligned (
  xspan<vec_value_type_t<V>> dst,
  V                          src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store_unaligned (dst.data(), src);
}
//------------------------------------------------------------------------------
// vec_set: broadcast a value to a vector
template <class V, enable_if_vec_or_scalar_t<V>* = nullptr>
constexpr inline void vec_set (V& dst, vec_value_type_t<V> v)
{
  constexpr auto traits = vec_traits<V>();
  dst                   = v - V {}; // convoluted but works.
}

template <class V, enable_if_vec_or_scalar_t<V>* = nullptr>
constexpr inline V vec_set (vec_value_type_t<V> v)
{
  V ret;
  vec_set (ret, v);
  return ret;
}

template <uint N, class T>
constexpr inline vec<T, N> vec_set (T v)
{
  vec<T, N> ret;
  vec_set (ret, v);
  return ret;
}
//------------------------------------------------------------------------------
// Frequently when having builtin types it's desired to conditionally transform
// to vec<x,1> for interfacting with functions that only take vectors. These
// overloads are for helping. The idea is that:
//
// - make_vec(x): will let a vector type pass or make a vector of 1 element of
//                from a builtin type.
//
// - make_vec<V>(x):
//    - when "V" is a vector and "x" a builtin: returns a vector of N elements
//      with the value in "x" broadcasted.
//    - when "V" is a builtin and "x" is a builtin returns a vector of size 1.
//    - when "V" is a vector and "x" is a vector it is passed through.

template <
  class T,
  std::enable_if_t<std::is_arithmetic_v<T> && !is_vec_v<T>>* = nullptr>
constexpr inline auto make_vec (T v)
{
  return vec_set<1> (v);
}

// casting overload, requires providing the type.
template <
  class Dst_vec,
  class T,
  std::enable_if_t<std::is_arithmetic_v<T> && is_vec_v<Dst_vec>>* = nullptr>
constexpr inline auto make_vec (T v)
{
  return vec_set<Dst_vec> (static_cast<vec_value_type_t<Dst_vec>> (v));
}

// passthrough
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V make_vec (V v)
{
  return v;
}

//------------------------------------------------------------------------------
// convenience functions. Created basically to deal with vectors of size 1
// conversion, so they should result in no-ops after optimization.
template <uint VecN, class T, size_t Size>
constexpr inline auto vec_array_wrap (std::array<T, Size> v)
{
  static_assert (Size % VecN == 0);

  std::array<vec<T, VecN>, Size / VecN> ret;
  for (uint i = 0; i < ret.size(); ++i) {
    for (uint j = 0; j < VecN; ++j) {
      ret[i][j] = v[i * VecN + j];
    }
  }
  return ret;
}

template <class VecType, size_t Size>
constexpr inline auto vec_array_unwrap (std::array<VecType, Size> v)
{
  static_assert (is_vec_v<VecType>);
  using traits     = vec_traits_t<VecType>;
  using value_type = typename traits::value_type;

  std::array<value_type, traits::size * Size> ret;

  for (uint i = 0; i < ret.size(); ++i) {
    for (uint j = 0; j < traits::size; ++j) {
      ret[i * traits::size + j] = v[i][j];
    }
  }
  return ret;
}

template <class T, size_t Size>
constexpr inline auto vec1_array_wrap (std::array<T, Size> v)
{
  return vec_array_wrap<1> (v);
}

template <class VecType, size_t Size>
constexpr inline auto vec1_array_unwrap (std::array<VecType, Size> v)
{
  return vec_array_unwrap (v);
}

//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_to_array (V v)
{
  using T = vec_value_type_t<V>;
  std::array<T, vec_traits_t<V>::size> arr;
  // no memcpy, memory ordering is unspecified
  for (uint i = 0; i < arr.size(); ++i) {
    arr[i] = v[i];
  }
  return arr;
}

template <class T, size_t Size>
constexpr inline auto vec_from_array (std::array<T, Size> arr)
{
  vec<T, Size> v;
  // no memcpy, memory ordering is unspecified
  for (uint i = 0; i < Size; ++i) {
    v[i] = arr[i];
  }
  return v;
}

//------------------------------------------------------------------------------
#if 1
// This is the best implementation, as "__builtin_shufflevector" requires
// constants as the indexes. This means that a function would require passing
// std::integral_constant and a lot uglyness to just achieve the same that this
// macro. The ifdefed function (under #else) below doesn't work. Kept for
// ilustration purposes
#define vec_shuffle __builtin_shufflevector
#else
template <class V, class... Ts>
inline V vec_shuffle (V a, V b, Ts... indexes)
{
  constexpr auto traits = vec_traits<V>();
  using common_t        = std::common_type_t<Ts...>;
  static_assert (std::is_integral_v<common_t>, "invalid index type");
  static_assert (sizeof...(Ts) == traits.size, "incorrect index count");

#if !defined(NDEBUG)
  auto indexes_tuple = hana::make<hana::tuple_tag> (indexes...);
  hana::for_each (indexes_tuple, [&] (auto&& idx) {
    assert (idx >= 0 && idx < (traits.size * 2));
  });
#endif
  return __builtin_shufflevector (a, b, indexes...);
}
#endif
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
using vec_array_type = std::array<vec_value_type_t<V>, vec_traits_t<V>::size>;
//------------------------------------------------------------------------------
// "vec_cat" from an array
template <class V, size_t N, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_cat (std::array<V, N> a)
{
  using T                 = vec_value_type_t<V>;
  using traits            = vec_traits_t<V>;
  constexpr uint vec_size = traits::size;
  constexpr uint size     = traits::size * N;
  static_assert ((size % 2) == 0);

  vec<T, size> ret;
  for (uint i = 0; i < size; ++i) {
    auto const vec_idx  = i / vec_size;
    auto const elem_idx = i % vec_size;
    ret[i]              = a[vec_idx][elem_idx];
  }
  return ret;
}

// "vec_cat" from an array
template <class V, class... Ts, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_cat (V v1, Ts&&... vecs)
{
  constexpr size_t n_elems = 1 + sizeof...(Ts);
  static_assert ((n_elems % 2) == 0);
  return vec_cat (std::array<V, n_elems> {{v1, std::forward<Ts> (vecs)...}});
}
//------------------------------------------------------------------------------
// split vector in an array of vectors of a smaller (divisible) vector type.
template <uint N, class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_split (V src)
{
  using T               = vec_value_type_t<V>;
  constexpr auto traits = vec_traits<V>();

  static_assert (traits.size >= N);
  static_assert ((traits.size % N) == 0);

  constexpr uint n_elems = traits.size / N;

  std::array<vec<T, N>, n_elems> dst;
  // no memcpy, ordering inside the vector not guaranteed.
  for (uint i = 0; i < traits.size; ++i) {
    dst[i / N][i % N] = src[i];
  }
  return dst;
}
//------------------------------------------------------------------------------
// init a vector from scalars, type deduced from the first element
template <class T, class... Ts>
constexpr inline auto vec_init (T&& a, Ts&&... b)
{
  constexpr uint n_elems = sizeof...(Ts) + 1;
  static_assert (!is_vec_v<T>);
  static_assert (std::is_arithmetic_v<T>);
  static_assert ((n_elems % 2) == 0);
  using V = vec<T, n_elems>;

  alignas (V) vec_array_type<V> dst;
  dst[0] = a;
  uint i = 1;
  hana::for_each (hana::make<hana::tuple_tag> (b...), [&] (auto&& idx) {
    using U = std::remove_reference_t<std::remove_cv_t<decltype (idx)>>;
    static_assert (!is_vec_v<U>);
    static_assert (std::is_arithmetic_v<U>);
    dst[i] = static_cast<T> (idx);
    ++i;
  });
  return vec_load<V> (dst);
};
//------------------------------------------------------------------------------
// Vector cast.
//
// - If "V" is a vector and "T" is scalar it will cast the input parameter to a
//   vector of the same size.
//
// - If "V" is scalar and "T" a vector, the function is equivalent to a
//   static_cast to the value_type of "V" followed by a vector broadcast
//  (setting all vector elements static_casted value).
//
// - If "V" and "T" are scalars, the function is equivalent to static_cast.
//
// - If "V" and "T" are vectors, both vectors have to be of the same size. In
//   that case the result will be a vector cased to the value_type of "V".
//
// Notice that when one vector is scaled another requiring bigger SIMD with
// "-Wpsabi" warnings might be generated. The ABI-related warnings are no
// problem if there are no vectors in shared library interfaces. Suppressing has
// to be done globally or in place unfortunately.
//
template <
  class T,
  class V,
  std::enable_if_t<is_vec_or_scalar_v<T> && is_vec_or_scalar_v<V>>* = nullptr>
constexpr inline auto vec_cast (V a)
{
  using src_traits = vec_traits_t<V>;
  using dst_traits = std::conditional_t<
    std::is_arithmetic_v<T>,
    typename src_traits::template rebind_traits_same_size<T>,
    vec_traits_t<T>>;

  if constexpr (src_traits::size != 0) {
    static_assert (src_traits::size == dst_traits::size, "sizes must match");
    return __builtin_convertvector (a, typename dst_traits::type);
  }
  else if constexpr (dst_traits::size != 0) {
    // scalar to vector
    return vec_set<T> (static_cast<vec_value_type_t<T>> (a));
  }
  else {
    // scalar to scalar (passthrough)
    return static_cast<T> (a);
  }
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
constexpr inline auto vec_abs (V&& x)
{
  return (x < (vec_value_type_t<V>) 0) ? -x : x;
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_min (V1&& x, V2&& y)
{
  return x < y ? x : y;
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_min (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_min (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_min (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_min (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_max (V1&& x, V2&& y)
{
  return x > y ? x : y;
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_max (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_max (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_sgn (V x)
{
  using T = vec_value_type_t<V>;
  return (x < (T) 0) ? vec_set<V> (-1) : vec_set<V> (1);
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  class V3,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2> && is_vec_v<V3>>* = nullptr>
constexpr inline auto vec_clamp (V1&& x, V2&& min, V3&& max)
{
  return vec_min (vec_max (x, min), max);
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_clamp (
  V&&                 x,
  vec_value_type_t<V> y,
  vec_value_type_t<V> z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_clamp (std::forward<V> (x), vec_set<Vv> (y), vec_set<Vv> (z));
}

template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_clamp (V1&& x, V2&& y, vec_value_type_t<V1> z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    std::forward<V1> (x), std::forward<V2> (y), vec_set<Vv> (z));
}

template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_clamp (V1&& x, vec_value_type_t<V1> y, V2&& z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    std::forward<V1> (x), vec_set<Vv> (y), std::forward<V2> (z));
}

template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
constexpr inline auto vec_clamp (vec_value_type_t<V1> x, V1&& y, V2&& z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    vec_set<Vv> (x), std::forward<V1> (y), std::forward<V2> (z));
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_clamp (
  vec_value_type_t<V> x,
  vec_value_type_t<V> y,
  V&&                 z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_clamp (vec_set<V> (x), vec_set<Vv> (y), std::forward<V> (z));
}

template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_clamp (
  vec_value_type_t<V> x,
  V&&                 y,
  vec_value_type_t<V> z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_clamp (vec_set<V> (x), std::forward<V> (y), vec_set<Vv> (z));
}

//------------------------------------------------------------------------------
template <
  class V1,
  class V2 = std::decay_t<V1>,
  class V3 = std::decay_t<V1>,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2> && is_vec_v<V3>>* = nullptr>
constexpr inline auto vec_sgn_no_zero (
  V1&& x,
  V2&& neg      = vec_set<std::decay_t<V1>> ((vec_value_type_t<V1>) -1.),
  V3&& pos_zero = vec_set<std::decay_t<V1>> ((vec_value_type_t<V1>) 1.))
{
  using Vv              = std::common_type_t<V1, V2, V3>;
  constexpr auto traits = vec_traits<Vv>();
  using value_type      = vec_value_type_t<Vv>;
  return (x < ((value_type) 0.)) ? neg : pos_zero;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_inner_and (V v)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;
  static_assert (std::is_integral_v<T>);

  T r {v[0]};
  for (uint i = 1; i < traits.size; ++i) {
    r &= v[i];
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_inner_or (V v)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;
  static_assert (std::is_integral_v<T>);

  T r {v[0]};
  for (uint i = 1; i < traits.size; ++i) {
    r |= v[i];
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_inner_add (V v)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  T r {v[0]};
  for (uint i = 1; i < traits.size; ++i) {
    r += v[i];
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline auto vec_inner_mul (V v)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  T r {v[0]};
  for (uint i = 1; i < traits.size; ++i) {
    r |= v[i];
  }
  return r;
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline bool vec_is_all_zeros (V v)
{
  return !vec_inner_or (v);
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline bool vec_is_all_ones (V v)
{
  return !(~vec_inner_and (v));
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
constexpr inline V zero_to_lowest (V v)
{
  using T         = vec_value_type_t<V>;
  constexpr T min = -std::numeric_limits<T>::lowest();
  return (v != (T) 0) ? v : vec_set<V> (min);
}
//------------------------------------------------------------------------------
} // namespace artv
