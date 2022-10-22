#pragma once

// A header to make a bit less painful to work with native GCC/clang vector
// types

// This was previously using xsimd direcly, but inteface-breaking changes
// motivated this intermediate wrapping (unfortunately, boring). As this code is
// only for Desktop PC and probably desktop ARM, Clang (and maybe GCC) will be
// used. Leveraging the builtin vector extensions instead of OO wrappers.
//
// Unfortunately at the time of writing GCC support for vector types is buggy
// too on C++. This code requires clang 11 or higher:
// https://gcc.gnu.org/bugzilla//show_bug.cgi?id=57572

#include <type_traits>

#include "artv-common/misc/short_ints.hpp"

namespace artv {

#if defined(__AVX512F__) || defined(__AVX512CD__) || defined(__AVX512DQ__) \
  || defined(__AVX512BW__)
static constexpr uint vec_max_bytes = 64;

#elif defined(__AVX__) || defined(__AVX2__)
static constexpr uint vec_max_bytes = 32;

#elif defined(__SSE2__) || defined(__SSE3__) || defined(__SSSE3__) \
  || defined(__SSE4_1__) || defined(__SSE4_2__) || defined(__ARM_NEON)
static constexpr uint vec_max_bytes = 16;

#else
static constexpr uint vec_max_bytes = 0;
#endif
//------------------------------------------------------------------------------
template <class T, uint N>
struct simd_vector_traits {
  static_assert (std::is_arithmetic_v<T>, "");
  static_assert (N > 0, "");

  using value_type            = T;
  static constexpr uint size  = N;
  static constexpr uint bytes = N * sizeof (T);

  // allows aligned loads by casting
  using type __attribute__ ((vector_size (bytes), __may_alias__)) = T;
  // To allow unaligned loads by casting
  using type_u
    __attribute__ ((vector_size (bytes), __may_alias__, __aligned__ (1)))
    = T;

  static constexpr uint alignment = alignof (type);

  template <class U>
  using rebind __attribute__ ((vector_size (bytes), __may_alias__)) = U;

  using same_size_int_type  = rebind<same_size_int<T>>;
  using same_size_uint_type = rebind<same_size_uint<T>>;

  template <class U>
  using rebind_traits_same_bytes = simd_vector_traits<U, bytes / sizeof (U)>;

  template <class U>
  using rebind_traits_same_size = simd_vector_traits<U, N>;

  // As attributes don't participate on template argument deduction or
  // specialization, only the type of the vector can be deduced. Fortunately
  // deducing the type is mostly enough when wrapping real SIMD instruction is
  // the goal, as the "vector_size" matches "sizeof(vec)" for vectors
  // containing a number of elements that is a power of two.
  //
  // As "sizeof" the vectors grow in powers of 2
  // (e.g. "sizeof (vec<int, 5>) = 32" (on x64)), we restrict our vectors to
  // power of two sizes. At the time of writing this we can't detect the vector
  // size attribute. At the time of writing this e.g. GCC doesn't allow sizes
  // that are not powers of 2 anyways.
  static_assert (bytes > 0);
  static_assert ((bytes & (bytes - 1)) == 0, "bytes is not a power of 2.");
};

namespace detail {
//------------------------------------------------------------------------------
template <class T, uint N>
using vec = typename simd_vector_traits<T, N>::type;
template <class T, uint N>
using vec_u = typename simd_vector_traits<T, N>::type_u;
//------------------------------------------------------------------------------
template <uint size_of, class T>
static inline constexpr auto get_traits (
  T __attribute__ ((vector_size (size_of))))
{
  return simd_vector_traits<T, size_of / sizeof (T)> {};
}

template <uint size_of, class T>
static inline constexpr auto get_traits (T)
{
  return nullptr;
}

template <class V>
static constexpr auto get_traits()
{
  return get_traits<sizeof (V)> (V {});
}
//------------------------------------------------------------------------------
// Template specialization doesn't play nice with attributes, using constexpr...
template <uint size_of, class T>
static constexpr bool is_vec (T)
{
  return false;
}

template <uint size_of, class T>
static constexpr bool is_vec (T __attribute__ ((vector_size (size_of))))
{
  return true;
}

// caveats: to implement "is_vec", aka attribute detection, only constexpr
// overloaded functions work on Clang13+ (totally broken on GCC). This puts a
// requirement on that every class passed to "is_vec" must be default constexpr
// constructible, otherwise the compiler will complain.
//
// This renders some of the SFINAE overload usages a bit broken when passing a
// random class. This is left as-is, as it works mostly as intended by catching
// user errors, it's that the error message will be about constant expressions
// instead of about non found vector functions, hence this comment to clarify.
template <class T>
static constexpr bool is_vec_v
  = is_vec<sizeof (T)> (std::remove_reference_t<std::remove_cv_t<T>> {});
//------------------------------------------------------------------------------
// for unevaluated contexts
template <class V>
using vec_traits_t
  = decltype (get_traits<std::remove_reference_t<std::remove_cv_t<V>>>());

template <class, bool>
struct vec_vt_enabler;

template <class V>
struct vec_vt_enabler<V, false> {
  using type = void;
};

template <class V>
struct vec_vt_enabler<V, true> {
  using type = typename vec_traits_t<V>::value_type;
};

template <class V>
using vec_value_type_t = typename vec_vt_enabler<V, is_vec_v<V>>::type;

} // namespace detail

//------------------------------------------------------------------------------
// TYPES
//------------------------------------------------------------------------------
template <class T, uint N>
using vec = detail::vec<T, N>;

// Gets the vector traits (see "simd_vector_traits") of a vector type "V".
// removes references and CV qualification on "V"
template <class V>
using vec_traits_t = detail::vec_traits_t<V>;

// Gets the value type (e.g float) of the vector type "V". removes references
// and CV qualification on "V".
template <class V>
using vec_value_type_t = detail::vec_value_type_t<V>;

// gets if "V" is a vector. removes references and CV qualification on "V".
template <class V>
static constexpr bool is_vec_v = detail::is_vec_v<V>;

// gets if "V" is a floating point vector. removes references and CV
// qualification on "V".
// clang-format off
template <class V>
static constexpr bool is_vec_of_float_type_v
  = is_vec_v<V> && std::is_floating_point_v<vec_value_type_t<V>>;
// clang-format on

using float_x1 = vec<float, 1>;
using f32_x1   = float_x1;

// Convenience half simd, normally used to avoid duplicating two channels
using float_x2  = vec<float, 2>;
using double_x1 = vec<double, 1>;
using f32_x2    = float_x2;
using f64_x1    = double_x1;
using s8_x8     = vec<s8, 8>;
using u8_x8     = vec<u8, 8>;
using s16_x4    = vec<s16, 4>;
using u16_x4    = vec<u16, 4>;
using s32_x2    = vec<s32, 2>;
using u32_x2    = vec<u32, 2>;

template <
  class T,
  std::enable_if_t<sizeof (T) <= 8 && std::is_arithmetic_v<T>>* = nullptr>
using vec8 = vec<T, 8 / sizeof (T)>;

// SSE or equivalent
using float_x4  = vec<float, 4>;
using double_x2 = vec<double, 2>;
using f32_x4    = float_x4;
using f64_x2    = double_x2;
using s8_x16    = vec<s8, 16>;
using u8_x16    = vec<u8, 16>;
using s16_x8    = vec<s16, 8>;
using u16_x8    = vec<u16, 8>;
using s32_x4    = vec<s32, 4>;
using u32_x4    = vec<u32, 4>;
using s64_x2    = vec<s64, 2>;
using u64_x2    = vec<u64, 2>;

template <
  class T,
  std::enable_if_t<sizeof (T) <= 16 && std::is_arithmetic_v<T>>* = nullptr>
using vec16 = vec<T, 16 / sizeof (T)>;

// AVX or equivalent.
using float_x8  = vec<float, 8>;
using double_x4 = vec<double, 4>;
using f32_x8    = float_x8;
using f64_x4    = double_x4;
using s8_x32    = vec<s8, 32>;
using u8_x32    = vec<u8, 32>;
using s16_x16   = vec<s16, 16>;
using u16_x16   = vec<u16, 16>;
using s32_x8    = vec<s32, 8>;
using u32_x8    = vec<u32, 8>;
using s64_x4    = vec<s64, 4>;
using u64_x4    = vec<u64, 4>;

template <
  class T,
  std::enable_if_t<sizeof (T) <= 32 && std::is_arithmetic_v<T>>* = nullptr>
using vec32 = vec<T, 32 / sizeof (T)>;

static constexpr uint sse_bytes = 16;
static constexpr uint avx_bytes = 32;

template <class V>
using enable_if_vec_t = std::enable_if_t<is_vec_v<V>>;

template <class V>
using enable_if_vec_of_float_point_t
  = std::enable_if_t<is_vec_of_float_type_v<V>>;

template <class V, class T>
using enable_if_vec_value_type_t
  = std::enable_if_t<is_vec_v<V> && std::is_same_v<vec_value_type_t<V>, T>>;

//------------------------------------------------------------------------------
// TODO: As of now this is only done for x64, but it can be easily ifdefed
//------------------------------------------------------------------------------
// Notice:
//
// - All these loads, stores and casts could be done with "memcpy" but I
//   lazily went with __may_alias__ and type punning to avoid verifying that
//   memcpy does the right thing.
// - No object orientation, as it would defeat the purpuse of using the
//   builtin compiler wrappers.
//------------------------------------------------------------------------------
// FUNCTIONS
//------------------------------------------------------------------------------
// returns a "simd_vector_traits"
template <class V, enable_if_vec_t<V>* = nullptr>
static constexpr auto vec_traits()
{
  return vec_traits_t<V> {};
}

template <class V, enable_if_vec_t<V>* = nullptr>
static constexpr auto vec_traits (V)
{
  return vec_traits<V>();
}
//------------------------------------------------------------------------------
} // namespace artv
