#pragma once

// This was previously using xsimd direcly, but inteface-breaking changes
// motivated this intermediate wrapping (unfortunately, boring). As this code is
// only for Desktop PC and probably desktop ARM, Clang (and maybe GCC) will be
// used. Leveraging the builtin vector extensions instead of OO wrappers.
//
// Unfortunately at the time of writing GCC support for vector types is buggy
// too on C++. This code requires clang 11 or higher:
// https://gcc.gnu.org/bugzilla//show_bug.cgi?id=57572

#define XSIMD_DISABLED 1
#define XSIMD_NEW_INTERFACE 1

#include <array>
#include <cmath>
#include <immintrin.h>
#include <tuple>
#include <type_traits>
#include <utility>
#include <xsimd/xsimd.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"
#if 0 // was required when maybe implementing __shufflevector wrapper.
#include "artv-common/misc/hana.hpp"
#endif

namespace artv {
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
  using rebind_traits = simd_vector_traits<U, bytes / sizeof (U)>;

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
  static_assert (is_pow2 (bytes), "See comment above this assertion.");
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
  T __attribute__ ((vector_size (size_of), __may_alias__)))
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
static constexpr bool is_vec (
  T __attribute__ ((vector_size (size_of), __may_alias__)))
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

//------------------------------------------------------------------------------
// TODO: As of now this is only done for x64, but it can be easily ifdefed
//------------------------------------------------------------------------------
template <uint size>
struct to_xsimd_arch;

#if XSIMD_NEW_INTERFACE

// TODO sse2 for now, but it could be 4.1 and come from compiler flags
template <>
struct to_xsimd_arch<16> {
  using type = xsimd::sse2;
};

template <>
struct to_xsimd_arch<32> {
  using type = xsimd::avx;
};

template <uint bytes>
using to_xsimd_arch_t = typename to_xsimd_arch<bytes>::type;

template <class T, uint bytes>
using to_xsimd_batch_t = xsimd::batch<T, to_xsimd_arch_t<bytes>>;

#else

template <class T, uint bytes>
using to_xsimd_batch_t = xsimd::batch<T, bytes / sizeof (T)>;

#endif

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

// SSE or equivalent
using float_x4  = vec<float, 4>;
using double_x2 = vec<double, 2>;
using s8_x16    = vec<s8, 16>;
using u8_x16    = vec<u8, 16>;
using s16_x8    = vec<s16, 8>;
using u16_x8    = vec<u16, 8>;
using s32_x4    = vec<s32, 4>;
using u32_x4    = vec<u32, 4>;
using s64_x2    = vec<s64, 2>;
using u64_x2    = vec<u64, 2>;

// AVX or equivalent.
using float_x8  = vec<float, 8>;
using double_x4 = vec<double, 4>;
using s8_x32    = vec<s8, 32>;
using u8_x32    = vec<u8, 32>;
using s16_x16   = vec<s16, 16>;
using u16_x16   = vec<u16, 16>;
using s32_x8    = vec<s32, 8>;
using u32_x8    = vec<u32, 8>;
using s64_x4    = vec<s64, 4>;
using u64_x4    = vec<u64, 4>;

static constexpr uint sse_bytes = 16;
static constexpr uint avx_bytes = 32;

// just an array rounded to a simd width.
template <class T, size_t N, size_t instr_set_bytes>
using simd_array
  = std::array<T, round_ceil<size_t> (N, (instr_set_bytes / sizeof (T)))>;

template <class V>
using enable_if_vec_t = std::enable_if_t<is_vec_v<V>>;

template <class V>
using enable_if_vec_of_float_point_t
  = std::enable_if_t<is_vec_of_float_type_v<V>>;

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
template <class T>
vec<T, 1> make_vec_x1 (T v)
{
  return {v};
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_to_intrin (V simdvec)
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
template <class V, enable_if_vec_t<V>* = nullptr>
static inline V vec_load (vec_value_type_t<V> const* src)
{
  vec_traits<V>(); // check type validity only
  return *reinterpret_cast<V const*> (src);
}

template <uint N, class T>
static inline vec<T, N> vec_load (T const* src)
{
  return vec_load<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_load (V& dst, vec_value_type_t<V> const* src)
{
  dst = vec_load<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline V vec_load (crange<const vec_value_type_t<V>> src)
{
  constexpr auto traits = vec_traits<V>();
  assert (src.size() >= traits.size);
  return vec_load<V> (src.data());
}

template <uint N, class T>
static inline vec<T, N> vec_load (crange<const T> src)
{
  return vec_load<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_load (V& dst, crange<const vec_value_type_t<V>> src)
{
  dst = vec_load<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline V vec_load_unaligned (vec_value_type_t<V> const* src)
{
  constexpr auto traits = vec_traits<V>();
  using vec_u_type      = typename decltype (traits)::type_u;

  return *reinterpret_cast<vec_u_type const*> (src);
}

template <uint N, class T>
static inline vec<T, N> vec_load_unaligned (T const* src)
{
  return vec_load_unaligned<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_load_unaligned (V& dst, vec_value_type_t<V> const* src)
{
  dst = vec_load_unaligned<V> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline V vec_load_unaligned (crange<const vec_value_type_t<V>> src)
{
  constexpr auto traits = vec_traits<V>();
  assert (src.size() >= traits.size);
  return vec_load_unaligned<V> (src.data());
}

template <uint N, class T>
static inline vec<T, N> vec_load_unaligned (crange<const T> src)
{
  return vec_load_unaligned<vec<T, N>> (src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_load_unaligned (
  V&                                dst,
  crange<const vec_value_type_t<V>> src)
{
  dst = vec_load_unaligned<V> (src);
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_store (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<V*> (dst) = src;
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_store (crange<vec_value_type_t<V>> dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store (dst.data(), src);
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_store_unaligned (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using vec_u_type      = typename decltype (traits)::type_u;
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<vec_u_type*> (dst) = src;
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_store_unaligned (crange<vec_value_type_t<V>> dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store_unaligned (dst.data(), src);
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
static inline void vec_set (V& dst, vec_value_type_t<V> v)
{
  constexpr auto traits = vec_traits<V>();
  dst                   = v - V {}; // convoluted but works.
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline V vec_set (vec_value_type_t<V> v)
{
  V ret;
  vec_set (ret, v);
  return ret;
}

template <uint N, class T>
static inline vec<T, N> vec_set (T v)
{
  vec<T, N> ret;
  vec_set (ret, v);
  return ret;
}
//------------------------------------------------------------------------------
#if 0
// Clang doesn't support "__builtin_shuffle".
template <class V, class IV>
static inline V vec_shuffle (V a, IV mask)
{
  constexpr auto traits      = vec_traits<V>();
  constexpr auto mask_traits = vec_traits<IV>();

  static_assert (traits.size == mask_traits.size, "sizes must match");
  static_assert (
    std::is_integral_v<typename decltype (mask_traits)::value_type>,
    "the mask vector must contain integrals");

  for (uint i = 0; i < traits.size; ++i) {
    assert (mask[i] >= 0 && mask[i] < traits.size);
  }
  return __builtin_shuffle (a, mask);
}

template <class V, class IV>
static inline V vec_shuffle (V a, V b, IV mask)
{
  constexpr auto traits      = vec_traits<V>();
  constexpr auto mask_traits = vec_traits<IV>();

  static_assert (traits.size == mask_traits.size, "sizes must match");
  static_assert (
    std::is_integral_v<typename decltype (mask_traits)::value_type>,
    "the mask vector must contain integrals");

  for (uint i = 0; i < traits.size; ++i) {
    assert (mask[i] >= 0 && mask[i] < (traits.size * 2));
  }

  return __builtin_shuffle (a, b, mask);
}
#endif

#if 1
// This is the best implementation, as "__builtin_shufflevector" requires
// constants as the indexes. This means that a function would require passing
// std::integral_constant and a lot uglyness to just achieve the same that this
// macro. The ifdefed function below doesn't work.
#define vec_shuffle __builtin_shufflevector
#else
template <class V, class... Ts>
static inline V vec_shuffle (V a, V b, Ts... indexes)
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
template <class T, class V>
static inline auto vec_cast (V a)
{
  constexpr auto src_traits = vec_traits<V>();
  using dst_traits = typename decltype (src_traits)::template rebind_traits<T>;

  static_assert (src_traits.size == dst_traits {}.size, "sizes must match");

  return __builtin_convertvector(a, typename dst_traits::type);
}
//------------------------------------------------------------------------------
// SIMD functions
//------------------------------------------------------------------------------
namespace detail {
template <class F1, class F2, class... Ts>
static inline auto call_vec_function_impl (
  F1&& simdf,
  F2&& scalarf,
  Ts&&... args)
{
  using V               = std::common_type_t<Ts...>;
  constexpr auto traits = vec_traits<V>();

  if constexpr (traits.size > 1) {
    // simd, As of now only XSIMD, but the type returned by "simdf" could be
    // detected.
    using batch  = to_xsimd_batch_t<vec_value_type_t<V>, traits.bytes>;
    using intrin = decltype (vec_to_intrin (V {}));

    // xsimd batches can be casted back to intrinsic types
    auto result = static_cast<intrin> (
      simdf (batch {vec_to_intrin (std::forward<Ts> (args))}...));

    return *reinterpret_cast<V*> (&result);
  }
  else {
    // vectors of size == 1, scalar.

    // Workaround. Doing the obvious:
    //  T result = scalarf (args[0]...);
    //
    // gives back:
    //  "error: non-const reference cannot bind to vector element".
    //
    // Vector types are a bit special.
    //
    // Workarounding via casts, as they are simpler than the  alternative I can
    // think of now, which is making copies via a tuple, and then writing an
    // ad-hoc custom function to unpack while accessing the first element of
    // the vector type. Fugly.
    //
    // It is know that the vector is of size 1 and that it has the
    // "__may_alias__" attribute set.
    using T  = vec_value_type_t<V>;
    T result = scalarf ((*((T*) &args))...);
    return make_vec_x1 (result);
  }
}
// NOTE: The argument reversing on these three overloads below could be made
// generic with some metaprogramming (e.g. boost hana). Going for the explicit
// solution, as the upper bound number of args is known and low
template <class V, class F1, class F2>
static inline auto call_vec_function (V&& x, F1&& simd, F2&& scalar)
{
  return call_vec_function_impl (
    std::forward<F1> (simd), std::forward<F2> (scalar), std::forward<V> (x));
}

template <class V1, class V2, class F1, class F2>
static inline auto call_vec_function (V1&& x, V2&& y, F1&& simd, F2&& scalar)
{
  return call_vec_function_impl (
    std::forward<F1> (simd),
    std::forward<F2> (scalar),
    std::forward<V1> (x),
    std::forward<V2> (y));
}

template <class V1, class V2, class V3, class F1, class F2>
static inline auto call_vec_function (
  V1&& x,
  V2&& y,
  V3&& z,
  F1&& simd,
  F2&& scalar)
{
  return call_vec_function_impl (
    std::forward<F1> (simd),
    std::forward<F2> (scalar),
    std::forward<V1> (x),
    std::forward<V2> (y),
    std::forward<V3> (z));
}

} // namespace detail
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_exp (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::exp (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return exp (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = exp (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_exp2 (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::exp2 (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return exp2 (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = exp2 (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<
    is_vec_of_float_type_v<V1> && is_vec_of_float_type_v<V2>>* = nullptr>
static inline auto vec_pow (V1&& x, V2&& y)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V1> (x),
    std::forward<V2> (y),
    [] (auto&& v1, auto&& v2) {
      return xsimd::pow (
        std::forward<decltype (v1)> (v1), std::forward<decltype (v2)> (v2));
    },
    [] (auto&& v1, auto&& v2) { return pow (v1, v2); });
#else
  constexpr auto                                traits = vec_traits<V1>();
  std::remove_reference_t<std::remove_cv_t<V1>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = pow (x[i], y[i]);
  }
  return ret;
#endif
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_pow (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_pow (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_sqrt (V&& x)
{
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::sqrt (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return sqrt (v); });
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_cbrt (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::cbrt (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return cbrt (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = cbrt (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_log (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::log (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return log (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = log (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_log2 (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::log2 (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return log2 (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = log2 (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_log10 (V&& x)
{
#if XSIMD_DISABLED == 0
  using T = vec_value_type_t<V>;
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) {
      return xsimd::log (std::forward<decltype (v)> (v)) * (T) M_LOG10E;
    },
    [] (auto&& v) { return log10 (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = log10 (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_sin (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::sin (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return sin (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = sin (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_cos (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::cos (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return cos (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = cos (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_tan (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::tan (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return tan (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = tan (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_sinh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::sinh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return sinh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = sinh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_cosh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::cosh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return cosh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = cosh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_tanh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::tanh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return tanh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = tanh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_asin (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::asin (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return asin (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = asin (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_acos (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::acos (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return acos (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = acos (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_atan (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::atan (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return atan (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = atan (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_asinh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::asinh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return asinh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = asinh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_acosh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::acosh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return acosh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = acosh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_atanh (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::atanh (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return atanh (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = atanh (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_atan2 (V&& x, V&& y)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    std::forward<V> (y),
    [] (auto&& v1, auto&& v2) {
      return xsimd::atan2 (
        std::forward<decltype (v1)> (v1), std::forward<decltype (v2)> (v2));
    },
    [] (auto&& v1, auto&& v2) { return atan2 (v1, v2); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = atan2 (x[i], y[i]);
  }
  return ret;
#endif
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_atan2 (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_atan2 (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_atan2 (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_atan2 (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_abs (V&& x)
{
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::abs (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return abs (v); });
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
static inline auto vec_min (V1&& x, V2&& y)
{
  return x < y ? x : y;
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_min (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_min (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_min (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_min (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
static inline auto vec_max (V1&& x, V2&& y)
{
  return x > y ? x : y;
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_max (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_max (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  class V3,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2> && is_vec_v<V3>>* = nullptr>
static inline auto vec_clamp (V1&& x, V2&& min, V3&& max)
{
  return detail::call_vec_function (
    std::forward<V1> (x),
    std::forward<V2> (min),
    std::forward<V3> (max),
    [] (auto&& v1, auto&& v2, auto&& v3) {
      return xsimd::clip (
        std::forward<decltype (v1)> (v1),
        std::forward<decltype (v2)> (v2),
        std::forward<decltype (v3)> (v3));
    },
    [] (auto&& v1, auto&& v2, auto&& v3) { return std::clamp (v1, v2, v3); });
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_clamp (
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
static inline auto vec_clamp (V1&& x, V2&& y, vec_value_type_t<V1> z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    std::forward<V1> (x), std::forward<V2> (y), vec_set<Vv> (z));
}

template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
static inline auto vec_clamp (V1&& x, vec_value_type_t<V1> y, V2&& z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    std::forward<V1> (x), vec_set<Vv> (y), std::forward<V2> (z));
}

template <
  class V1,
  class V2,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2>>* = nullptr>
static inline auto vec_clamp (vec_value_type_t<V1> x, V1&& y, V2&& z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V1>>;
  return vec_clamp (
    vec_set<Vv> (x), std::forward<V1> (y), std::forward<V2> (z));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_clamp (
  vec_value_type_t<V> x,
  vec_value_type_t<V> y,
  V&&                 z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_clamp (vec_set<V> (x), vec_set<Vv> (y), std::forward<V> (z));
}

template <class V, enable_if_vec_t<V>* = nullptr>
static inline auto vec_clamp (
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
static inline auto vec_sgn_no_zero (
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
bool vec_is_all_zeros (V v)
{
  using Vv              = std::remove_reference_t<std::remove_cv_t<V>>;
  constexpr auto traits = vec_traits<Vv>();

  if constexpr (traits.size == 1) {
    // scalar, vector of size 1. Fast case
    using uint_ssv = typename decltype (traits)::same_size_uint_type;
    using uint_ss  = decltype (std::declval<uint_ssv>()[0]);
    return (*reinterpret_cast<uint_ssv*> (&v))[0] == (uint_ss) 0ull;
  }
  else {
    using V64 = vec<u64, sizeof (Vv) / sizeof (u64)>;

    if constexpr (sizeof (Vv) == 16) {
      // hoping that this translates to e.g. PSET on SSE4.1
      V64 u = *reinterpret_cast<V64*> (&v);
      return (u[0] | u[1]) == 0;
    }
    else if constexpr (sizeof (Vv) == 32) {
      V64 u = *reinterpret_cast<V64*> (&v);
      return (u[0] | u[1] | u[2] | u[3]) == 0;
    }
    else {
      static_assert (!std::is_same_v<V, V>, "To be implemented");
      return false;
    }
  }
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_t<V>* = nullptr>
bool vec_is_all_ones (V v)
{
  using Vv              = std::remove_reference_t<std::remove_cv_t<V>>;
  constexpr auto traits = vec_traits<Vv>();

  if constexpr (traits.size == 1) {
    // scalar, vector of size 1. Fast case
    using uint_ssv = typename decltype (traits)::same_size_uint_type;
    using uint_ss  = decltype (std::declval<uint_ssv>()[0]);
    return (*reinterpret_cast<uint_ssv*> (&v))[0] == (uint_ss) -1ull;
  }
  else {
    using V64 = vec<u64, traits.bytes / sizeof (u64)>;

    if constexpr (sizeof (Vv) == 16) {
      V64 u = *reinterpret_cast<V64*> (&v);
      return (u[0] & u[1]) == (u64) -1ull;
    }
    else if constexpr (sizeof (Vv) == 32) {
      V64 u = *reinterpret_cast<V64*> (&v);
      return (u[0] & u[1] & u[2] & u[3]) == (u64) -1ull;
    }
    else {
      static_assert (!std::is_same_v<V, V>, "To be implemented");
      return false;
    }
  }
}
//------------------------------------------------------------------------------
// approximations.
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static inline auto vec_tanh_approx_vaneev (V&& x)
{
  using Vv              = std::remove_reference_t<std::remove_cv_t<V>>;
  constexpr auto traits = vec_traits<Vv>();
  using T               = vec_value_type_t<Vv>;
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=388650&start=45
  auto ax = vec_abs (x);
  auto x2 = x * x;

  // clang-format off
  return (
    x  * ((T) 2.45550750702956 +
    (T) 2.45550750702956 * ax +
    ((T) 0.893229853513558 + (T) 0.821226666969744 * ax) * x2)
    / ((T) 2.44506634652299 + ((T) 2.44506634652299 + x2) *
    vec_abs (x + (T) 0.814642734961073 * x * ax)));
  // clang-format on
}
//------------------------------------------------------------------------------
} // namespace artv
