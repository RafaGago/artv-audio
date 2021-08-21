#pragma once

// This was previously using xsimd direcly, but inteface-breaking changes
// motivated this intermediate wrapping (unfortunately, boring). As this code is
// only for Desktop PC and probably desktop ARM, Clang (and maybe GCC) will be
// used. Leveraging the builtin vector extensions instead of OO wrappers.
//
// Unfortunately at the time of writing GCC support for vector types is buggy
// too on C++. This code requires clang 11 or higher:
// https://gcc.gnu.org/bugzilla//show_bug.cgi?id=57572

#define XSIMD_BROKEN_W_FAST_MATH 1
#define XSIMD_NEW_INTERFACE 1

#include <array>
#include <immintrin.h>
#include <type_traits>
#include <xsimd/xsimd.hpp>
#if defined(XSIMD_BROKEN_W_FAST_MATH)
#include <cmath>
#endif

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
  // size attribute.
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
static inline constexpr bool is_vec (T)
{
  return false;
}

template <uint size_of, class T>
static inline constexpr bool is_vec (
  T __attribute__ ((vector_size (size_of), __may_alias__)))
{
  return true;
}

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
static constexpr bool is_floating_point_vec_v
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
//------------------------------------------------------------------------------
// TODO: As of now this is only done for x64, but it can be easily ifdefed
//------------------------------------------------------------------------------
// Notice:
//
// - All these loads, stores and casts could be done with "memcpy" but I
//   lazily went with __may_alias__ and type punning to avoid verifying that
//   memcpy does the right thing.
// - No object orientation, as it would defeat the purpuse of using the
// builtin
//   compiler wrappers.
//------------------------------------------------------------------------------
// FUNCTIONS
//------------------------------------------------------------------------------
// returns a "simd_vector_traits"
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static constexpr auto vec_traits()
{
  return vec_traits_t<V> {};
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static constexpr auto vec_traits (V)
{
  return vec_traits<V>();
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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
  else if constexpr (std::is_integral_v<T> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256i*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, double> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256d*> (&simdvec);
  }
  else if constexpr (std::is_same_v<T, float> && traits.bytes == avx_bytes) {
    return *reinterpret_cast<__m256*> (&simdvec);
  }
  else {
    static_assert (!std::is_same_v<V, V>, "Unkown intrinsic conversion");
    return 0;
  }
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_load (V& dst, vec_value_type_t<V> const* src)
{
  dst = vec_load<V> (src);
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_load (V& dst, crange<const vec_value_type_t<V>> src)
{
  dst = vec_load<V> (src);
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_load_unaligned (V& dst, vec_value_type_t<V> const* src)
{
  dst = vec_load_unaligned<V> (src);
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_load_unaligned (
  V&                                dst,
  crange<const vec_value_type_t<V>> src)
{
  dst = vec_load_unaligned<V> (src);
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_store (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<V*> (dst) = src;
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_store (crange<vec_value_type_t<V>> dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store (dst.data(), src);
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_store_unaligned (vec_value_type_t<V>* dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  using vec_u_type      = typename decltype (traits)::type_u;
  using T               = vec_value_type_t<V>;

  *reinterpret_cast<vec_u_type*> (dst) = src;
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_store_unaligned (crange<vec_value_type_t<V>> dst, V src)
{
  constexpr auto traits = vec_traits<V>();
  assert (dst.size() >= traits.size);
  vec_store_unaligned (dst.data(), src);
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline void vec_set (V& dst, vec_value_type_t<V> v)
{
  constexpr auto traits = vec_traits<V>();
  dst                   = v - V {}; // convoluted but works.
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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
// XSIMD functions
//------------------------------------------------------------------------------
namespace detail {
template <class F, class... Ts>
static inline auto call_xsimd_function_impl (F&& func, Ts&&... args)
{
  using V               = std::common_type_t<Ts...>;
  constexpr auto traits = vec_traits<V>();
  using batch           = to_xsimd_batch_t<vec_value_type_t<V>, traits.bytes>;
  using intrin          = decltype (vec_to_intrin (V {}));

  // xsimd batches can be casted back to intrinsic types
  auto result = static_cast<intrin> (
    func (batch {vec_to_intrin (std::forward<Ts> (args))}...));
  return *reinterpret_cast<V*> (&result);
}
// NOTE: The argument reversing on these three overloads below could be made
// generic with some metaprogramming (e.g. boost hana). Going for the explicit
// solution, as the upper bound number of args is known and low
template <class V, class F>
static inline auto call_xsimd_function (V&& x, F&& func)
{
  return call_xsimd_function_impl (std::forward<F> (func), std::forward<V> (x));
}

template <class V1, class V2, class F>
static inline auto call_xsimd_function (V1&& x, V2&& y, F&& func)
{
  return call_xsimd_function_impl (
    std::forward<F> (func), std::forward<V1> (x), std::forward<V2> (y));
}

template <class V1, class V2, class V3, class F>
static inline auto call_xsimd_function (V1&& x, V2&& y, V3&& z, F&& func)
{
  return call_xsimd_function_impl (
    std::forward<F> (func),
    std::forward<V1> (x),
    std::forward<V2> (y),
    std::forward<V3> (z));
}

} // namespace detail
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_exp (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::exp (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = exp (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  std::enable_if_t<
    is_floating_point_vec_v<V1> && is_floating_point_vec_v<V2>>* = nullptr>
static inline auto vec_pow (V1&& x, V2&& y)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (
    std::forward<V1> (x), std::forward<V2> (y), [] (auto&& v1, auto&& v2) {
      return xsimd::pow (
        std::forward<decltype (v1)> (v1), std::forward<decltype (v2)> (v2));
    });
#else
  constexpr auto traits = vec_traits<V1>();
  auto           ret    = x;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = pow (x[i], y[i]);
  }
  return ret;
#endif
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_pow (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_pow (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_sqrt (V&& x)
{
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::sqrt (std::forward<decltype (v)> (v));
  });
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_log (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::log (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = log (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_sin (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::sin (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = sin (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_cos (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::cos (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = cos (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_tan (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::tan (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = tan (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_sinh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::sinh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = sinh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_cosh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::cosh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = cosh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_tanh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::tanh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = tanh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_asin (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::asin (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = asin (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_acos (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::acos (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = acos (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_atan (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::atan (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = atan (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_asinh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::asinh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = asinh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_acosh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::acosh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = acosh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_atanh (V&& x)
{
#if XSIMD_BROKEN_W_FAST_MATH == 0
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::atanh (std::forward<decltype (v)> (v));
  });
#else
  constexpr auto traits = vec_traits<V>();
  for (uint i = 0; i < traits.size; ++i) {
    x[i] = atanh (x[i]);
  }
  return x;
#endif
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_floating_point_vec_v<V>>* = nullptr>
static inline auto vec_abs (V&& x)
{
  return detail::call_xsimd_function (std::forward<V> (x), [] (auto&& v) {
    return xsimd::abs (std::forward<decltype (v)> (v));
  });
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_min (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_min (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_max (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_max (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_max (vec_set<Vv> (x), std::forward<V> (x));
}
//------------------------------------------------------------------------------
template <
  class V1,
  class V2,
  class V3,
  std::enable_if_t<is_vec_v<V1> && is_vec_v<V2> && is_vec_v<V3>>* = nullptr>
static inline auto vec_clamp (V1&& x, V2&& min, V3&& max)
{
  return detail::call_xsimd_function (
    std::forward<V1> (x),
    std::forward<V2> (min),
    std::forward<V3> (max),
    [] (auto&& v1, auto&& v2, auto&& v3) {
      return xsimd::clip (
        std::forward<decltype (v1)> (v1),
        std::forward<decltype (v2)> (v2),
        std::forward<decltype (v3)> (v3));
    });
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
static inline auto vec_clamp (
  vec_value_type_t<V> x,
  vec_value_type_t<V> y,
  V&&                 z)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_clamp (vec_set<V> (x), vec_set<Vv> (y), std::forward<V> (z));
}

template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
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
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
bool vec_is_all_zeros (V v)
{
  using Vv              = std::remove_reference_t<std::remove_cv_t<V>>;
  constexpr auto traits = vec_traits<Vv>();
  using V64             = vec<u64, sizeof (Vv) / sizeof (u64)>;

  if constexpr (sizeof (Vv) == 16) {
    // hoping that this translates to e.g. PSET on SSE4.1
    V64 u = *reinterpret_cast<V64*> (&v);
    return !(u[0] & u[1]);
  }
  else if constexpr (sizeof (Vv) == 32) {
    V64 u = *reinterpret_cast<V64*> (&v);
    return !(u[0] & u[1] & u[2] & u[3]);
  }
  else {
    static_assert (!std::is_same_v<V, V>, "To be implemented");
    return false;
  }
}
//------------------------------------------------------------------------------
template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
bool vec_is_all_ones (V v)
{
  using Vv              = std::remove_reference_t<std::remove_cv_t<V>>;
  constexpr auto traits = vec_traits<Vv>();
  using V64             = vec<u64, traits.bytes / sizeof (u64)>;

  if constexpr (sizeof (Vv) == 16) {
    V64 u = *reinterpret_cast<V64*> (&v);
    return (u[0] & u[1]) == (u64) -1;
  }
  else if constexpr (sizeof (Vv) == 32) {
    V64 u = *reinterpret_cast<V64*> (&v);
    return !!(u[0] & u[1] & u[2] & u[3]) == (u64) -1;
  }
  else {
    static_assert (!std::is_same_v<V, V>, "To be implemented");
    return false;
  }
}
//------------------------------------------------------------------------------
} // namespace artv
