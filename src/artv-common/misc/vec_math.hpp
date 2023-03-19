

#pragma once

// math.h fuctions to help working with vectors. Using XSIMD or scalars as of
// now
//
// At some point XSIMD might be dropped, as Clang supports libmvec and SVML, so
// it is a matter of trying if it works for Windows and doing the refactoring
// chores. Given my (good)interactions with XSIMD devs, it seems it is not to
// rely on with --ffast-math enabled, as it is not their use case. I'd trust
// more something built-in on the compiler.

#define XSIMD_DISABLED      1
#define XSIMD_NEW_INTERFACE 1

#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <xsimd/xsimd.hpp>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/hana.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_util.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

namespace detail {

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
// SIMD math functions
//------------------------------------------------------------------------------
namespace detail {
template <class F1, class F2, class... Ts>
inline auto call_vec_function_impl (F1&& simdf, F2&& scalarf, Ts&&... args)
{
  using V               = std::common_type_t<Ts...>;
  using T               = vec_value_type_t<V>;
  constexpr auto traits = vec_traits<V>();

  if constexpr (traits.size > 1) {
    if constexpr (traits.bytes <= vec_max_bytes) {
      // simd, As of now only XSIMD, but the type returned by "simdf" could be
      // detected.
      using batch  = to_xsimd_batch_t<T, traits.bytes>;
      using intrin = decltype (vec_to_intrin (V {}));

      // xsimd batches can be casted back to intrinsic types
      auto result = static_cast<intrin> (
        simdf (batch {vec_to_intrin (std::forward<Ts> (args))}...));
      // vector types are __may_alias__, so they can be casted back from
      // intrinsic types
      return *reinterpret_cast<V*> (&result);
    }
    else {
      // The virtual vector simd type is bigger than the hardware one, so doing
      // it in batches of the biggest SIMD instruction available and rejoining.
      V ret;

      static constexpr uint n_elems = vec_max_bytes / sizeof (T);
      using V_small                 = vec<T, n_elems>;
      static constexpr uint n_parts = sizeof (V) / sizeof (V_small);
      using batch                   = to_xsimd_batch_t<T, vec_max_bytes>;
      using intrin                  = decltype (vec_to_intrin (V_small {}));

      for (uint i = 0; i < n_parts; ++i) {
        // xsimd batches can be casted back to intrinsic types
        auto result = static_cast<intrin> (simdf (batch {vec_to_intrin (
          vec_slice<n_elems> (std::forward<Ts> (args), i * n_elems))}...));
        // vector types are __may_alias__, so they can be casted back from
        // intrinsic types
        vec_cp_slice (ret, *reinterpret_cast<V_small*> (&result), i * n_elems);
      }
      return ret;
    }
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
    T result = scalarf ((*((T*) &args))...);
    return make_vec (result);
  }
}
// NOTE: The argument reversing on these three overloads below could be made
// generic with some metaprogramming (e.g. boost hana). Going for the explicit
// solution, as the upper bound number of args is known and low
template <class V, class F1, class F2>
inline auto call_vec_function (V&& x, F1&& simd, F2&& scalar)
{
  return call_vec_function_impl (
    std::forward<F1> (simd), std::forward<F2> (scalar), std::forward<V> (x));
}

template <class V1, class V2, class F1, class F2>
inline auto call_vec_function (V1&& x, V2&& y, F1&& simd, F2&& scalar)
{
  return call_vec_function_impl (
    std::forward<F1> (simd),
    std::forward<F2> (scalar),
    std::forward<V1> (x),
    std::forward<V2> (y));
}

template <class V1, class V2, class V3, class F1, class F2>
inline auto call_vec_function (V1&& x, V2&& y, V3&& z, F1&& simd, F2&& scalar)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_exp (V&& x)
{
#if 0 // TBI
  return __builtin_elementwise_exp (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_exp2 (V&& x)
{
#if 0 // TBI
  return __builtin_elementwise_exp2 (x);
#elif XSIMD_DISABLED == 0
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
  std::enable_if_t<is_floatpt_vec_v<V1> && is_floatpt_vec_v<V2>>* = nullptr>
inline auto vec_pow (V1&& x, V2&& y)
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
inline auto vec_pow (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
inline auto vec_pow (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_pow (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_sqrt (V&& x)
{
#if XSIMD_DISABLED == 0
  return detail::call_vec_function (
    std::forward<V> (x),
    [] (auto&& v) { return xsimd::sqrt (std::forward<decltype (v)> (v)); },
    [] (auto&& v) { return sqrt (v); });
#else
  constexpr auto                               traits = vec_traits<V>();
  std::remove_reference_t<std::remove_cv_t<V>> ret;
  for (uint i = 0; i < traits.size; ++i) {
    ret[i] = sqrt (x[i]);
  }
  return ret;
#endif
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_cbrt (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_log (V&& x)
{
#if 0 // TBI
  return __builtin_elementwise_log (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_log2 (V&& x)
{
#if 0 // TBI
  return __builtin_elementwise_log2 (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_log10 (V&& x)
{
#if 0 // TBI
  return __builtin_elementwise_log10 (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_sin (V&& x)
{
#if 1
  return __builtin_elementwise_sin (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_cos (V&& x)
{
#if 1
  return __builtin_elementwise_cos (x);
#elif XSIMD_DISABLED == 0
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_tan (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_sinh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_cosh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_tanh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_asin (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_acos (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_atan (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_asinh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_acosh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_atanh (V&& x)
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
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline auto vec_atan2 (V&& x, V&& y)
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
inline auto vec_atan2 (V&& x, vec_value_type_t<V> y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_atan2 (std::forward<V> (x), vec_set<Vv> (y));
}

template <class V, enable_if_vec_t<V>* = nullptr>
inline auto vec_atan2 (vec_value_type_t<V> x, V&& y)
{
  using Vv = std::remove_reference_t<std::remove_cv_t<V>>;
  return vec_atan2 (vec_set<Vv> (x), std::forward<V> (y));
}
//------------------------------------------------------------------------------
// common functions
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline V vec_db_to_gain (V db)
{
  return vec_exp (db * vec_set<V> (M_LN10 / 20.));
}
//------------------------------------------------------------------------------
template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
inline V vec_gain_to_db (V gain)
{
  using T           = vec_value_type_t<V>;
  constexpr T small = -(std::numeric_limits<T>::min() * 1e6f);

  V absv = vec_abs (gain) + small; // add ultra small signal to avoid log(0);
  return vec_log (absv) * (20. / M_LN10);
}
//------------------------------------------------------------------------------

} // namespace artv
