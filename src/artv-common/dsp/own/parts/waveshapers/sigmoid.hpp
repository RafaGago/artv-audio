#pragma once

#include <limits>

#include "artv-common/misc/approximations.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv { namespace sigmoid {
//------------------------------------------------------------------------------
namespace detail {
template <class T>
static constexpr auto epsilon = (T) 1e-20;
}
//------------------------------------------------------------------------------
struct passthrough {
  static constexpr uint n_params = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x)
  {
    return x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit()
  {
    using T = vec_value_type_t<V>;
    return vec_set<V> (std::numeric_limits<T>::infinity());
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return true; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <bool Has_drive = false>
struct rsqrt;
//------------------------------------------------------------------------------
template <>
struct rsqrt<false> {
  static constexpr uint n_params = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x)
  {
    using T = vec_value_type_t<V>;
    return x / vec_sqrt (x * x + (T) 1);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf()
  {
    return vec_set<V> (1);
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return false; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <>
struct rsqrt<true> {
  static constexpr uint n_params = 1;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V drive)
  {
    using T = vec_value_type_t<V>;
    auto nz = (drive != 0.f);
    if (!vec_is_all_zeros (nz)) {
      return nz ? x / vec_sqrt (x * x * drive + (T) 1) : x;
    }
    return x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf (V drive) // drive has to be positive
  {
    using T         = vec_value_type_t<V>;
    constexpr T inf = std::numeric_limits<T>::infinity();
    auto        nz  = (drive != 0.f);
    if (!vec_is_all_zeros (nz)) {
      return nz ? ((T) 1 / vec_sqrt (drive)) : vec_set<V> (inf);
    }
    // TODO: check if this branch is worth for skipping a division
    return vec_set<V> (inf);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static bool is_linear (V drive)
  {
    return vec_is_all_ones (drive == 0.f);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <bool Has_drive = false>
struct tanh;
//------------------------------------------------------------------------------
template <>
struct tanh<false> {
  static constexpr uint n_params = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x)
  {
    return vec_tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf()
  {
    return vec_set<V> (1);
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return false; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <>
struct tanh<true> {
  static constexpr uint n_params = 1;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V drive)
  {
    using T = vec_value_type_t<V>;

    auto nz = (drive != 0.f);
    if (!vec_is_all_zeros (nz)) {
      // avoid 0 div. epsilon is small enough for the audio range to be
      // mostly unnafected.
      drive += detail::epsilon<T>;
      return nz ? (vec_tanh (x * drive) / drive) : x;
    }
    return x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf (V drive)
  {
    using T         = vec_value_type_t<V>;
    constexpr T inf = std::numeric_limits<T>::infinity();
    auto        nz  = (drive != 0.f);
    if (!vec_is_all_zeros (nz)) {
      // avoid 0 div. epsilon is small enough for the audio range to be
      // mostly unnafected.
      drive += detail::epsilon<T>;
      return nz ? ((T) 1 / drive) : vec_set<V> (inf);
    }
    // TODO: check if this branch is worth for skipping a division
    return vec_set<V> (inf);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static bool is_linear (V drive)
  {
    using T = vec_value_type_t<V>;
    return vec_is_all_ones (drive == (T) 0);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <bool Has_drive = false>
struct tanh_vaneev;
//------------------------------------------------------------------------------
template <>
struct tanh_vaneev<false> {
  static constexpr uint n_params = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x)
  {
    using T = vec_value_type_t<V>;
    // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=388650&start=45
    auto ax = vec_abs (x);
    auto x2 = x * x;

    V num = x
      * ((T) 2.45550750702956 + (T) 2.45550750702956 * ax
         + ((T) 0.893229853513558 + (T) 0.821226666969744 * ax) * x2);

    V den = (T) 2.44506634652299
      + ((T) 2.44506634652299 + x2)
        * vec_abs (x + (T) 0.814642734961073 * x * ax);
    return num / den;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf()
  {
    return vec_set<V> (1.00808198701850);
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return false; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <>
struct tanh_vaneev<true> {
  static constexpr uint n_params = 1;
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr auto epsilon = detail::epsilon<T>;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V drive)
  {
    using T = vec_value_type_t<V>;
    auto nz = (drive != 0.f);
    if (!vec_is_all_zeros (nz)) {
      // tanh(x*h)/h
      // avoid 0 div. epsilon is small enough for the audio range to be
      // mostly unnafected.
      drive += detail::epsilon<T>;
      return nz ? tanh_vaneev<false> (x * drive) / drive : x;
    }
    // TODO: check if this branch is worth
    return x;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf (V drive) // drive has to be positive
  {
    using T         = vec_value_type_t<V>;
    constexpr T inf = std::numeric_limits<T>::infinity();
    auto        nz  = drive > epsilon<T>;
    if (!vec_is_all_zeros (nz)) {
      // avoid 0 div. epsilon is small enough for the audio range to be
      // mostly unnafected.
      drive += detail::epsilon<T>;
      return nz ? ((T) 1.00808198701850 / drive) : vec_set<V> (inf);
    }
    // TODO: check if this branch is worth for trading a division
    return vec_set<V> (inf);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static bool is_linear (V drive)
  {
    using T = vec_value_type_t<V>;
    return vec_is_all_ones (drive == (T) 0);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// based on the thread
// https://www.kvraudio.com/forum/viewtopic.php?t=521377
// This one has the advantage of being well-behaved with zero drive, as there is
// no zero division
//------------------------------------------------------------------------------
template <uint N_terms = 3, bool Has_drive = false>
struct tanh_mystran;
//------------------------------------------------------------------------------
template <uint N_terms>
struct tanh_mystran<N_terms, false> {
  static constexpr uint n_params = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x)
  {
    return sigmoid::rsqrt<false>::tick (
      taylor::get<taylor::sinh_zero, N_terms> (x));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf()
  {
    return sigmoid::rsqrt<false>::limit_inf<V>();
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return false; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint N_terms>
struct tanh_mystran<N_terms, true> {
  static constexpr uint n_params = 1;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V drive)
  {
    return sigmoid::rsqrt<true>::tick (
      taylor::get<taylor::sinh_zero, N_terms> (x), drive);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf (V drive)
  {
    return sigmoid::rsqrt<true>::limit_inf (drive);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static bool is_linear (V drive)
  {
    using T = vec_value_type_t<V>;
    return vec_is_all_ones (drive == (T) 0);
  }
  //----------------------------------------------------------------------------
};

namespace detail {
// x / sqrt (x^2 + 1), absolute value of the coeffs
// see
// https://www.kvraudio.com/forum/viewtopic.php?p=7334972&sid=a26443d2f2d2434f8fd2decd8c66d238#p7334972
//
// Demo on desmos, showing how it morphs into a hard clipper.
// some more coeffs "https://www.desmos.com/calculator/qmlvfhsrad".
struct mystran_pre_shape {
  template <uint N_terms, class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static constexpr V tick (V x, V hardness)
  {
    using T            = vec_value_type_t<V>;
    constexpr auto tbl = make_array<T> (
      (T) 0.5, // order 3
      (T) 0.375, // order 5
      (T) 0.3125, // order 7, from now on even and odd
      (T) 0.2734375,
      (T) 0.24609375,
      (T) 0.2255859375,
      (T) 0.20947265625,
      (T) 0.196380615234375,
      (T) 0.1854705810546875,
      (T) 0.176197052001953125,
      (T) 0.1681880950927734375,
      (T) 0.1611802577972412109375,
      (T) 0.15498101711273193359375,
      (T) 0.1494459807872772216796875);

    static_assert (N_terms <= tbl.size());

    V ret = x;
    V pwr = x;
    V h   = hardness;
    for (uint i = 0; i < N_terms; ++i) {
      pwr *= x * x;
      ret += pwr * h * tbl[i];
      h *= hardness * hardness;
    };
    return ret;
  }
  //----------------------------------------------------------------------------
};
} // namespace detail

//------------------------------------------------------------------------------
// based on the thread
// https://www.kvraudio.com/forum/viewtopic.php?t=521377
// This one is a sigmoid with variable hardness and drive.
//------------------------------------------------------------------------------
template <uint N_terms = 3, bool Has_drive = false>
struct mystran;
//------------------------------------------------------------------------------
template <uint N_terms>
struct mystran<N_terms, false> {
  static constexpr uint n_params = 1;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V hardness)
  {
    return sigmoid::rsqrt<false>::tick (
      detail::mystran_pre_shape::tick<N_terms> (x, hardness));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf()
  {
    return sigmoid::rsqrt<false>::limit_inf<V>();
  }
  //----------------------------------------------------------------------------
  static bool is_linear() { return false; }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint N_terms>
struct mystran<N_terms, true> {
  static constexpr uint n_params = 2;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (V x, V drive, V hardness)
  {
    return sigmoid::rsqrt<true>::tick (
      detail::mystran_pre_shape::tick<N_terms> (x, hardness), drive);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V limit_inf (V drive, V)
  {
    return sigmoid::rsqrt<true>::limit_inf (drive);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static bool is_linear (V drive, V)
  {
    using T = vec_value_type_t<V>;
    return vec_is_all_ones (drive == (T) 0);
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
}} // namespace artv::sigmoid
