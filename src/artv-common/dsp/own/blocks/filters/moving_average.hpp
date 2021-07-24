#pragma once

#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

// Quick and dirty. As of now just using the variant with M/L = 2. No need to
// care about a recursive/non-recursive implementation (the recursive variant is
// very efficient but has cummulative errors. Implementing as a template.
//------------------------------------------------------------------------------
namespace artv {

// Aka box/boxcar
template <uint L>
struct moving_average;

template <>
struct moving_average<2> {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum state { z1, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static void lowpass (crange<T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (
    crange<const T>, // coeffs
    crange<T> z, // state
    T         in)
  {
    static_assert (std::is_floating_point<T>::value, "");
    assert (z.size() >= n_states);

    T ret = z[z1];
    z[z1] = in;
    ret += in;
    ret *= (T) 0.5;
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint simd_bytes, class T>
  static simd_reg<T, simd_bytes> tick_multi_aligned (
    crange<const T>,
    crange<T>               z, // state 'z1' 1 to N
    simd_reg<T, simd_bytes> in) // in' 1 to N
  {
    static_assert (std::is_floating_point<T>::value, "");
    using simdreg             = simd_reg<T, simd_bytes>;
    constexpr auto n_builtins = simdreg::size;

    assert (z.size() >= n_builtins * n_states);

    auto    z1_ptr = &z[z1 * n_builtins];
    simdreg ret {z1_ptr, xsimd::aligned_mode {}};
    ret += in;
    ret *= (T) 0.5;

    in.store_aligned (z1_ptr);
    return ret;
  }
  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} // namespace artv
