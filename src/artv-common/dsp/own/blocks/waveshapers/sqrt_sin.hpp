#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// clang-format on
struct sqrt_sin_sigmoid_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return (x / sqrt (x * x + (T) 1.)) + (sin ((T) 16 * x)) * (T) (1. / 16.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;
    return (x / vec_sqrt (x * x + (T) 1.))
      + (vec_sin ((T) 16 * x)) * (T) (1. / 16.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return sqrt (x * x + T (1.)) - (cos ((T) 16 * x)) * (T) (1. / 256.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;
    return vec_sqrt (x * x + T (1.)) - (vec_cos ((T) 16 * x)) * (T) (1. / 256.);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return ((x * sqrt (x * x + (T) 1.)) + asinh (x)) * (T) 0.5
      - (sin ((T) 16 * x)) * (T) (1. / 4096.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;
    return ((x * vec_sqrt (x * x + (T) 1.)) + vec_asinh (x)) * (T) 0.5
      - (vec_sin ((T) 16 * x)) * (T) (1. / 4096.);
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_sin_sigmoid_adaa
  = adaa::waveshaper<sqrt_sin_sigmoid_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
