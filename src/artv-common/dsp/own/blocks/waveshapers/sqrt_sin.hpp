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
    // return tanh (x);
    return (x / sqrt (x * x + T (1.))) + (sin (16 * x)) * (1. / 16.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    V sq, sn;
    sq = x / vec_sqrt (x * x + (T) 1.);
    sn = vec_sin (x * (T) 16.) * (T) (1. / 16.);
    return sq + sn;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    return sqrt (x * x + T (1.)) - (cos (16 * x)) * (1. / 256.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;

    V sq, cs;
    sq = vec_sqrt (x * x + (T) 1.);
    cs = vec_cos (x * (T) 16.) * ((T) 1. / 256.);
    return sq + cs;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    return ((x * sqrt (x * x + T (1.))) + asinh (x)) * T (0.5)
      - (sin (16 * x)) * (1. / 4096.);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    V sq, as, sn;
    sq = vec_sqrt (x * x + (T) 1.) * x;
    as = vec_asinh (x) * (T) 0.5;
    sn = vec_sin (x * (T) 16.) * (T) (1. / 4096.);
    return sq + as - sn;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using sqrt_sin_sigmoid_adaa
  = adaa::waveshaper<sqrt_sin_sigmoid_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
