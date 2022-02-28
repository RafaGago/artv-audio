#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

// clang-format on
struct sqrt_sin_sigmoid_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;
    return (x / vec_sqrt (x * x + (T) 1.))
      + (vec_sin ((T) 16 * x)) * (T) (1. / 16.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T = vec_value_type_t<V>;
    return vec_sqrt (x * x + T (1.)) - (vec_cos ((T) 16 * x)) * (T) (1. / 256.);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
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
