#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

struct hardclip_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return std::clamp<T> (x, -1.0, 1.0);
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_clamp (x, (T) -1.0, (T) 1.0);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    auto ax        = abs (x);
    bool unclipped = ax <= 1.0;
    return unclipped ? (x * x) * 0.5 : ax - 0.5;
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int_fn (V x)
  {
    using T      = vec_value_type_t<V>;
    auto ax      = vec_abs (x);
    V    noclip  = x * x * (T) 0.5;
    V    yesclip = ax - (T) 0.5;
    return (ax <= (T) 1.0) ? noclip : yesclip;
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    bool unclipped = abs (x) <= 1.0;
    if (unclipped) {
      return (x * x * x) * (1.0 / 6.0);
    }
    else {
      return ((x * x * 0.5) + (1.0 / 6.0)) * artv::sgn_no_zero (x) - (x * 0.5);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, std::enable_if_t<is_vec_v<V>>* = nullptr>
  static V int2_fn (V x)
  {
    using T = vec_value_type_t<V>;

    V sixth  = vec_set<V> ((T) 1. / 6.);
    V half_x = x * (T) 0.5;

    V noclip  = x * x * x * sixth;
    V yesclip = (x * half_x + sixth) * vec_sgn_no_zero (x) - half_x;

    return (vec_abs (x) <= (T) 1.0) ? noclip : yesclip;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using hardclip_adaa = adaa::waveshaper<hardclip_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
