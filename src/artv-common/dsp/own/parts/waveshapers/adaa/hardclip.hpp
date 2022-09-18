#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

struct hardclip_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    using T = vec_value_type_t<V>;

    return vec_clamp (x, (T) -1.0, (T) 1.0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T      = vec_value_type_t<V>;
    auto ax      = vec_abs (x);
    V    noclip  = x * x * (T) 0.5;
    V    yesclip = ax - (T) 0.5;
    return (ax <= (T) 1.0) ? noclip : yesclip;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
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
