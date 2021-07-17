#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
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
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return xsimd::clip (x, batch {(T) -1.0}, batch {(T) 1.0});
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
    bool unclipped = abs (x) <= 1.0;
    return unclipped ? (x * x) * 0.5 : x * artv::sgn_no_zero (x) - 0.5;
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;

    batch noclip  = x * x * batch {(T) 0.5};
    batch yesclip = x * sgn_no_zero (x) - batch {(T) 0.5};
    return xsimd::select (xsimd::abs (x) <= batch {(T) 1.0}, noclip, yesclip);
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
  template <class T, size_t N>
  static simd_batch<T, N> int_fn2 (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;

    batch sixth {(T) 1. / 6.};
    batch half_x = x * batch {(T) 0.5};

    batch noclip  = x * x * x * sixth;
    batch yesclip = (x * half_x + sixth) * sgn_no_zero (x) - half_x;

    return xsimd::select (xsimd::abs (x) <= batch {(T) 1.0}, noclip, yesclip);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> sgn_no_zero (simd_batch<T, N> x)
  {
    using batch = simd_batch<T, N>;
    return xsimd::select (x < batch {(T) 0.}, batch {(T) -1.}, batch {(T) 1.});
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using hardclip_waveshaper_adaa = adaa::waveshaper<hardclip_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
