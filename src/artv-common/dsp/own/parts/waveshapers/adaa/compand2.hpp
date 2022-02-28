#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/parts/waveshapers/adaa.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

// based on x^(N+1/N) and to revert x^(N/N+1). The special case of N=2 is
// optimizable with sqrt, and hence done separately as "compand1".

template <uint Num, uint Den>
struct power_functions {
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V fn (V x)
  {
    using T          = vec_value_type_t<V>;
    constexpr T frac = (T) Num / (T) Den;
    return vec_sgn_no_zero (x) * vec_pow (vec_abs (x), frac);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int_fn (V x)
  {
    using T              = vec_value_type_t<V>;
    constexpr T frac     = (T) Num / (T) Den;
    constexpr T int_frac = frac + (T) 1.;

    return vec_pow (vec_abs (x), int_frac) * (1. / int_frac);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V int2_fn (V x)
  {
    using T               = vec_value_type_t<V>;
    constexpr T frac      = (T) Num / (T) Den;
    constexpr T int_frac  = frac + (T) 1.;
    constexpr T int2_frac = int_frac + (T) 1.;

    return vec_sgn_no_zero (x) * vec_pow (vec_abs (x), int2_frac)
      * (1. / (int_frac * int2_frac));
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint N>
using compand_2a_functions = power_functions<N, N + 1>;

template <uint N>
using compand_2b_functions = power_functions<N + 1, N>;
//------------------------------------------------------------------------------
template <uint order, class N>
using compand_2a_adaa = adaa::waveshaper<compand_2a_functions<N::value>, order>;

template <uint order, class N>
using compand_2b_adaa = adaa::waveshaper<compand_2b_functions<N::value>, order>;
//------------------------------------------------------------------------------
} // namespace artv
