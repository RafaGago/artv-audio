#pragma once

#include <algorithm>
#include <cmath>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

struct hardclip_functions {
  //----------------------------------------------------------------------------
  template <class T>
  static T fn (T x)
  {
    return std::clamp<T> (x, -1.0, 1.0);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int_fn (T x)
  {
    bool unclipped = abs (x) <= 1.0;
    return unclipped ? (x * x) / 2. : x * sgn_no_zero (x) - 0.5;
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T int2_fn (T x)
  {
    bool unclipped = abs (x) <= 1.0;
    if (unclipped) {
      return (x * x * x) / 6.0;
    }
    else {
      return ((x * x / 2.0) + (1.0 / 6.0)) * sgn_no_zero (x) - (x / 2.0);
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using hardclip_waveshaper_adaa = adaa::waveshaper<hardclip_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
