#pragma once

#include "artv-common/misc/util.hpp"

namespace artv {

// Seen on Saike TanhAA saturator
struct allpass_interpolator {
  //----------------------------------------------------------------------------
  enum coeffs { nu, n_coeffs };
  enum state { y0, d0, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void init (crange<T> c, T frac)
  {
    c[nu] = (1.0 - frac) / (1.0 + frac);
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void repair_unsmoothable_coeffs (crange<T>, crange<const T>)
  {}
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr T tick (
    crange<const T> c, // coeffs
    crange<T>       z, // state
    T               in)
  {
    z[y0] = c[nu] * in + z[d0] - c[nu] * z[y0];
    z[d0] = in;
    return z[y0];
  }
  //----------------------------------------------------------------------------
};

} // namespace artv
