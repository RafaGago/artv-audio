#pragma once

#include <cmath>
#include <type_traits>

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class T>
static T sinc_function (T x)
{
  static_assert (std::is_floating_point_v<T>);
  auto pix = (T) M_PI * x;

  if (abs (x) < (T) 0.01) {
    return cos (pix * (T) 0.5) * cos (pix * (T) 0.25) * cos (pix * (T) 0.125);
  }
  else {
    return sin (pix) / pix;
  }
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class T>
static T log_gamma_function (T z)
{
  static_assert (std::is_floating_point_v<T>);
  if (z <= (T) 0) {
    return (T) -1;
  }

  T correction = (T) 0;

  while (z < (T) 10) {
    correction += log (z);
    z += (T) 1;
  }
  static constexpr T ln2pi = gcem::log (2 * M_PI);

  T g = (T) 0.5 * (ln2pi - log (z));
  g += z * (log (z + ((T) 1 / ((T) 12 * z - (T) 0.1 / z))) - (T) 1);
  g -= correction;
  return g;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class T>
static T bessel_i0_function (T z)
{
  static_assert (std::is_floating_point_v<T>);
  constexpr uint n_iter = sizeof (z) * 16; // original is 32 for float

  if (z == (T) 0) {
    return (T) 1;
  }

  T y = (T) 0;

  for (uint k = 0; k < n_iter; ++k) {
    T t = (T) k * log (0.5 * z) - log_gamma_function ((T) k + (T) 1);
    y += exp ((T) 2 * t);
  }
  return y;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class T>
static T kaiser (uint n, uint n_count, T beta, T mu)
{
  static_assert (std::is_floating_point_v<T>);

  T t  = (T) n - ((T) n_count - (T) 1) * (T) 0.5 + mu;
  T r  = (T) 2 * t / (T) n_count;
  T r2 = r * r;
  T a  = (r2 <= (T) 1) ? bessel_i0_function (beta * sqrt ((T) 1 - r2)) : (T) 1;
  T b  = bessel_i0_function (beta);
  return a / b;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp

template <class T>
static T kaiser_beta_estimate (T att_db)
{
  static_assert (std::is_floating_point_v<T>);

  att_db = fabs (att_db);

  T beta = (T) 0;

  if (att_db > (T) 50) {
    beta = (T) 0.1102 * (att_db - (T) 8.7);
  }
  else if (att_db > (T) 21) {
    beta = (T) 0.5842 * pow (att_db - (T) 21, (T) 0.4)
      + (T) 0.07886 * (att_db - (T) 21);
  }
  return beta;
}
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class T>
static void apply_kaiser_window (crange<T> x, T beta, T mu)
{
  static_assert (std::is_floating_point_v<T>);

  for (uint i = 0; i < x.size(); ++i) {
    x[i] *= kaiser (i, x.size(), beta, mu);
  }
}
//------------------------------------------------------------------------------

} // namespace artv
