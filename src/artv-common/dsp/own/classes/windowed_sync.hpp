#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/classes/window.hpp"

namespace artv {
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class T>
static void get_sinc_lowpass (
  crange<T> dst_kernel,
  T         fc, // 0 to 0.5 (nyquist)
  T         mu)
{
  static_assert (std::is_floating_point_v<T>);

  for (uint n = 0; n < dst_kernel.size(); ++n) {
    T t           = (T) n - ((T) dst_kernel.size() - (T) 1) * (T) 0.5 + mu;
    dst_kernel[n] = sinc_function ((T) 2 * fc * t);
  }
}
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class T>
static void fir_kernel_normalize (crange<T> kernel)
{
  static_assert (std::is_floating_point_v<T>);

  T gain = (T) 1 / std::accumulate (kernel.begin(), kernel.end(), (T) 0);
  for (auto& v : kernel) {
    v *= gain;
  }
}
//------------------------------------------------------------------------------
// fc normalized freq, 0 to 0.5 (nyquist)
template <class T>
static void kaiser_lp_kernel_2 (crange<T> dst_kernel, T fc, T beta, T mu)
{
  static_assert (std::is_floating_point_v<T>);

  get_sinc_lowpass (dst_kernel, fc, mu);
  apply_kaiser_window (dst_kernel, beta, mu);
  fir_kernel_normalize (dst_kernel);
}
//------------------------------------------------------------------------------
template <class T>
static void kaiser_lp_kernel (
  crange<T> dst_kernel,
  T         fc,
  T         att_db, // positive
  T         mu)
{
  static_assert (std::is_floating_point_v<T>);

  kaiser_lp_kernel_2 (dst_kernel, fc, kaiser_beta_estimate (att_db), mu);
}
//------------------------------------------------------------------------------
} // namespace artv
