#pragma once

#include <algorithm>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/classes/window.hpp"

namespace artv {
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static void get_sinc_lowpass (
  crange<V> dst_kernel,
  V         fc, // 0 to 0.5 (nyquist)
  V         mu)
{
  using T = vec_value_type_t<V>;

  T halfsize = (T) dst_kernel.size() * (T) 0.5;

  for (uint n = 0; n < dst_kernel.size(); ++n) {
    V t           = (T) n - ((T) dst_kernel.size() - (T) 1) / (T) 2 + mu;
    dst_kernel[n] = sinc_function (2 * fc * t);
  }
}
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static void normalize_fir_kernel (crange<V> kernel)
{
  using T = vec_value_type_t<V>;

  V gain
    = (T) 1 / std::accumulate (kernel.begin(), kernel.end(), vec_set<V> (0));
  for (auto& v : kernel) {
    v *= gain;
  }
}
//------------------------------------------------------------------------------
// fc normalized freq, 0 to 0.5 (nyquist)
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static void kaiser_lp_kernel_2 (crange<V> dst_kernel, V fc, V beta, V mu)
{
  get_sinc_lowpass (dst_kernel, fc, mu);
  apply_kaiser_window (dst_kernel, beta, mu);
  normalize_fir_kernel (dst_kernel);
}
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static void kaiser_lp_kernel (
  crange<V> dst_kernel,
  V         fc,
  V         att_db, // positive
  V         mu)
{
  kaiser_lp_kernel_2 (dst_kernel, fc, kaiser_beta_estimate (att_db), mu);
}
//------------------------------------------------------------------------------
} // namespace artv
