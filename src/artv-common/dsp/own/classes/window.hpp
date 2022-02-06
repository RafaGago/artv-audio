#pragma once

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V sinc_function (V x)
{
  using T = vec_value_type_t<V>;

  auto pi_x = (T) M_PI * x;

  auto cmp = vec_abs (x) < (T) 0.01;
  if (vec_is_all_ones (cmp)) {
    return vec_cos (pi_x * (T) 0.5) * vec_cos (pi_x * (T) 0.25)
      * vec_cos (pi_x * (T) 0.125);
  }
  else if (vec_is_all_zeros (cmp)) {
    return vec_sin (pi_x) / pi_x;
  }
  else {
    V ret;
    for (uint i = 0; i < vec_traits<V>().size; ++i) {
      if (x[i] < (T) 0.01) {
        ret[i] = cos (pi_x[i] * (T) 0.5) * cos (pi_x[i] * (T) 0.25)
          * cos (pi_x[i] * (T) 0.125);
      }
      else {
        ret[i] = sin (pi_x[i]) / pi_x[i];
      }
    }
    return ret;
  }
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V log_gamma_function (V z)
{
  using T = vec_value_type_t<V>;

  V g = vec_set<V> ((T) -1);

  for (uint i = 0; i < vec_traits<V>().size; ++i) {
    T correction = (T) 0;
    T z_v        = z[i];

    if (z_v <= (T) 0) {
      continue;
    }
    while (z_v < (T) 10) {
      correction += log (z_v);
      z_v += (T) 1;
    }
    static constexpr T ln2pi = gcem::log (2 * M_PI);

    g[i] = (T) 0.5 * (ln2pi - log (z_v));
    g[i]
      += z_v * (log (z_v + ((T) 1 / ((T) 12 * z_v - (T) 0.1 / z_v))) - (T) 1);
    g[i] -= correction;
  }
  return g;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V bessel_i0_function (V z)
{
  constexpr uint n_iter = sizeof (z[0]) * 16; // original is 32 for float
  using T               = vec_value_type_t<V>;

  V y = vec_set<V> ((T) 1);

  for (uint i = 0; i < vec_traits<V>().size; ++i) {
    y[i] = (T) 0;

    if (z[i] == 0) {
      continue;
    }
    for (uint k = 0; k < n_iter; ++k) {
      T t = (T) k * log (0.5 * z[i])
        - log_gamma_function (vec_set<1> ((T) k + (T) 1))[0];
      y[i] += exp ((T) 2 * t);
    }
  }
  return y;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V kaiser (uint n, uint n_count, V beta, V mu)
{
  using T = vec_value_type_t<V>;

  V t = (T) n - ((T) n_count - (T) 1) * (T) 0.5 + mu;
  V r = (T) 2 * t / (T) n_count;
  V a = bessel_i0_function (beta * vec_sqrt (1 - r * r));
  V b = bessel_i0_function (beta);
  return a / b;
}
//------------------------------------------------------------------------------
// Adapted from liquid DSP
// https://github.com/jgaeddert/liquid-dsp

template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static V kaiser_beta_estimate (V att_db)
{
  using T = vec_value_type_t<V>;

  att_db = vec_abs (att_db);

  V beta = vec_set<V> ((T) 0);

  for (uint i = 0; i < vec_traits<V>().size; ++i) {
    if (att_db[i] > (T) 50) {
      beta[i] = (T) 0.1102 * (att_db[i] - (T) 8.7);
    }
    else if (att_db[i] > (T) 50) {
      beta[i] = (T) 0.5842 * pow (att_db[i] - (T) 21, (T) 0.4)
        + (T) 0.07886 * (att_db[i] - (T) 21);
    }
  }
  return beta;
}
//------------------------------------------------------------------------------
// "mu" is a fractional delay in samples, from -0.5 to 0.5
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
static void apply_kaiser_window (crange<V> x, V beta, V mu)
{
  for (uint i = 0; i < x.size(); ++i) {
    x[i] *= kaiser (i, x.size(), beta, mu);
  }
}
//------------------------------------------------------------------------------

} // namespace artv
