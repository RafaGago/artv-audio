#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <type_traits>

#include "artv-common/misc/bits.hpp"
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
// log cepstrum method (Mystran's post):
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=556692&start=45
template <class T, class U, class FFT>
static void fir_kernel_to_minphase (
  crange<T> kernel,
  FFT&      fft,
  crange<U> fft_work_buffer)
{
  static_assert (std::is_floating_point_v<T>);
  static_assert (std::is_floating_point_v<U>);
  static_assert (std::is_same_v<U, typename FFT::value_type>);

  auto       fftwb    = fft_work_buffer;
  const uint fftsize  = fftwb.size() / 2;
  const U    fftscale = (U) 1 / (U) fftsize;

  assert (is_pow2 (fftsize));
  assert (fftsize > (kernel.size() * 8));

  crange_memset (fftwb, 0);

  // regular fft of the impulse
  for (uint i = 0; i < kernel.size(); ++i) {
    fftwb[i * 2] = (U) kernel[i];
  }
  fft.forward_transform (fftwb.data(), fftsize);

  // replace real part by log of magnitude, leave imaginary part to 0.
  for (uint i = 0; i < fftwb.size(); i += 2) {
    U magnitude  = std::abs (std::complex<U> {fftwb[i], fftwb[i + 1]});
    fftwb[i]     = (U) log (std::max (magnitude, 1e-24));
    fftwb[i + 1] = (U) 0;
  }
  // FFT back from log spectrum to log cepstrum (and rescaling)
  fft.backward_transform (fftwb.data(), fftsize);
  for (auto& v : fftwb) {
    v *= fftscale;
  }
  //  multiply elements in range [1, N/2-1] by 2 and zero the elements in range
  // [N/2+1, N-1], bins 0 and N/2 are left "as-is".
  for (uint i = 2; i < fftsize; i += 2) {
    fftwb[i] *= (U) 2;
    fftwb[i + 1] = (U) 0; // should be real already...
  }
  memset (&fftwb[fftsize + 2], 0, sizeof (U) * (fftsize - 2));
  // FFT again
  fft.forward_transform (fftwb.data(), fftsize);
  // complex "exp" of each bin
  for (uint i = 0; i < fftwb.size(); i += 2) {
    auto expv = std::exp (std::complex<U> {fftwb[i], fftwb[i + 1]});
    fftwb[i] *= expv.real();
    fftwb[i + 1] = expv.imag();
  }
  // IFFT and truncating back to the kernel.
  fft.backward_transform (fftwb.data(), fftsize);
  for (auto& v : fftwb) {
    v *= fftscale;
  }

  for (uint i = 0; i < kernel.size(); ++i) {
    kernel[i] = (T) fftwb[i * 2];
  }
}
//------------------------------------------------------------------------------
// https://www.dsprelated.com/freebooks/filters/Numerical_Computation_Group_Delay.html
template <class T, class U, class FFT>
static void fir_kernel_group_delay (
  crange<T> kernel,
  T         fc, // normalized, from 0 to 0.5
  FFT&      fft,
  crange<U> fft_work_buffer1,
  crange<U> fft_work_buffer2)
{
  static_assert (std::is_floating_point_v<T>);
  static_assert (std::is_floating_point_v<U>);
  static_assert (std::is_same_v<U, typename FFT::value_type>);

  auto b  = fft_work_buffer1;
  auto br = fft_work_buffer2;

  const uint fftsize = b.size() / 2;

  assert (br.size() == b.size());
  assert (is_pow2 (fftsize));

  crange_memset (b, 0);
  crange_memset (br, 0);

  for (uint i = 0; i < kernel.size(); ++i) {
    b[i * 2]  = (U) kernel[i];
    br[i * 2] = (U) kernel[i] * (U) i; // ramped polynomial
  }

  fft.forward_transform (b.data(), fftsize);
  fft.forward_transform (br.data(), fftsize);

  fft.forward_permute (b.data(), fftsize);
  fft.forward_permute (br.data(), fftsize);

  uint bucket = (uint) (((T) 2 * fc * (T) 2 * (T) (fftsize - 1)) + (T) 0.5);
  auto div    = std::complex<double> {br[bucket], br[bucket + 1]};
  div /= std::complex<double> {b[bucket], b[bucket + 1]};
  return (T) real (div);
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