#pragma once

#include <cassert>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

#error "UNTESTED! UNOPTIMIZED! DRAFT!"

#if 0
// Taken from "antti" on KVR. Converts from linear phase FIR to Min phase.
function minfir = fminphase (fir)
  x2 = [fir zeros(1, length(fir)*15)];
  [y, ym] = rceps(x2);
  minfir = ym(1:length(fir));
endfunction
#endif
//------------------------------------------------------------------------------
// Naive but cache friendly Polyphase FIR Downsampler
template <class T, uint Order, uint N, uint Channels = 1>
class fir_decimator {
public:
  using value_type               = T;
  static constexpr uint order    = Order;
  static constexpr uint n        = N;
  static constexpr uint channels = Channels;

  static_assert (Order % N == 0, "Not polyphase");
  //----------------------------------------------------------------------------
  void init (contiguous_range<const T> coeffs)
  {
    assert (coeffs.size() <= Order && "Not enough taps");
    // "get_sample" will store the samples on the delay line in the order the
    // N sample batches are passed (newest last). Reorder the coefficients to
    // match the delay line ordering.
    for (uint i = 0; i < coeffs.size(); i += N) {
      for (uint j = 0; j < N; ++j) {
        _h[i + ((N - 1) - j)] = coeffs[i + j];
      }
    }
    for (uint i = coeffs.size(); i < order; ++i) {
      _h[i] = decltype (_h[i]) {};
    }
    reset (0.);
  }
  //----------------------------------------------------------------------------
  void reset (float) { _z = decltype (_z) {}; }
  //----------------------------------------------------------------------------
  inline std::array<T, channels> process_sample (
    contiguous_range<contiguous_range<const T>> in)
  {
    assert (in.size() == channels && "in must match the number of channels");
    for (uint ch = 0; ch < in.size(); ++ch) {
      assert (
        in[ch].size() == N && "n samples must match the downsampling factor");
    }
    // place all the samples on the queue. There is only one array with the
    // delay lines of each channel interleaved, so it is more cache-friendly.
    _z_head -= N * channels;
    _z_head %= _z.size();
    for (uint ch = 0; ch < channels; ++ch) {
      memcpy (
        &_z[_z_head + (ch * N)],
        in[ch].data(),
        in[ch].size() * sizeof in[0][0]);
    }

    // process the filter only once
    std::array<T, channels> out {};
    for (uint i = 0; i < _h.size(); ++i) {
      uint z_pos_base = (_z_head + i) % _z.size();
      for (uint ch = 0; ch < channels; ++ch) {
        out[ch] += _h[i] * _z[z_pos_base + (ch * N)];
      }
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  size_t                          _z_head = 0;
  std::array<T, order * channels> _z; // delay line(s)
  std::array<T, order>            _h;
};

template <uint N, uint Channels, class T, uint C>
auto make_fir_decimator (std::array<T, C> const& coeffs)
{
  fir_decimator<T, div_ceil (C, N) * N, N, Channels> r {};
  r.init (coeffs);
  return r;
}
//------------------------------------------------------------------------------
// Naive but cache friendly Polyphase FIR Upsampler
template <class T, uint Order, uint N, uint Channels>
class fir_interpolator {
public:
  using value_type               = T;
  static constexpr uint order    = Order;
  static constexpr uint n        = N;
  static constexpr uint channels = Channels;

  static_assert (Order % N == 0, "Not polyphase");
  //----------------------------------------------------------------------------
  void init (contiguous_range<const T> coeffs)
  {
    assert (coeffs.size() <= order && "Not enough taps");
    memcpy (_h.data(), coeffs.data(), coeffs.size() * sizeof coeffs[0]);
    for (uint i = coeffs.size(); i < order; ++i) {
      _h[i] = decltype (_h[i]) {};
    }
    reset (0.);
  }
  //----------------------------------------------------------------------------
  void reset (float) { _z = decltype (_z) {}; }
  //----------------------------------------------------------------------------
  std::array<std::array<T, N>, channels> process_sample (
    contiguous_range<const T> in)
  {
    assert (in.size() == channels && "in must match the number of channels");

    // place all the samples on the queue. There is only one array with the
    // delay lines of each channel interleaved, so it is more cache-friendly.
    _z_head -= channels;
    _z_head %= _z.size();
    for (uint ch = 0; ch < channels; ++ch) {
      _z[_z_head + ch] = in;
    }

    // receive many samples on one pass of the full filter
    std::array<std::array<T, N>, channels> out {};
    for (uint i = 0; i < order; ++i) {
      uint z_pos_base = (_z_head + i) % _z.size();
      for (uint ch = 0; ch < channels; ++ch) {
        out[ch][i % N] += _h[i] * _z[z_pos_base + ch];
      }
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  size_t                          _z_head = 0;
  std::array<T, order * channels> _z; // delay line(s)
  std::array<T, order>            _h;
};

template <uint N, uint Channels, class T, uint C>
auto make_fir_interpolator (std::array<T, C> const& coeffs)
{
  fir_interpolator<T, div_ceil (C, N) * N, N, Channels> r {};
  r.init (coeffs);
  return r;
}
//------------------------------------------------------------------------------
} // namespace artv
