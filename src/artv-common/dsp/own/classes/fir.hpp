#pragma once

// Possible (known) optimizations left:
//
// -Handrolled SSE version. Some performance might be gained if some dev time is
//    spent optimizing this. Worth?

#include <cassert>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
namespace detail {

// A interleaved delay line that always returns a contiguous chunk of previous
// samples without needing to reshuffle the memory on each sample. Basically the
// head pointer slides backwards and then jumps to the center. It requires
// double the memory for when the pointer jumps from the head to the center, but
// each insertion is only two (predictable) writes instead of a full reshuffle
// of the delay line.
//
// The last inserted sample comes first on the samples memory chunk. The delay
// memory is not owned by the class, this is so because the most optimal memory
// layout is not known at this abstraction level.
template <class T, uint channels = 1>
class convolution_delay_line {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> delay_line_mem, uint delay_line_size)
  {
    assert (delay_line_mem.size() >= (delay_line_size * 2 * n_channels));
    _head = 0;
    _size = delay_line_size;
    _z    = delay_line_mem.data();
    memset (_z, 0, sizeof _z[0] * _size * 2 * n_channels);
  }
  //----------------------------------------------------------------------------
  void push (std::array<T, n_channels> v)
  {
    _head = _head == 0 ? (_size * n_channels) : _head;
    _head -= n_channels;
    memcpy (_z + _head, v.data(), sizeof v[0] * n_channels);
    memcpy (
      _z + _head + (_size * n_channels), v.data(), sizeof v[0] * n_channels);
  }
  //----------------------------------------------------------------------------
  crange<const T> samples() const { return {_z + _head, _size * n_channels}; }
  //----------------------------------------------------------------------------
  uint size() const { return _size; }
  //----------------------------------------------------------------------------
private:
  T*   _z;
  uint _head = 0;
  uint _size = 0;
};
//------------------------------------------------------------------------------
// A convolution where the delay line is externally controlled. The kernel is
// not owned by the class. This doesn't own memory because the most optimal
// memory layout is not known at this abstraction depth.
template <class T, uint channels = 1>
class convolution_block {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem) { _h = kernel_mem.data(); }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (
    convolution_delay_line<T, n_channels> const& dl,
    uint                                         start_idx = 0)
  {
    std::array<T, n_channels> out {};
    auto                      z = dl.samples();
    for (uint i = 0; i < dl.size(); ++i) {
      for (uint c = 0; c < n_channels; ++c) {
        out[c] += _h[i] * z[(i * n_channels) + c];
      }
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  T const* _h = nullptr; // kernel
};

} // namespace detail
//------------------------------------------------------------------------------
// A non-FFT convolution/FIR filter. The kernel and delay line memory are not
// owned by the class. This is so because the most optimal memory layout is not
// known at this abstraction level.
template <class T, uint channels = 1>
class convolution {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem, crange<T> delay_line_mem)
  {
    _impl.reset (kernel_mem);
    _dl.reset (delay_line_mem, kernel_mem.size());
  }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (std::array<T, n_channels> in)
  {
    _dl.push (in);
    return _impl.tick (in, _dl);
  }
  //----------------------------------------------------------------------------
  uint order() const { return _dl.size(); }
  //----------------------------------------------------------------------------
private:
  detail::convolution_block<T, n_channels>      _impl;
  detail::convolution_delay_line<T, n_channels> _dl;
};
//------------------------------------------------------------------------------
// class for verifying if a filter kernel (set of coefficients) can be used
// on a lth-band filter.
struct lth_band {
  template <class T>
  static bool verify_coeffs (crange<const T> kernel, uint ratio)
  {
    for (uint i = 0; i < kernel.size(); ++i) {
      if (i % ratio) {
        continue;
      }
      if (i != (kernel.size() / 2)) {
        if (kernel[i] != (T) 0) {
          return false;
        }
      }
      else {
        T expected = ((T) 1) / ((T) ratio);
        if (std::abs (kernel[i] - expected) > 0.0001) {
          return false;
        }
      }
    }
    return true;
  }
};
//------------------------------------------------------------------------------
// Polyphase FIR decimator
template <class T, uint channels = 1>
class fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel, uint ratio)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio". This is done by adding zeros coefficients that don't contribute
    // to the response.
    uint ksize      = round_ceil<uint> (kernel.size(), ratio);
    uint kernel_mem = ksize;
    // magic "2" reminder: delay line buffer has to be double the required size
    uint delay_lines_mem = 2 * ksize * n_channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _ratio = ratio;
    // On the decimator the zero coefficients only added for the filter to be
    // polyphase decomposable have to be at the front.
    memcpy (
      _mem.data() + (ksize - kernel.size()),
      kernel.data(),
      sizeof kernel[0] * kernel.size());

    _filter.reset (make_crange (&_mem[0], ksize));
    _delay.reset (make_crange (&_mem[ksize], delay_lines_mem), ksize);
  }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (std::array<crange<const T>, n_channels> in)
  {
    for (auto& r : in) {
      assert (r.size() == _ratio);
    }
    // interleaving channel samples and pushing them to the delay lines (inverts
    // order).
    for (uint i = 0; i < _ratio; ++i) {
      std::array<T, n_channels> interleaved_in;
      for (uint c = 0; c < n_channels; ++c) {
        interleaved_in[c] = in[c][i];
      }
      _delay.push (interleaved_in);
    }
    // filtering
    return _filter.tick (_delay);
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint order() const { return _filter.order(); }
  //----------------------------------------------------------------------------
private:
  detail::convolution_block<T, n_channels>      _filter;
  detail::convolution_delay_line<T, n_channels> _delay;
  std::vector<T>                                _mem;
  uint                                          _ratio;
};
//------------------------------------------------------------------------------
// Polyphase FIR decimator optimized for L-th band decimation (coefficients with
// a regular pattern of zeros).
template <class T, uint channels = 1>
class lth_band_fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel, uint ratio)
  {
    assert (lth_band::verify_coeffs (kernel, ratio));

    uint ksize      = round_ceil<uint> (kernel.size(), ratio);
    uint subk_size  = ksize / ratio;
    uint kernel_mem = ksize - subk_size;
    // magic "2" reminder: delay line buffer has to be double the required size
    uint delay_lines_mem  = 2 * ksize * n_channels;
    uint zero_dl_size     = 2 * subk_size * n_channels;
    uint non_zero_dl_size = delay_lines_mem - zero_dl_size;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);

    // Copy skipping the coefficients multiple of the ratio: the (zeros or
    // the central coefficient).
    //
    // On the decimator the zero coefficients only added for the filter to be
    // polyphase decomposable have to be at the front.
    uint wr_pos = ratio - 1;
    for (uint i = 0; i < kernel.size(); ++i) {
      if ((i % ratio) == 0) {
        continue;
      }
      _mem[wr_pos] = kernel[i];
      ++wr_pos;
    }

    _ratio        = ratio;
    _center_coeff = kernel[kernel.size() / 2];

    _filter.reset (make_crange (&_mem[0], kernel_mem));

    _delays[0].reset (
      make_crange (&_mem[kernel_mem], non_zero_dl_size), ksize - subk_size);
    _delays[1].reset (
      make_crange (&_mem[kernel_mem + non_zero_dl_size], zero_dl_size),
      subk_size);
  }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (std::array<crange<const T>, n_channels> in)
  {
    for (auto& r : in) {
      assert (r.size() == _ratio);
    }
    // interleaving channel samples and pushing them to the delay lines (inverts
    // order).
    for (uint i = 0; i < _ratio; ++i) {
      std::array<T, n_channels> interleaved_in;
      for (uint c = 0; c < n_channels; ++c) {
        interleaved_in[c] = in[c][i];
      }
      _delays[i == 0].push (interleaved_in);
    }
    // filtering branches with non zero coefficients
    uint                      n_trailing_zero_coeffs = _ratio - 1;
    std::array<T, n_channels> ret
      = _filter.tick (_delays[0], n_trailing_zero_coeffs);

    // adding the only non zero coeff at the center of the remaining branch.
    T const* center = _delays[1].samples().data();
    center += (_delays[1].size() / 2) * channels;
    for (uint c = 0; c < n_channels; ++c) {
      ret[c] += center[c] * _center_coeff;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint order() const { return (_filter.order() / (_ratio - 1)) * _ratio; }
  //----------------------------------------------------------------------------
private:
  detail::convolution_block<T, n_channels>                     _filter;
  std::array<detail::convolution_delay_line<T, n_channels>, 2> _delays;
  std::vector<T>                                               _mem;
  T                                                            _center_coeff;
  uint                                                         _ratio;
};
//------------------------------------------------------------------------------
// Polyphase FIR interpolator with selectable L-th band optimization (in case
// the coefficients are suitable for Lth band optimization (skipping zeroes)).
//
// As the interpolator is always decomposed in branches, adding the Lth band
// optimization doesn't add too much bloat to the class, so it is done in place.
template <class T, uint channels = 1>
class fir_interpolator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  // If the filter cutoff frequency is Fs/4 and the number of coefficients minus
  // one is a multiple of 4 then this is Lth-band filter. Lth band filters have
  // 1/ratio of the coefficients as zeros, so they can be skipped. A polyphase
  // FIR interpolator is structured in a way that the skipping the zeroed branch
  // is trivial, so it doesn't have a separate implementation as the decimator
  // does.
  void reset (
    crange<const T> kernel,
    uint            ratio,
    bool            enable_lth_band_optimization = false)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio" the unused parts of the kernel will have zeros that don't
    // contribute to the response.

    _is_lth_band = enable_lth_band_optimization;
    _ratio       = ratio;

    assert (!_is_lth_band || lth_band::verify_coeffs (kernel, ratio));

    uint ksize      = round_ceil<uint> (kernel.size(), ratio);
    uint subk_size  = ksize / ratio;
    uint kernel_mem = ksize;
    kernel_mem -= (_is_lth_band) ? subk_size : 0;
    // magic "2" reminder: delay line buffer has to be double the required size
    uint delay_lines_mem = 2 * subk_size * n_channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _filters.clear();
    _filters.resize (ratio - _is_lth_band);

    // One delay line is shared between all subfilters.
    _delay.reset (make_crange (&_mem[kernel_mem], delay_lines_mem), subk_size);

    for (uint subk = _is_lth_band, offset = 0; subk < ratio;
         ++subk, offset += subk_size) {
      // polyphase subkernel interleaving.
      // The interpolator has the inactive parts of the filtering kernel at the
      // back.
      for (uint src = subk, dst = 0; src < kernel.size(); src += ratio, ++dst) {
        _mem[offset + dst] = kernel[src];
        _mem[offset + dst] *= (T) _ratio; // gain loss compensation
      }
      // filter initialization
      _filters[subk - _is_lth_band].reset (
        make_crange (&_mem[offset], subk_size));
    }
  }
  //----------------------------------------------------------------------------
  void tick (
    std::array<crange<T>, n_channels> out,
    std::array<T, n_channels>         in)
  {
    for (auto& r : out) {
      assert (r.size() >= _ratio);
    }
    _delay.push (in);
    uint is_lth_band = _is_lth_band;
    for (uint i = is_lth_band; i < _ratio; ++i) {
      auto smpl_arr = _filters[i - is_lth_band].tick (_delay, is_lth_band);
      for (uint c = 0; c < channels; ++c) {
        out[c][i] = smpl_arr[c];
      }
    }
    if (is_lth_band) {
      // Handle the central coefficient, as this class is normalizing the
      // kernel to compensate for the energy spread on the upsampling imaging,
      // the center coefficient on a L-th band linear phase lowpass is always 1,
      // so it is a direct copy of the sample, as multiplying by one has no
      // effect.
      const T* center = _delay.samples().data();
      center += (_delay.size() / 2) * channels;
      for (uint c = 0; c < channels; ++c) {
        out[c][0] = center[c];
      }
    }
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint order() const { return _delay.size() * _ratio; }
  //----------------------------------------------------------------------------
private:
  std::vector<detail::convolution_block<T, n_channels>> _filters;
  detail::convolution_delay_line<T, n_channels>         _delay;
  std::vector<T>                                        _mem;
  uint                                                  _ratio;
  bool                                                  _is_lth_band;
};
//------------------------------------------------------------------------------

} // namespace artv
