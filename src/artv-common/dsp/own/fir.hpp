#pragma once

// Possible optimizations:
//
// -Zero coeff skipping: When setting the FC at FS/4 (2x) every other
//    coefficient is zero. At FC at FS/8 one of every 4 coefficients is
//    zero, etc. Interesting?
//
// -Handrolled SSE version. Some performance might be gained if some time is
//    spent optimizing this.
//
// - An up/downsampler class that reuses the filtering kernels and allocates all
//   memory in one chunk.

#include <cassert>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
namespace detail {

// A interleaved delay line that always return a contiguous chunk of previous
// samples at the expense of using the double of the minimum required memory.
// The delay memory is not owned by the class.
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
    _head += _head == 0 ? (_size * n_channels) : 0;
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
// not owned by the class.
template <class T, uint channels = 1>
class convolution_block {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem) { _h = kernel_mem.data(); }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (
    convolution_delay_line<T, n_channels> const& dl)
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
// A convolution/FIR filter. The kernel and delay line memory are not owned by
// the class.
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
// Naive-ish polyphase implementation of a FIR decimator
template <class T, uint channels = 1>
class fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel, uint ratio)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio" the unused parts of the kernel will have zeros that don't
    // contribute.
    uint ksize           = round_ceil<uint> (kernel.size(), ratio);
    uint kernel_mem      = ksize;
    uint delay_lines_mem = 2 * ksize * n_channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _ratio = ratio;
    // On the decimator the inactive parts of the filtering kernel sits at the
    // front.
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
    for (uint i = 0; i < _ratio; ++i) {
      std::array<T, n_channels> interleaved_in;
      for (uint c = 0; c < n_channels; ++c) {
        interleaved_in[c] = in[c][i];
      }
      _delay.push (interleaved_in);
    }
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
// Naive-ish polyphase implementation of a polyphase FIR interpolator
template <class T, uint channels = 1>
class fir_interpolator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel, uint ratio)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio" the unused parts of the kernel will have zeros that don't
    // contribute.
    uint ksize           = round_ceil<uint> (kernel.size(), ratio);
    uint subk_size       = ksize / ratio;
    uint kernel_mem      = ksize;
    uint delay_lines_mem = 2 * subk_size * n_channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _filters.clear();
    _filters.resize (ratio);
    _ratio = ratio;

    // One delay line is shared between all subfilters.
    _delay.reset (make_crange (&_mem[kernel_mem], delay_lines_mem), subk_size);

    for (uint subk = 0, offset = 0; subk < ratio; ++subk, offset += subk_size) {
      // polyphase subkernel interleaving.
      // The interpolator has the inactive parts of the filtering kernel at the
      // back.
      for (uint src = subk, dst = 0; src < kernel.size(); src += ratio, ++dst) {
        _mem[offset + dst] = kernel[src];
      }
      // filter initialization
      _filters[subk].reset (make_crange (&_mem[offset], subk_size));
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
    for (uint i = 0; i < _ratio; ++i) {
      auto smpl_arr = _filters[i].tick (_delay);
      for (uint c = 0; c < channels; ++c) {
        smpl_arr[c] *= (T) _ratio; // gain loss compensation
        out[c][i] = smpl_arr[c];
      }
    }
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _filters.size(); }
  //----------------------------------------------------------------------------
  uint order() const { return _delay.size() * _ratio; }
  //----------------------------------------------------------------------------
private:
  std::vector<detail::convolution_block<T, n_channels>> _filters;
  detail::convolution_delay_line<T, n_channels>         _delay;
  std::vector<T>                                        _mem;
  uint                                                  _ratio;
};
//------------------------------------------------------------------------------

} // namespace artv
