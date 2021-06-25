#pragma once

// Possible optimizations:
//
// -Zero coeff skipping: When setting the FC at FS/4 (2x) every other
//    coefficient is zero. At FC at FS/8 one of every 4 coefficients is
//    zero, etc. Interesting?
//
// - Interleaved delay line. Less coefficient memory loads.
//
// -Handrolled SSE version. Some performance might be gained if some time is
//    spent optimizing this.

#include <cassert>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
namespace detail {

// A delay line that always return a contiguous chunk of previous samples at
// the expense of using the double of the minimum required memory. The delay
// memory is not owned by the class.
template <class T>
class convolution_delay_line {
public:
  //----------------------------------------------------------------------------
  void reset (crange<T> delay_line_mem, uint delay_line_size)
  {
    assert (delay_line_mem.size() >= delay_line_size * 2);
    _head = 0;
    _size = delay_line_size;
    _z    = delay_line_mem.data();
    memset (_z, 0, sizeof _z[0] * _size * 2);
  }
  //----------------------------------------------------------------------------
  void push (T v)
  {
    _head += _head == 0 ? _size : 0;
    _head -= 1;
    _z[_head]         = v;
    _z[_head + _size] = v;
  }
  //----------------------------------------------------------------------------
  crange<const T> samples() const { return {_z + _head, _size}; }
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
template <class T>
class convolution_block {
public:
  using value_type = T;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem) { _h = kernel_mem.data(); }
  //----------------------------------------------------------------------------
  T tick (convolution_delay_line<T> const& dl)
  {
    T    out {};
    auto z = dl.samples();
    for (uint i = 0; i < z.size(); ++i) {
      out += _h[i] * z[i];
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  T const* _h = nullptr; // kernel
};

} // namespace detail
//------------------------------------------------------------------------------
// A convolution/FIR filter. The kernel and delay line memory are  not owned by
// the class.
template <class T>
class convolution {
public:
  using value_type = T;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem, crange<T> delay_line_mem)
  {
    _impl.reset (kernel_mem);
    _dl.reset (delay_line_mem, kernel_mem.size());
  }
  //----------------------------------------------------------------------------
  T tick (T in)
  {
    _dl.push (in);
    return _impl.tick (in, _dl);
  }
  //----------------------------------------------------------------------------
  uint order() const { return _dl.size(); }
  //----------------------------------------------------------------------------
private:
  detail::convolution_block<T>      _impl;
  detail::convolution_delay_line<T> _dl;
};

//------------------------------------------------------------------------------
template <uint FilterOrder, uint Ratio, uint Channels>
class fir_sr_change_traits {
public:
  static constexpr uint ratio    = Ratio;
  static constexpr uint order    = FilterOrder;
  static constexpr uint channels = Channels;
};
//------------------------------------------------------------------------------
// Polyphase implementation of a FIR decimator
template <class... T>
class fir_decimator;

template <class T>
class fir_decimator<T> {
public:
  using value_type = T;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel, uint ratio, uint channels)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio" the unused parts of the kernel will have zeros that don't
    // contribute.
    uint ksize           = round_ceil<uint> (kernel.size(), ratio);
    uint kernel_mem      = ksize;
    uint delay_lines_mem = 2 * ksize * channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _filters.clear();
    _filters.resize (channels);
    _delays.clear();
    _delays.resize (channels);
    _ratio = ratio;
    // On the decimator the inactive part of the filtering kernel sits at the
    // front.
    memcpy (
      _mem.data() + (ksize - kernel.size()),
      kernel.data(),
      sizeof kernel[0] * kernel.size());

    for (uint i = 0; i < channels; ++i) {
      _filters[i].reset (make_crange (&_mem[0], ksize));
      _delays[i].reset (
        make_crange (&_mem[ksize + (ksize * 2 * i)], ksize * 2), ksize);
    }
  }
  //----------------------------------------------------------------------------
  T tick (crange<const T> in, uint channel)
  {
    assert (in.size() == _ratio);
    assert (channel < _filters.size());
    for (uint i = 0; i < _ratio; ++i) {
      _delays[channel].push (in[i]);
    }
    return _filters[channel].tick (_delays[channel]);
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint channels() const { return _filters.size(); }
  //----------------------------------------------------------------------------
  uint order() const
  {
    assert (_filters.size());
    return _filters[0].order();
  }
  //----------------------------------------------------------------------------
private:
  std::vector<detail::convolution_block<T>>      _filters;
  std::vector<detail::convolution_delay_line<T>> _delays;
  std::vector<T>                                 _mem;
  uint                                           _ratio;
};

template <class T, uint Order, uint Ratio, uint Channels>
class fir_decimator<T, fir_sr_change_traits<Order, Ratio, Channels>> {
public:
  static constexpr uint ratio    = Ratio;
  static constexpr uint order    = Order;
  static constexpr uint channels = Channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel) { _impl.reset (kernel, ratio, channels); }
  //----------------------------------------------------------------------------
  T tick (std::array<T, ratio> in, uint channel)
  {
    return _impl.tick (in, channel);
  }
  //----------------------------------------------------------------------------
private:
  fir_decimator<T> _impl;
};

//------------------------------------------------------------------------------
// Polyphase implementation of a FIR interpolator
template <class... T>
class fir_interpolator;

template <class T>
class fir_interpolator<T> {
public:
  using value_type = T;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel, uint ratio, uint channels)
  {
    // To make it polyphase, we ensure that the coefficients are multiples of
    // "ratio" the unused parts of the kernel will have zeros that don't
    // contribute.
    uint ksize           = round_ceil<uint> (kernel.size(), ratio);
    uint subk_size       = ksize / ratio;
    uint kernel_mem      = ksize;
    uint delay_lines_mem = 2 * subk_size * channels;

    _mem.clear();
    _mem.resize (kernel_mem + delay_lines_mem);
    _filters.clear();
    _filters.resize (channels * ratio);
    _delays.clear();
    _delays.resize (channels);
    _ratio = ratio;

    // One delay line is shared between all subfilters.
    for (uint c = 0; c < channels; ++c) {
      _delays[c].reset (
        make_crange (&_mem[kernel_mem + (c * 2 * subk_size)], 2 * subk_size),
        subk_size);
    }
    for (uint subk = 0, offset = 0; subk < ratio; ++subk, offset += subk_size) {
      // polyphase subkernel interleaving.
      // The interpolator has the inactive parts of the filtering kernel at the
      // front.
      for (uint src = subk, dst = 0; src < kernel.size(); src += ratio, ++dst) {
        _mem[offset + dst] = kernel[src];
      }
      // filter initialization
      for (uint c = 0; c < channels; ++c) {
        _filters[(subk * channels) + c].reset (
          make_crange (&_mem[offset], subk_size));
      }
    }
  }
  //----------------------------------------------------------------------------
  void tick (crange<T> out, T in, uint channel)
  {
    assert (out.size() >= _ratio);
    assert (channel < n_channels());

    auto& z = _delays[channel];
    z.push (in);
    for (uint i = 0; i < _ratio; ++i) {
      out[i] = _filters[(n_channels() * i) + channel].tick (z);
    }
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint n_channels() const { return _delays.size(); }
  //----------------------------------------------------------------------------
  uint order() const
  {
    assert (_delays.size());
    return _delays[0].size() * _ratio;
  }
  //----------------------------------------------------------------------------
private:
  std::vector<detail::convolution_block<T>>      _filters;
  std::vector<detail::convolution_delay_line<T>> _delays;
  std::vector<T>                                 _mem;
  uint                                           _ratio;
};

template <class T, uint Order, uint Ratio, uint Channels>
class fir_interpolator<T, fir_sr_change_traits<Order, Ratio, Channels>> {
public:
  static constexpr uint ratio    = Ratio;
  static constexpr uint order    = Order;
  static constexpr uint channels = Channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> kernel) { _impl.reset (kernel, ratio, channels); }
  //----------------------------------------------------------------------------
  std::array<T, ratio> tick (T in, uint channel)
  {
    std::array<T, ratio> ret;
    _impl.tick (ret, in, channel);
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  fir_interpolator<T> _impl;
};
//------------------------------------------------------------------------------

} // namespace artv
