#pragma once

// Possible (known) optimizations left:
//
// -Handrolled SSE version. Some performance might be gained if some dev time is
//    spent optimizing this. Worth?

#include <cassert>

#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
namespace detail {

// A delay line that always returns a contiguous chunk of previous samples
// without needing to reshuffle the memory on each sample. Basically the head
// pointer slides backwards and then jumps to the center. It requires double the
// memory for when the pointer jumps from the head to the center, but each
// insertion is only two (predictable) writes per channel instead of a full
// reshuffle of the delay line.
//
// The last inserted sample comes first on the samples memory chunk. The delay
// memory is not owned by the class, this is so because the most optimal memory
// layout is not known at this abstraction level.
//
template <class T, uint channels = 1>
class convolution_delay_line {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> delay_line_mem, uint channel_size)
  {
    assert (delay_line_mem.size() >= (channel_size * n_channels * 2));
    _head = 0;
    _size = channel_size;
    _z    = delay_line_mem.data();
    memset (_z, 0, sizeof _z[0] * _size * 2 * n_channels);
  }
  //----------------------------------------------------------------------------
  void push (std::array<T, n_channels> v)
  {
    _head = _head == 0 ? _size : _head;
    --_head;

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head = _z + _head + (chnl * _size * 2);
      *chnl_head = *(chnl_head + _size) = v[chnl];
    });
  }
  //----------------------------------------------------------------------------
  std::array<const T * artv_restrict, n_channels> samples() const
  {
    std::array<const T * artv_restrict, n_channels> ret;

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      ret[chnl] = _z + _head + (chnl * _size * 2);
    });
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr uint size() const { return _size; }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  T*   _z;
  uint _head = 0;
  uint _size = 0;
};
//------------------------------------------------------------------------------
// A convolution where the delay line is externally controlled. The kernel is
// not owned by the class. This doesn't own memory because the most optimal
// memory layout is not known at this abstraction depth.
//
// the "alignment" template parameter is an external guarantee not enforced by
// the class. It specifies that the memory for both the kernel will be aligned
// to this value.
template <class T, uint channels = 1, uint alignment = alignof (T)>
class convolution_block {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel_mem) { _h = kernel_mem.data(); }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (
    convolution_delay_line<T, n_channels> const& dl,
    uint                                         trailing_zero_coeffs_hint = 0)
  {
    // The pointer of the kernel should be aligned to "alignment".
    alignas (alignment) std::array<T, n_channels> out {};

    auto                   z = dl.samples();
    T const* artv_restrict h = assume_aligned_hint<alignment> (_h);
    assert (is_aligned_to (alignment, h));

    uint size = dl.size();

    // TODO: More hints for clang?
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < size; ++i) {
      mp_foreach_idx<n_channels> ([&] (auto chnl) {
        out[chnl] += h[i] * z[chnl][i];
      });
      // static_assert? if there are a lots of channels this won't run on
      // registers.
    }
    return out;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  T const* _h = nullptr; // kernel
};

static constexpr uint fir_cache_line_bytes = 128;

template <class T>
using fir_std_vector
  = std::vector<T, overaligned_allocator<T, fir_cache_line_bytes>>;

constexpr bool fir_verify_alignment (uint v)
{
  return is_pow2 (v) && v <= fir_cache_line_bytes;
}

} // namespace detail
//------------------------------------------------------------------------------
// A non-FFT convolution/FIR filter. The kernel and delay line memory are not
// owned by the class. This is so because the most optimal memory layout is not
// known at this abstraction level.
//
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
// overview of the FIR classes.
//
// There is "fir_decimator", "fir_interpolator" and "lth_band_fir_decimator".
//
// Each of them has an "alignment" template parameter that will make:
//
// - The kernel size in bytes to be rounded (and zero padded) up to the data
//   alignment.
// - The delay line size(s) in bytes to be rounded up to the data alignment.
//
// The intent is to avoid a possible compiler-generated slow path to deal with
// scalars by filling with zeros. It might matter on some machines. On a Ryzen
// 5800x it is a slight but measureable pessimization.
//
// Then there is an Lth variant of the Decimator and a Lth band optimization
// flag on the interpolator.
//
// Non strict measurements on a Ryzen 5800x have shown the optimization to be
// slightly working at 2x-4x (which was surprising, I was expecting at least
// 25%) and slightly detrimental equal and over 8x.

template <class T, uint channels = 1, uint alignment = alignof (T)>
class fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel, uint ratio)
  {
    constexpr uint overalign = alignment / sizeof (T);

    uint ksize_ratio     = round_ceil<uint> (kernel.size(), ratio);
    uint ksize_overalign = round_ceil<uint> (ksize_ratio, overalign);

    uint ksize   = ksize_overalign;
    uint delsize = 2 * ksize * n_channels;

    _mem.clear();
    _mem.resize (ksize + delsize); // 0 filled
    _ratio = ratio;
    // On the decimator the zero coefficients only added for the filter to be
    // polyphase decomposable have to be at the front.
    memcpy (
      _mem.data() + (ksize_ratio - kernel.size()),
      kernel.data(),
      sizeof kernel[0] * kernel.size());

    _filter.reset (make_crange (&_mem[0], ksize));
    _delay.reset (make_crange (&_mem[ksize], delsize), ksize);
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
  detail::convolution_block<T, n_channels, alignment> _filter;
  detail::convolution_delay_line<T, n_channels>       _delay;
  detail::fir_std_vector<T>                           _mem;
  uint                                                _ratio;
};
//------------------------------------------------------------------------------
// Polyphase FIR decimator optimized for L-th band decimation (coefficients with
// a regular pattern of zeros).
template <class T, uint channels = 1, uint alignment = alignof (T)>
class lth_band_fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<const T> kernel, uint ratio)
  {
    assert (lth_band::verify_coeffs (kernel, ratio));

    constexpr uint overalign = alignment / sizeof (T);

    uint ksize_ratio     = round_ceil<uint> (kernel.size(), ratio);
    uint ksize_overalign = round_ceil<uint> (ksize_ratio, overalign * ratio);

    uint ksize      = ksize_overalign;
    uint subk_size  = ksize / ratio;
    uint kernel_mem = ksize - subk_size; // remove zeroes
    // magic "2" reminder: delay line buffer has to be double the required size

    // position on the delay line handling the zero coefficients.
    _center_sample_pos = ksize_ratio / (ratio * 2);
    // The delay line with the zeroed coefficients only has a coefficient at the
    // center, so the delay line can be shortened on half + 1.
    uint zero_dl_size = round_ceil<uint> (
      (_center_sample_pos + 1) * (2 * n_channels), overalign);
    uint zero_dl_elems = zero_dl_size / (2 * n_channels);

    uint non_zero_dl_elems = ksize - subk_size;
    uint non_zero_dl_size  = non_zero_dl_elems * (2 * n_channels);

    _mem.clear();
    _mem.resize (kernel_mem + zero_dl_size + non_zero_dl_size); // 0 filled

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

    T* dst = &_mem[0];
    _filter.reset (make_crange (dst, kernel_mem));
    dst += kernel_mem;

    _delays[0].reset (make_crange (dst, non_zero_dl_size), non_zero_dl_elems);
    dst += non_zero_dl_size;

    _delays[1].reset (make_crange (dst, zero_dl_size), zero_dl_elems);
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
    // skipping positions with padded zero coefficients
    uint n_trailing_zero_coeffs = _ratio - 1;
    alignas (alignment) std::array<T, n_channels> ret
      = _filter.tick (_delays[0], n_trailing_zero_coeffs);

    // adding the only non zero coeff at the center of the remaining branch.
    auto z_ptrs = _delays[1].samples();
    for (uint c = 0; c < n_channels; ++c) {
      ret[c] += z_ptrs[c][_center_sample_pos] * _center_coeff;
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
  detail::fir_std_vector<T>                                    _mem;
  uint _center_sample_pos;
  T    _center_coeff;
  uint _ratio;
};
//------------------------------------------------------------------------------
// Polyphase FIR interpolator with selectable L-th band optimization (in case
// the coefficients are suitable for Lth band optimization (skipping zeroes)).
template <class T, uint channels = 1, uint alignment = alignof (T)>
class fir_interpolator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  // If the filter is a Lth-band filter it has "1/ratio" of the coefficients as
  // zeros, so they can be skipped. A polyphase FIR interpolator is structured
  // in a way that the skipping the zeroed branch is trivial, so it doesn't have
  // a separate implementation as the decimator does.
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

    constexpr uint overalign = alignment / sizeof (T);

    uint ksize_ratio     = round_ceil<uint> (kernel.size(), ratio);
    uint ksize_overalign = round_ceil<uint> (ksize_ratio, overalign * ratio);

    uint ksize      = ksize_overalign;
    uint subk_size  = ksize / ratio;
    uint kernel_mem = ksize;
    kernel_mem -= (_is_lth_band) ? subk_size : 0;
    // magic "2" reminder: delay line buffer has to be double the required size
    uint delay_lines_mem = 2 * subk_size * n_channels;

    // the center coefficient position before padding the delay line.
    _center_sample_pos = ksize_ratio / (2 * ratio);

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
      auto z_ptrs = _delay.samples();
      for (uint c = 0; c < channels; ++c) {
        out[c][0] = z_ptrs[c][_center_sample_pos];
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
  detail::fir_std_vector<T>                             _mem;
  uint                                                  _center_sample_pos;
  uint                                                  _ratio;
  bool                                                  _is_lth_band;
};
//------------------------------------------------------------------------------

} // namespace artv
