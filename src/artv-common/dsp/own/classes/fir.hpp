#pragma once

// Possible (known) optimizations left:
//
// -Handrolled SSE version. Some performance might be gained if some dev time is
//    spent optimizing this. Worth?

#include <algorithm>
#include <cassert>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

#include "artv-common/dsp/own/classes/windowed_sync.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/math.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

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
// TODO: maybe this is generic enough to be moved to delay_line.hpp?
template <class T, uint channels = 1>
class convolution_delay_line {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (xspan<T> delay_line_mem, uint size)
  {
    assert (delay_line_mem.size() >= (size * n_channels * 2));
    _head = 0;
    _size = size;
    _z    = delay_line_mem.data();
    memset (_z, 0, sizeof _z[0] * _size * 2 * n_channels);
  }
  //----------------------------------------------------------------------------
  // interleaved input
  void push (xspan<T const> v)
  {
    assert (v.size() >= n_channels);

    _head = _head == 0 ? _size : _head;
    --_head;

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head1 = _z + _head + (chnl * _size * 2);
      T* chnl_head2 = chnl_head1 + _size;
      *chnl_head1 = *chnl_head2 = v[chnl];
    });
  }

  void push (xspan<T const> v, uint count)
  {
    assert (v.size() >= (n_channels * count));
    for (uint i = 0; i < count; ++i) {
      push (xspan {&v[i * n_channels], n_channels});
    }
  }
  //---------------------------------------------------------------------------
  // instead of "push" "push_leading" and "copy_trailing" can be used. This is
  // for probable cache-friendliness, so the tail is written near the last
  // touched memory location.
  void push_leading (xspan<T const> v)
  {
    assert (v.size() >= n_channels);

    _head = _head == 0 ? _size : _head;
    --_head;

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head1 = _z + _head + (chnl * _size * 2);
      *chnl_head1   = v[chnl];
    });
  }

  void prepare_trailing()
  {
    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head1 = _z + _head + (chnl * _size * 2);
      T* chnl_head2 = chnl_head1 + _size;
      *chnl_head2   = *chnl_head1;
    });
  }

  void prepare_trailing (xspan<T const> v)
  {
    assert (v.size() >= n_channels);

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head2 = _z + _head + (chnl * _size * 2) + _size;
      *chnl_head2   = v[chnl];
    });
  }
  //----------------------------------------------------------------------------
  // sparse input (all the leading/trailing functions TBD if required)
  void push (std::array<xspan<T const>, n_channels> v, uint count = 1)
  {
    for (auto& r : v) {
      assert (r.size() >= count);
    }
    for (uint i = 0; i < count; ++i) {
      _head = _head == 0 ? _size : _head;
      --_head;

      mp_foreach_idx<n_channels> ([&] (auto chnl) {
        T* chnl_head1 = _z + _head + (chnl * _size * 2);
        T* chnl_head2 = chnl_head1 + _size;
        *chnl_head1 = *chnl_head2 = v[chnl][i];
      });
    }
  }
  //----------------------------------------------------------------------------
  std::array<T const * artv_restrict, n_channels> samples() const
  {
    std::array<T const * artv_restrict, n_channels> ret;

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
  void reset (xspan<T const> kernel_mem) { _h = kernel_mem.data(); }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (
    convolution_delay_line<T, n_channels> const& dl)
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
  void reset (xspan<T const> kernel_mem, xspan<T> delay_line_mem)
  {
    _impl.reset (kernel_mem);
    _dl.reset (delay_line_mem, kernel_mem.size());
  }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (std::array<T, n_channels> in)
  {
    _dl.push (in);
    return _impl.tick (_dl);
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
  static bool verify_coeffs (xspan<T> const kernel, uint ratio)
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
//------------------------------------------------------------------------------
// overview of the FIR classes for resampling power of 2 ratios .
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
//

template <class T, uint channels = 1, uint alignment = alignof (T)>
class fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (xspan<T const> kernel, uint ratio)
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

    _filter.reset (xspan {&_mem[0], ksize});
    _delay.reset (xspan {&_mem[ksize], delsize}, ksize);
  }
  //----------------------------------------------------------------------------
  // inputs spread.
  std::array<T, n_channels> tick (std::array<xspan<T const>, n_channels> in)
  {
    for (auto& r : in) {
      assert (r.size() >= _ratio);
    }
    _delay.push (in, _ratio);
    return _filter.tick (_delay);
  }
  //----------------------------------------------------------------------------
  // inputs interleaved
  std::array<T, n_channels> tick (xspan<T const> in)
  {
    assert (in.size() >= _ratio * n_channels);
    _delay.push (in, _ratio);
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
//------------------------------------------------------------------------------
// Polyphase FIR decimator optimized for L-th band decimation (coefficients with
// a regular pattern of zeros).
template <class T, uint channels = 1, uint alignment = alignof (T)>
class lth_band_fir_decimator {
public:
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (xspan<T const> kernel, uint ratio)
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
    _filter.reset (xspan {dst, kernel_mem});
    dst += kernel_mem;

    _delays[0].reset (xspan {dst, non_zero_dl_size}, non_zero_dl_elems);
    dst += non_zero_dl_size;

    _delays[1].reset (xspan {dst, zero_dl_size}, zero_dl_elems);
  }
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick (std::array<xspan<T const>, n_channels> in)
  {
    _delays[1].push (in);
    for (auto& r : in) {
      r.cut_head (1);
    }
    _delays[0].push (in, _ratio - 1);
    return tick_after_input();
  }
  //----------------------------------------------------------------------------
  // interleaved version
  std::array<T, n_channels> tick (xspan<T const> in)
  {
    _delays[1].push (in.get_head (1));
    _delays[0].push (in.advanced (n_channels), _ratio - 1);
    return tick_after_input();
  }
  //----------------------------------------------------------------------------
  uint ratio() const { return _ratio; }
  //----------------------------------------------------------------------------
  uint order() const { return (_filter.order() / (_ratio - 1)) * _ratio; }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  std::array<T, n_channels> tick_after_input()
  {
    // skipping positions with padded zero coefficients
    alignas (alignment) std::array<T, n_channels> ret;

    uint n_trailing_zero_coeffs = _ratio - 1;
    ret = _filter.tick (_delays[0], n_trailing_zero_coeffs);

    // adding the only non zero coeff at the center of the remaining branch.
    auto z_ptrs = _delays[1].samples();
    for (uint c = 0; c < n_channels; ++c) {
      ret[c] += z_ptrs[c][_center_sample_pos] * _center_coeff;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  detail::convolution_block<T, n_channels>                     _filter;
  std::array<detail::convolution_delay_line<T, n_channels>, 2> _delays;
  detail::fir_std_vector<T>                                    _mem;
  uint _center_sample_pos;
  T    _center_coeff;
  uint _ratio;
};
//------------------------------------------------------------------------------
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
    xspan<T const> kernel,
    uint           ratio,
    bool           enable_lth_band_optimization = false)
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
    _delay.reset (xspan {&_mem[kernel_mem], delay_lines_mem}, subk_size);

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
      _filters[subk - _is_lth_band].reset (xspan {&_mem[offset], subk_size});
    }
  }
  //----------------------------------------------------------------------------
  void tick (std::array<xspan<T>, n_channels> out, xspan<T const> in)
  {
    for (auto& r : out) {
      assert (r.size() >= _ratio);
    }
    assert (in.size() >= n_channels);

    _delay.push (in);
    uint is_lth_band = _is_lth_band;
    for (uint i = is_lth_band; i < _ratio; ++i) {
      auto interleaved = _filters[i - is_lth_band].tick (_delay);
      for (uint c = 0; c < n_channels; ++c) {
        out[c][i] = interleaved[c];
      }
    }
    if (unlikely (is_lth_band)) {
      // Handle the central coefficient, as this class is normalizing the
      // kernel to compensate for the energy spread on the upsampling imaging,
      // the center coefficient on a L-th band linear phase lowpass is always 1,
      // so it is a direct copy of the sample, as multiplying by one has no
      // effect.
      auto z_ptrs = _delay.samples();
      for (uint c = 0; c < n_channels; ++c) {
        out[c][0] = z_ptrs[c][_center_sample_pos];
      }
    }
  }
  //----------------------------------------------------------------------------
  // interleaved out
  void tick (xspan<T> out, xspan<T const> in)
  {
    assert (out.size() >= (_ratio * n_channels));
    assert (in.size() >= n_channels);

    _delay.push (in);
    uint is_lth_band = _is_lth_band;
    for (uint i = is_lth_band; i < _ratio; ++i) {
      auto interleaved = _filters[i - is_lth_band].tick (_delay);
      memcpy (&out[n_channels * i], interleaved.data(), sizeof interleaved);
    }
    if (unlikely (is_lth_band)) {
      // Handle the central coefficient, as this class is normalizing the
      // kernel to compensate for the energy spread on the upsampling imaging,
      // the center coefficient on a L-th band linear phase lowpass is always 1,
      // so it is a direct copy of the sample, as multiplying by one has no
      // effect.
      auto z_ptrs = _delay.samples();
      for (uint c = 0; c < channels; ++c) {
        out[c] = z_ptrs[c][_center_sample_pos];
      }
    }
  }
  // interleaved version
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
namespace resampling {
//------------------------------------------------------------------------------
static constexpr bool is_downsampler (fraction<uint> ratio)
{
  return ratio.num < ratio.den;
}
//------------------------------------------------------------------------------
static constexpr bool is_upsampler (fraction<uint> ratio)
{
  return ratio.num > ratio.den;
}
//------------------------------------------------------------------------------
static constexpr bool is_no_resampler (fraction<uint> ratio)
{
  return ratio.num == ratio.den;
}
//------------------------------------------------------------------------------
// max samples in one iteration.
static constexpr uint max_samples (fraction<uint> ratio)
{
  return div_ceil (ratio.num, ratio.den);
}
//------------------------------------------------------------------------------
// gets how many outputs for a given number of inputs gets generated by a
// resampler.
static constexpr uint get_n_new_outs_for_n_ticks (
  fraction<uint> ratio,
  fraction<uint> pos,
  uint           n_ticks)
{
  uint in_rate    = ratio.num; // rate in units of time (less = faster)
  uint out_rate   = ratio.den; // rate in units of time (less = faster)
  uint n_spls_in  = pos.den;
  uint n_spls_out = pos.num;

  uint new_t          = (n_spls_in + n_ticks) * in_rate;
  uint n_spls_out_new = new_t / ratio.den;
  uint n_spls         = n_spls_out_new - n_spls_out;
  return n_spls;
}
//------------------------------------------------------------------------------
// number of ticks/input samples to generate the nearest sample amount of output
// samples (n_spls).
//
// More samples than "n_spls" might be generated by the returned number of
// ticks/input samples. It will never be exceeded by more than
// "max_n_samples(ratio)".
//------------------------------------------------------------------------------
static constexpr uint get_n_ticks_for_n_new_spls_ceil (
  fraction<uint> ratio, // resampling ratio
  fraction<uint> pos, // position on the cycle (normalized to the ratio)
  uint           n_spls) // desired number of samples
{
  uint in_rate    = ratio.num; // rate in units of time (less = faster)
  uint out_rate   = ratio.den; // rate in units of time (less = faster)
  uint n_spls_in  = pos.den;
  uint n_spls_out = pos.num;

  uint t_desired = (n_spls_out + n_spls) * out_rate;
  // round to the ceiling, as e.g. for upsamplers a regular division might not
  // advance the input position if e.g. "n_spls" is "1"
  uint n_in_desired = div_ceil (t_desired, in_rate);
  return n_in_desired - n_spls_in;
}
//------------------------------------------------------------------------------
// number of ticks/input samples to generate the nearest sample amount of output
// samples (n_spls).
//
// "n_spls" will never be surpassed, so sometimes less samples than requested
// might be generated by the returned number of ticks. Corollary: this function
// might return 0 samples if n_spls is small compared to the ratio.
//------------------------------------------------------------------------------
static constexpr uint get_n_ticks_for_n_new_spls_floor (
  fraction<uint> ratio, // resampling ratio
  fraction<uint> pos, // position on the cycle (normalized to the ratio)
  uint           n_spls) // desired number of samples
{
  uint in_rate    = ratio.num; // rate in units of time (less = faster)
  uint out_rate   = ratio.den; // rate in units of time (less = faster)
  uint n_spls_in  = pos.den;
  uint n_spls_out = pos.num;

  uint n_ticks    = get_n_ticks_for_n_new_spls_ceil (ratio, pos, n_spls);
  uint n_new_spls = get_n_new_outs_for_n_ticks (ratio, pos, n_ticks);
  uint ret        = n_ticks - (uint) (n_new_spls > n_spls);
  return ret;
}
} // namespace resampling
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// A fractional resampler for integer samplerates. Useful for e.g. 44100 to
// 48000 conversions.
//
// The usual approach in many libraries is either to calculate the sinc and
// window at runtime, which restricts the choice of windows to those that are
// easier to compute (and perform so-so) or to have a high(ish) number of tables
// with precomputed windowed sincs using a proper window and to interpolate
// between the closest.
//
// This resampler calculates the periodicity of both sample rates and stores
// a single table for each sample that has to be output on the full cycle.
//
// This approach has both advantages and caveats vs interpolating tables.
// Caveats:
//
// - It doesn't allow to support any rate with a given memory requirement. Some
//   rates might require a lot of memory.
//
// Advantages:
//
// - The memory access pattern is always ascending. Interpolating might jump
//   around between memory locations depending on the ratio, which might be less
//   friendly to the cache.
//
// - There is no interpolation to do.
//
// - For some serious resampler 128-256 tables aren't unheard of, for the main
//   use case of this 44100/48000, 88200/96000 conversions the number of tables
//   required is 147/160 depending on the direction (gcd is 300).
//
//   If this is used to set the internal frequency at wich a DSP process
//   operates then some  favorable numbers for the internal samplerate can be
//   chosen, as e.g. 54000, which has very low memory requirements for multiples
//   of 48KHz (gcd is 6000) and low for multiples of 44100 (gcd is 900).
//
// Along with "resampling::get_n_ticks_for_n_new_spls_floor/ceil" it is possible
// to precompute how many input samples should be fed to get N output samples.
// E.g.:
//
// ticks = resampling::get_n_ticks_for_n_new_spls_floor(
//    r.ratio(), r.corrected_pos(), x);

template <class T, uint channels = 1, uint alignment = alignof (T)>
class fractional_resampler {
public:
  //----------------------------------------------------------------------------
  using value_type                 = T;
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (
    uint  tgt_srate,
    uint  src_srate,
    uint  taps,
    float cutoff_hz,
    float kaiser_beta,
    bool  minphase)
  {
    uint gcd_v    = gcd (tgt_srate, src_srate);
    uint rate_tgt = tgt_srate / gcd_v;
    uint rate_src = src_srate / gcd_v;

    _ratio.den = rate_src;
    _ratio.num = rate_tgt;
    // using u8 for "_max_samples"
    assert (_ratio.den < 256 && _ratio.num < 256);
    _max_samples = resampling::max_samples (_ratio);
    _kernel_size = taps;
    memory_reset (rate_tgt);

    assert ((cutoff_hz * 2.f) <= (float) tgt_srate);
    auto  ratio = (float) rate_tgt / (float) rate_src;
    float fc;
    if (ratio <= 1.f) {
      fc = cutoff_hz / src_srate;
      fc /= rate_tgt;
    }
    else {
      fc = cutoff_hz / tgt_srate;
      fc /= rate_src;
    }

    std::vector<value_type> bigkernel {};
    bigkernel.resize (taps * rate_tgt);

    kaiser_lp_kernel_2<value_type> (
      bigkernel, fc, kaiser_beta, rate_tgt, minphase);

    // t1 = time of output sample
    // t2 = time of previous output sample
    // n1 = number of input samples at t1
    // n2 = number of input samples at t2

    auto t1 = rate_src * rate_tgt;
    auto t2 = t1 - rate_tgt;
    auto n2 = t2 / rate_src;

    uint subk = 0;
    for (uint i = 0; i < rate_src; ++i) {
      auto n1     = t1 / rate_src;
      auto n_spls = n1 - n2;
      assert (n_spls <= std::numeric_limits<u8>::max());

      for (uint spl = 1; spl < (n_spls + 1); ++spl) {
        auto t_in_spl = (n2 + spl) * rate_src;
        auto offset   = t1 - t_in_spl;
        offset        = rate_tgt - offset - 1;

        auto subkernel = get_subkernel (subk);
        for (uint j = 0; j < _kernel_size; ++j) {
          subkernel[j] = bigkernel[offset];
          offset += rate_tgt;
        }
        fir_kernel_normalize (subkernel);
        ++subk;
      }
      n2 = n1;
      t1 += rate_tgt;
    }
    _n_in  = 0;
    _n_out = 0;
  }
  //----------------------------------------------------------------------------
  uint next_tick_n_spls_out() const
  {
    return resampling::get_n_new_outs_for_n_ticks (ratio(), corrected_pos(), 1);
  }
  //----------------------------------------------------------------------------
  // ratio as a fraction, the first element is the numerator
  fraction<uint> ratio() const { return _ratio; };
  //----------------------------------------------------------------------------
  // position in the resampling cycle, normalized to the ratio (read comment for
  // function below).
  fraction<uint> pos() const { return {_n_out, _n_in}; };
  //----------------------------------------------------------------------------
  // On this class the output cycle the cycle is shifted forwards by one, so the
  // maximum amount of samples possible for a ratio is delivered on the first
  // sample. get_n_ticks_for_n_new_spls_ceil/floor need to use this corrected
  // function instead of "pos()".
  //
  // E.g. normally on a 5/4 resampler the cycle of output samples is: 1+1+1+2,
  // this class returns 2+1+1+1. For a 1/4 resampler the cycle is 0+0+0+1, this
  // class returns 1+0+0+0
  //
  // This function doesn't return a normalized fraction against the ratio
  // because of the shifting.
  //----------------------------------------------------------------------------
  fraction<uint> corrected_pos() const
  {
    auto ret = pos();
    ret.den += _ratio.den - 1;
    ret.num += _ratio.num - _max_samples;
    return ret;
  }
  //----------------------------------------------------------------------------
  // outputs and takes interleaved samples, "out" has to contain enough space
  // to contain "next_tick_n_spls_out" multichannel elements.
  uint tick (
    xspan<value_type>       out,
    xspan<value_type const> ins,
    uint                    block_size = 1)
  {
    auto in     = ins.get_head (block_size * n_channels);
    uint n_spls = 0;

    while (in.size()) {
      auto n_spls_now = next_tick_n_spls_out();
      auto end        = _n_out + n_spls_now;
      n_spls += n_spls_now;

      // not using "push_leading/prepare_trailing" instead of "push" so "in"
      // and "out" can alias.
      _delay.push (in.cut_head (n_channels));

      while (_n_out < end) {
        detail::convolution_block<value_type, n_channels, alignment> dotprod;
        dotprod.reset (get_subkernel (_n_out));
        auto spls = dotprod.tick (_delay);
        auto dst  = out.cut_head (n_channels);
        static_assert (sizeof spls == sizeof dst[0] * n_channels);
        memcpy (dst.data(), spls.data(), sizeof spls);
        ++_n_out;
      }

      ++_n_in;
      if (_n_in >= _ratio.den) {
        assert (_n_out == _ratio.num);
        _n_in = _n_out = 0;
      }
    }

    return n_spls;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  xspan<value_type> get_subkernel (uint pos)
  {
    return {(T*) &_kernels[_kernel_bytes * pos], _kernel_size};
  }
  //----------------------------------------------------------------------------
  void memory_reset (uint n_kernel_tables)
  {
    _kernel_bytes = round_ceil<uint> (_kernel_size * sizeof (T), alignment);

    uint mem_align      = std::max<uint> (alignment, sizeof (T));
    uint del_elems      = n_channels * _kernel_size * 2;
    uint del_size_bytes = round_ceil<uint> (del_elems * sizeof (T), mem_align);
    uint k_size_bytes   = _kernel_bytes * n_kernel_tables;

    _mem.resize (del_size_bytes + k_size_bytes);
    _delay.reset (xspan {(T*) &_mem[0], del_elems}, _kernel_size);
    _kernels = &_mem[del_size_bytes];
  }
  //----------------------------------------------------------------------------
  void write_subkernel (
    xspan<T> bigkernel,
    uint     offset,
    uint     pos,
    uint     rate_tgt)
  {
    auto kern = (T*) _kernels[_kernel_bytes * pos];
    for (uint i = 0; i < _kernel_size; ++i) {
      kern[i] = bigkernel[offset];
      offset += rate_tgt;
    }
    fir_kernel_normalize ({kern, _kernel_size});
  }
  //----------------------------------------------------------------------------
  detail::fir_std_vector<u8>                    _mem;
  detail::convolution_delay_line<T, n_channels> _delay;
  u8*                                           _kernels;
  uint                                          _kernel_size;
  uint                                          _kernel_bytes;
  uint                                          _n_in;
  uint                                          _n_out;
  fraction<uint>                                _ratio;
  u8                                            _max_samples;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// A polyphase decimator or interpolator taking care of power of two factors
// combined with a fractional resampler. Be sure to measure performance, as
// the latency is increased.
//
// This also provides is blockwise processing of both resampling stages.
//
// Along with "resampling::get_n_ticks_for_n_new_spls_floor/ceil" it is possible
// to precompute how many input samples should be fed to get N output samples.
// E.g.:
//
// ticks = resampling::get_n_ticks_for_n_new_spls_floor(
//    r.ratio(), r.corrected_pos(), x);
//------------------------------------------------------------------------------
template <class T, uint channels = 1, uint alignment = alignof (T)>
class resampler {
public:
  //----------------------------------------------------------------------------
  using value_type                           = T;
  static constexpr uint n_channels           = channels;
  static constexpr uint max_samples          = 128;
  static constexpr uint max_stack_size_bytes = 4 * 1024;
  // The magic "2" accounts for fractional resampler ratios of 1.9999999
  static constexpr uint max_block_size
    = max_stack_size_bytes / (2 * n_channels * sizeof (value_type));
  //----------------------------------------------------------------------------
  void reset (
    uint  tgt_srate,
    uint  src_srate,
    uint  taps_int,
    uint  taps_frac,
    float cutoff_hz,
    float kaiser_beta,
    bool  minphase)
  {
    bool is_downsampler = tgt_srate <= src_srate;
    _integer            = std::monostate {};
    _fractional.reset();
    _int_ratio.reset (1, 1);

    uint small     = std::min (tgt_srate, src_srate);
    uint big       = std::max (tgt_srate, src_srate);
    uint int_ratio = last_bit_set (big / small);
    // 0 is no bit set, 1 is a ratio of 1, 2 a ratio of 2, 3 a ratio of 4, etc
    int_ratio = (int_ratio >= 2) ? 1 << (int_ratio - 1) : 1;

    if (is_downsampler) {
      uint frac_srate = src_srate;

      if (int_ratio != 1) {
        frac_srate /= int_ratio;
        std::vector<value_type> tmp_kernel;
        tmp_kernel.resize (taps_int * int_ratio);

        auto fc = cutoff_hz / (float) src_srate;

        // TODO: fractional delay compensation?
        kaiser_lp_kernel_2<value_type> (
          tmp_kernel, fc, kaiser_beta, 0, minphase);

        _integer        = decimator_type {};
        auto& decim_rsc = std::get<decimator_type> (_integer);

        decim_rsc.decimator.reset (tmp_kernel, int_ratio);
        decim_rsc.input.reserve (int_ratio * n_channels);
        // Insert zeros: synchronize so the first sample input generates
        // output.
        decim_rsc.input.resize ((int_ratio - 1) * n_channels);

        _int_ratio.den = int_ratio;
      }

      if (tgt_srate != frac_srate) {
        _fractional.emplace();
        _fractional->reset (
          tgt_srate, frac_srate, taps_frac, cutoff_hz, kaiser_beta, minphase);
      }
    }
    else {
      uint frac_srate = tgt_srate;

      if (int_ratio != 1) {
        frac_srate /= int_ratio;
        std::vector<value_type> tmp_kernel;
        tmp_kernel.resize (taps_int * int_ratio);

        auto fc = cutoff_hz / (float) tgt_srate;
        // TODO: fractional delay compensation?
        kaiser_lp_kernel_2<value_type> (
          tmp_kernel, fc, kaiser_beta, 0, minphase);

        _integer     = interpolator_type {};
        auto& interp = std::get<interpolator_type> (_integer);
        interp.reset (tmp_kernel, int_ratio);

        _int_ratio.num = int_ratio;
      }

      if (src_srate != frac_srate) {
        _fractional.emplace();
        _fractional->reset (
          frac_srate, src_srate, taps_frac, cutoff_hz, kaiser_beta, minphase);
      }
    }
    _max_samples = resampling::max_samples (ratio());
  }
  //----------------------------------------------------------------------------
  // ratio as a fraction, the first element is the target rate, the second the
  // source.
  fraction<uint> ratio() const
  {
    if (_fractional) {
      auto ratio = _fractional->ratio();
      ratio.num *= _int_ratio.num;
      ratio.den *= _int_ratio.den;
      return ratio;
    }
    else {
      return _int_ratio;
    }
  }
  //----------------------------------------------------------------------------
  // position in the resampling cycle, normalized to the ratio (read comment for
  // function below).
  fraction<uint> pos() const
  {
    uint int_pos = 0;
    if (std::holds_alternative<decimator_type> (_integer)) {
      int_pos = std::get<decimator_type> (_integer).input.size();
    }
    if (_fractional) {
      auto pos = _fractional->pos();
      pos.num *= _int_ratio.num;
      pos.den *= _int_ratio.den;
      return pos;
    }
    else {
      // the integer ratio interpolators have a cycle of 1 so we don't have to
      // account them on ther denominator.
      // the integer ratio decimators always output 1 sample at most in the
      // whole cycle, so they don't need to be accounted on the numerator;
      return {0, int_pos};
    }
  };
  //----------------------------------------------------------------------------
  // On this class the output cycle the cycle is shifted forwards by one, so the
  // maximum amount of samples possible for a ratio is delivered on the first
  // sample. get_n_ticks_for_n_new_spls_ceil/floor need to use this corrected
  // function instead of "pos()".
  //
  // E.g. normally on a 5/4 resampler the cycle of output samples is: 1+1+1+2,
  // this class returns 2+1+1+1. For a 1/4 resampler the cycle is 0+0+0+1, this
  // class returns 1+0+0+0
  //
  // This function doesn't return a normalized fraction against the ratio
  // because of the shifting.
  //----------------------------------------------------------------------------
  fraction<uint> corrected_pos() const
  {
    uint int_pos = 0;
    if (std::holds_alternative<decimator_type> (_integer)) {
      int_pos = std::get<decimator_type> (_integer).input.size();
      int_pos /= n_channels;
    }
    if (_fractional) {
      // the fractional resample already accounts for the shift
      auto pos = _fractional->corrected_pos();
      pos.num *= _int_ratio.num;
      pos.den *= _int_ratio.den;
      return pos;
    }
    else {
      return {0, int_pos};
    }
  };
  //----------------------------------------------------------------------------
  uint tick_upsampler (
    xspan<value_type>       out,
    xspan<value_type const> in,
    uint                    block_size)
  {
    assert (_max_samples != 1);
    assert (in.size() >= (n_channels * block_size));
    // Using VLA for locality, this limits the "block_size" on this call. .
    // Normally "block_size" should be a fixed internal size of around 32/64.
    assert (block_size < max_block_size);
    // Magic "2", acounting for ratios of 1.99999999...
    uint              frac_mem_elems = 2 * block_size * n_channels;
    value_type        frac_mem[frac_mem_elems];
    xspan<value_type> fracbf {&frac_mem[0], frac_mem_elems};

    // the no-src change case is not handled.
    uint n_spls = 0;

    if (_fractional) {
      auto in_frac  = xspan {in.data(), block_size * n_channels};
      auto out_frac = fracbf;
      n_spls        = _fractional->tick (out_frac, in_frac, block_size);
    }

    if (std::holds_alternative<interpolator_type> (_integer)) {
      auto& interpolator = std::get<interpolator_type> (_integer);
      auto  ratio        = interpolator.ratio();

      xspan<value_type const> in_int;
      if (n_spls != 0) {
        in_int = xspan {fracbf.data(), n_spls * n_channels};
      }
      else {
        in_int = xspan {in.data(), block_size * n_channels};
        n_spls = block_size;
      }
      while (in_int.size()) {
        interpolator.tick (
          out.cut_head (n_channels * ratio), in_int.cut_head (n_channels));
      }
      n_spls *= ratio;
    }
    else if (n_spls) {
      assert (out.size() >= n_spls * n_channels);
      memcpy (out.data(), fracbf.data(), n_spls * n_channels * sizeof out[0]);
    }
    return n_spls;
  }
  //----------------------------------------------------------------------------
  uint tick_downsampler (
    xspan<value_type>       out,
    xspan<value_type const> in,
    uint                    block_size)
  {
    assert (_max_samples == 1);
    assert (block_size < max_block_size);
    assert (in.size() >= (block_size * n_channels));

    // Magic "2" because if there is an integer resampler, its more
    // unfavorable ratio will be 2. This buffer is smaller than the one used
    // when upsampling.
    uint       int_mem_elems = div_ceil<uint> (block_size, 2) * n_channels;
    value_type int_mem[int_mem_elems];
    xspan<value_type> intbf {&int_mem[0], int_mem_elems};

    uint n_spls    = 0;
    bool processed = false;

    if (std::holds_alternative<decimator_type> (_integer)) {
      processed          = true;
      auto& decimator    = std::get<decimator_type> (_integer).decimator;
      auto& decimator_in = std::get<decimator_type> (_integer).input;
      auto  ready_n_spls = decimator.ratio() * n_channels;

      auto in_int  = xspan {in.data(), block_size * n_channels};
      auto out_int = intbf;

      while (in_int.size()) {
        uint n_insert
          = std::min (ready_n_spls - decimator_in.size(), in_int.size());
        decimator_in.insert (
          decimator_in.end(), in_int.begin(), in_int.begin() + n_insert);
        in_int.cut_head (n_insert);

        if (decimator_in.size() == ready_n_spls) {
          auto spls_arr = decimator.tick (decimator_in);
          decimator_in.clear();
          memcpy (out_int.data(), spls_arr.data(), sizeof spls_arr);
          out_int.cut_head (n_channels);
          n_spls += 1;
        }
      }
    }
    if (_fractional) {
      if (n_spls || !processed) {
        xspan<value_type const> in_frac;
        uint                    bsz;

        if (processed) {
          bsz     = n_spls;
          in_frac = intbf.get_head (bsz * n_channels);
        }
        else {
          bsz     = block_size;
          in_frac = xspan {in.data(), bsz * n_channels};
        }
        processed = true;
        n_spls    = _fractional->tick (out, in_frac, bsz);
      }
    }
    else if (processed) {
      assert (out.size() >= n_spls * n_channels);
      memcpy (out.data(), intbf.data(), n_spls * n_channels * sizeof out[0]);
    }
    else if (unlikely (!processed)) { // this is an "else" but with the hint.
      // no src change.
      assert (out.size() >= block_size * n_channels);
      n_spls = block_size;
      memcpy (out.data(), in.data(), n_spls * n_channels * sizeof out[0]);
    }
    return n_spls;
  }
  //----------------------------------------------------------------------------
  uint tick (xspan<value_type> out, xspan<value_type const> in, uint block_size)
  {
    if (_max_samples == 1) {
      return tick_downsampler (out, in, block_size);
    }
    else {
      return tick_upsampler (out, in, block_size);
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  using interpolator_type = fir_interpolator<T, channels, alignment>;

  struct decimator_type {
    fir_decimator<T, channels, alignment> decimator;
    std::vector<value_type>               input;
  };

  std::variant<std::monostate, decimator_type, interpolator_type> _integer;
  std::optional<fractional_resampler<T, channels, alignment>>     _fractional;
  fraction<uint>                                                  _int_ratio;
  u8                                                              _max_samples;
};

//------------------------------------------------------------------------------

} // namespace artv
