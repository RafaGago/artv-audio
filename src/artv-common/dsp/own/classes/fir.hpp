#pragma once

// Possible (known) optimizations left:
//
// -Handrolled SSE version. Some performance might be gained if some dev time is
//    spent optimizing this. Worth?

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
// TODO: maybe this is generic enough to be moved to delay_line.hpp?
template <class T, uint channels = 1>
class convolution_delay_line {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = channels;
  //----------------------------------------------------------------------------
  void reset (crange<T> delay_line_mem, uint size)
  {
    assert (delay_line_mem.size() >= (size * n_channels * 2));
    _head = 0;
    _size = size;
    _z    = delay_line_mem.data();
    memset (_z, 0, sizeof _z[0] * _size * 2 * n_channels);
  }
  //----------------------------------------------------------------------------
  // interleaved input
  void push (crange<const T> v)
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

  void push (crange<const T> v, uint count)
  {
    assert (v.size() >= (n_channels * count));
    for (uint i = 0; i < count; ++i) {
      push (make_crange (&v[i * n_channels], n_channels));
    }
  }
  //---------------------------------------------------------------------------
  // instead of "push" "push_leading" and "copy_trailing" can be used. This is
  // for probable cache-friendliness, so the tail is written near the last
  // touched memory location.
  void push_leading (crange<const T> v)
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

  void prepare_trailing (crange<const T> v)
  {
    assert (v.size() >= n_channels);

    mp_foreach_idx<n_channels> ([&] (auto chnl) {
      T* chnl_head2 = _z + _head + (chnl * _size * 2) + _size;
      *chnl_head2   = v[chnl];
    });
  }
  //----------------------------------------------------------------------------
  // sparse input (all the leading/trailing functions TBD if required)
  void push (std::array<crange<const T>, n_channels> v, uint count = 1)
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
  void reset (crange<const T> kernel_mem, crange<T> delay_line_mem)
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
  static bool verify_coeffs (const crange<T> kernel, uint ratio)
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
  // inputs spread.
  std::array<T, n_channels> tick (std::array<crange<const T>, n_channels> in)
  {
    for (auto& r : in) {
      assert (r.size() >= _ratio);
    }
    _delay.push (in, _ratio);
    return _filter.tick (_delay);
  }
  //----------------------------------------------------------------------------
  // inputs interleaved
  std::array<T, n_channels> tick (const crange<T> in)
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
    _delays[1].push (in);
    for (auto& r : in) {
      r = r.shrink_head (1);
    }
    _delays[0].push (in, _ratio - 1);
    return tick_after_input();
  }
  //----------------------------------------------------------------------------
  // interleaved version
  std::array<T, n_channels> tick (crange<const T> in)
  {
    _delays[1].push (in);
    _delays[0].push (in.shrink_head (n_channels), _ratio - 1);
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
  void tick (std::array<crange<T>, n_channels> out, const crange<T> in)
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
  void tick (crange<T> out, const crange<T> in)
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

    _period      = rate_src;
    _kernel_size = taps;
    _max_samples = div_ceil (rate_tgt, rate_src);
    memory_reset (rate_tgt);

    assert ((cutoff_hz * 2.f) <= (float) tgt_srate);
    auto  ratio = (float) rate_tgt / (float) rate_src;
    float fc    = cutoff_hz / (float) tgt_srate;
    if (ratio < 1.f) {
      fc *= ratio;
    }

    std::vector<value_type> bigkernel {};
    bigkernel.resize (taps * rate_tgt);

    kaiser_lp_kernel_2<value_type> (
      bigkernel,
      fc / rate_tgt, // TODO: 0 to 0.5 or 0 to 1?
      kaiser_beta,
      rate_tgt,
      minphase);

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
      _n_spls_tbl[i] = n_spls;

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
  uint max_n_samples() const { return _max_samples; }
  //----------------------------------------------------------------------------
  // outputs and takes interleaved samples, "out" has to contain enough space
  // for "max_n_samples() * n_channels" elements.
  //
  uint tick (crange<value_type> out, const crange<value_type> in)
  {
    assert (in.size() >= n_channels);
    assert (out.size() >= (max_n_samples() * n_channels));

    auto n_spls = _n_spls_tbl[_n_in];
    // not using "push_leading/prepare_trailing" instead of "push" so "in" and
    // "out" can alias.
    _delay.push (in);

    for (uint i = 0; i < n_spls; ++i, ++_n_out) {
      detail::convolution_block<value_type, n_channels, alignment> dotprod;
      dotprod.reset (get_subkernel (_n_out));
      auto spls = dotprod.tick (_delay);
      memcpy (&out[i * n_channels], spls.data(), sizeof spls);
    }

    ++_n_in;
    if (_n_in >= _period) {
      _n_in = _n_out = 0;
    }
    return n_spls;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  crange<value_type> get_subkernel (uint pos)
  {
    return {(T*) &_kernels[_kernel_bytes * pos], _kernel_size};
  }
  //----------------------------------------------------------------------------
  void memory_reset (uint n_kernel_tables)
  {
    _kernel_bytes = round_ceil<uint> (_kernel_size * sizeof (T), alignment);

    uint mem_align        = std::max<uint> (alignment, sizeof (T));
    uint n_spls_tbl_bytes = round_ceil (_period, mem_align);
    uint del_elems        = n_channels * _kernel_size * 2;
    uint del_size_bytes = round_ceil<uint> (del_elems * sizeof (T), mem_align);
    uint k_size_bytes   = _kernel_bytes * n_kernel_tables;

    _mem.resize (n_spls_tbl_bytes + del_size_bytes + k_size_bytes);

    _n_spls_tbl = &_mem[0];
    _delay.reset (
      make_crange ((T*) &_mem[n_spls_tbl_bytes], del_elems), _kernel_size);
    _kernels = &_mem[n_spls_tbl_bytes + del_size_bytes];
  }
  //----------------------------------------------------------------------------
  void write_subkernel (
    crange<T> bigkernel,
    uint      offset,
    uint      pos,
    uint      rate_tgt)
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
  u8*                                           _n_spls_tbl;
  u8*                                           _kernels;
  uint                                          _period;
  uint                                          _kernel_size;
  uint                                          _kernel_bytes;
  uint                                          _max_samples;
  uint                                          _n_in;
  uint                                          _n_out;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// A polyphase decimator or interpolator taking care of power of two factors
// combined with a fractional resampler. Be sure to measure performance, as
// the latency is increased.
template <class T, uint channels = 1, uint alignment = alignof (T)>
class resampler {
public:
  //----------------------------------------------------------------------------
  using value_type                  = T;
  static constexpr uint n_channels  = channels;
  static constexpr uint max_samples = 128;
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

    uint small = std::min (tgt_srate, src_srate);
    uint big   = std::max (tgt_srate, src_srate);
    uint ratio = last_bit_set (big / small);
    // 0 is no bit set, 1 is a ratio of 1, 2 a ratio of 2, 3 a ratio of 4, etc
    ratio = (ratio >= 2) ? 1 << (ratio - 1) : 1;

    if (is_downsampler) {
      uint frac_srate = src_srate;
      _max_samples    = 1;

      if (ratio != 1) {
        frac_srate /= ratio;
        std::vector<value_type> tmp_kernel;
        tmp_kernel.resize (taps_int * ratio);

        auto fc = cutoff_hz / (float) src_srate;

        // TODO: fractional delay compensation?
        kaiser_lp_kernel_2<value_type> (
          tmp_kernel, fc, kaiser_beta, 0, minphase);

        _integer        = decimator_type {};
        auto& decim_rsc = std::get<decimator_type> (_integer);

        decim_rsc.decimator.reset (tmp_kernel, ratio);
        decim_rsc.input.reserve (ratio * n_channels);
        // Insert zeros: synchronize so the first sample input generates output.
        decim_rsc.input.resize ((ratio - 1) * n_channels);
      }

      if (tgt_srate != frac_srate) {
        _fractional.emplace();
        _fractional->reset (
          tgt_srate, frac_srate, taps_int, cutoff_hz, kaiser_beta, minphase);
        assert (_fractional->max_n_samples() == 1);
      }
    }
    else {
      uint frac_srate = tgt_srate;
      _max_samples    = ratio;

      if (ratio != 1) {
        frac_srate /= ratio;
        std::vector<value_type> tmp_kernel;
        tmp_kernel.resize (taps_int * ratio);

        auto fc = cutoff_hz / (float) tgt_srate;
        // TODO: fractional delay compensation?
        kaiser_lp_kernel_2<value_type> (
          tmp_kernel, fc, kaiser_beta, 0, minphase);

        _integer     = interpolator_type {};
        auto& interp = std::get<interpolator_type> (_integer);
        interp.reset (tmp_kernel, ratio);
      }

      if (src_srate != frac_srate) {
        _fractional.emplace();
        _fractional->reset (
          frac_srate, src_srate, taps_int, cutoff_hz, kaiser_beta, minphase);
        _max_samples *= _fractional->max_n_samples();
      }
    }
    // this uses the stack for intermediate memory, the implementation probably
    // neeeds dynamic memory for ratios that big.
    assert (_max_samples <= max_samples);
  }
  //----------------------------------------------------------------------------
  uint max_n_samples() const { return _max_samples; }
  //----------------------------------------------------------------------------
  uint tick_upsampler (crange<value_type> out, const crange<value_type> in)
  {
    assert (max_n_samples() > 1);
    assert (out.size() >= (max_n_samples() * n_channels));
    assert (in.size() >= n_channels);

    // the fractional resampler will output at most 2 samples, so we touch the
    // stack instead of the internal "_buffer" std::vector. Then "out" should
    // have space to hold all the data.
    std::array<value_type, 2 * n_channels> tmp;

    // handling the no-src change case.
    uint n_spls = 1;
    memcpy (tmp.data(), in.data(), sizeof tmp[0] * n_channels);

    if (_fractional) {
      n_spls = _fractional->tick (tmp, tmp);
      assert (n_spls <= 2);
    }

    if (std::holds_alternative<interpolator_type> (_integer)) {
      auto&     interpolator = std::get<interpolator_type> (_integer);
      crange<T> in           = tmp;
      auto      ratio        = interpolator.ratio();

      for (uint i = 0; i < n_spls; ++i) {
        interpolator.tick (
          out.cut_head (n_channels * ratio), in.cut_head (n_channels));
      }
      n_spls *= ratio;
    }
    else {
      memcpy (out.data(), tmp.data(), sizeof tmp[0] * n_channels * n_spls);
    }
    return n_spls;
  }
  //----------------------------------------------------------------------------
  uint tick_downsampler (crange<value_type> out, const crange<value_type> in)
  {
    assert (max_n_samples() == 1);
    assert (out.size() >= n_channels);
    assert (in.size() >= n_channels);

    // handling the no-src change case.
    uint n_spls = 1;
    memcpy (out.data(), in.data(), sizeof in[0] * n_channels);

    if (std::holds_alternative<decimator_type> (_integer)) {
      auto& dec_rsc = std::get<decimator_type> (_integer);
      dec_rsc.input.insert (
        dec_rsc.input.end(), in.begin(), in.begin() + n_channels);
      n_spls = 0;

      if (dec_rsc.input.size() == (dec_rsc.decimator.ratio() * n_channels)) {
        auto spls_arr = dec_rsc.decimator.tick (dec_rsc.input);
        dec_rsc.input.clear();
        memcpy (out.data(), spls_arr.data(), sizeof spls_arr);
        n_spls = 1;
      }
    }
    if (_fractional && n_spls == 1) {
      n_spls = _fractional->tick (out, out);
    }
    return n_spls;
  }
  //----------------------------------------------------------------------------
  uint tick (crange<value_type> out, const crange<value_type> in)
  {
    if (max_n_samples() == 1) {
      return tick_downsampler (out, in);
    }
    else {
      return tick_upsampler (out, in);
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
  uint                                                            _max_samples;
};
//------------------------------------------------------------------------------

} // namespace artv
