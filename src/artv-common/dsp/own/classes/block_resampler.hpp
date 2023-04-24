#pragma once

#include "artv-common/dsp/own/classes/circular_queue.hpp"
#include "artv-common/dsp/own/classes/fir.hpp"
#include "artv-common/dsp/own/classes/window.hpp"

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
//------------------------------------------------------------------------------
// Resampling at non-integer rates and processing blockwise has enough
// complexity to guarantee an abstraction. This class abstracts the SRC
// conversion and rate change of a FX.
template <class T, uint N_channels = 2>
class block_resampler {
public:
  static constexpr uint n_channels = N_channels;

  using value_type  = T;
  using sample_type = std::array<value_type, n_channels>;
  //----------------------------------------------------------------------------
  // "block_processor_fn" has to have the signature below:
  //
  // void (xspan<sample_type>)
  //
  // The block of samples to process is passed as a read and write contiguous
  // range.
  template <class U, class Functor>
  void process (
    xspan<U*>       outs,
    xspan<U const*> ins,
    uint            block_samples,
    Functor         block_processor_fn)
  {
    // Using the VLA extension (non-portable) to avoid touching a lot of memory
    // whenever it's possible. It may inhibit inlining but I guess it's not that
    // important at this level.
    sample_type stack_tgt_buff[_tgt_rate_bf ? 0 : _tgt_rate_bf_spls];
    sample_type stack_src_buff[_src_rate_bf ? 0 : _src_rate_bf_spls];

    auto tgt_rate_bf
      = xspan {_tgt_rate_bf ? _tgt_rate_bf : stack_tgt_buff, _tgt_rate_bf_spls};

    auto src_rate_bf
      = xspan {_src_rate_bf ? _src_rate_bf : stack_src_buff, _src_rate_bf_spls};

    uint spls_out  = 0;
    uint block_rem = block_samples;

    while (block_rem) {
      // ensure that the number of ticks never exceeds the internal buffer
      // sizes, so client code can happily declare its variables related to the
      // block size on the stack.
      auto block_spls_in = resampling::get_n_ticks_for_n_new_spls_floor (
        _converter_in.ratio(), _converter_in.corrected_pos(), _max_block_size);
      assert (block_spls_in && "block extremely small compared with the ratio");
      block_spls_in = std::min (block_rem, block_spls_in);
      block_rem -= block_spls_in;

      // interleave
      auto interleaved = src_rate_bf;
      for (uint i = 0; i < block_spls_in; ++i) {
        for (uint c = 0; c < n_channels; ++c) {
          interleaved[i][c] = ins[c][spls_out + i];
        }
      }
      // resample
      auto block     = tgt_rate_bf;
      uint blocksize = _converter_in.tick (
        block.cast (value_type {}),
        interleaved.cast (value_type {}),
        block_spls_in);
      assert (blocksize <= _max_block_size && "Bug!!!");

      xspan<sample_type> converted {};
      if (likely (blocksize > 0)) {
        //  process
        block = block.get_head (blocksize);
        block_processor_fn (block); // It this doesn't inline, then CRTP

        // resample
        converted        = src_rate_bf;
        uint n_converted = _converter_out.tick (
          converted.cast (value_type {}),
          block.cast (value_type {}),
          blocksize);
        converted = converted.get_head (n_converted);
      }

      // try to clear remaining samples from past runs
      uint allowed_spls_out = spls_out + block_spls_in;
      while (_remainder.size() && spls_out < allowed_spls_out) {
        auto spl = _remainder.pop();
        for (uint c = 0; c < n_channels; ++c) {
          outs[c][spls_out] = spl[c];
        }
        ++spls_out;
      }
      // generate as many output samples as input samples where fed
      while (spls_out < allowed_spls_out && converted.size()) {
        for (uint c = 0; c < n_channels; ++c) {
          outs[c][spls_out] = converted[0][c];
        }
        ++spls_out;
        converted.cut_head (1);
      }
      // save the new remainder (if any)
      while (converted.size()) {
        _remainder.push (converted[0]);
        converted.cut_head (1);
      }
    }
  }
  //----------------------------------------------------------------------------
  // "max_block_size": The maximum block length.
  //
  // "block_buffers_max_stack_bytes": The implementation tries to keep the cache
  //  hot by using the stack for as many intermediate buffers as possible. This
  //  sets a limit. This parameter is not defaulted to make clear that this is
  //  done.
  void reset (
    uint  tgt_srate,
    uint  src_srate,
    float in_cutoff_hz,
    float out_cutoff_hz,
    uint  n_taps_int,
    uint  n_taps_frac,
    float kaiser_att_db,
    bool  minphase,
    uint  max_block_size,
    uint  block_buffers_max_stack_bytes)
  {
    auto beta = kaiser_beta_estimate (kaiser_att_db);
    _converter_in.reset (
      tgt_srate,
      src_srate,
      n_taps_int,
      n_taps_frac,
      in_cutoff_hz,
      beta,
      minphase);
    _converter_out.reset (
      src_srate,
      tgt_srate,
      n_taps_int,
      n_taps_frac,
      out_cutoff_hz,
      beta,
      minphase);

    _max_block_size = max_block_size;

    uint remainder_n_spls;

    const uint max_stack_spls
      = block_buffers_max_stack_bytes / sizeof (sample_type);

    // precompute intermediate buffer sizes
    if (tgt_srate > src_srate) {
      // upsampler first
      uint ratio        = div_ceil (tgt_srate, src_srate);
      _tgt_rate_bf_spls = max_block_size;
#if 0
      // A downsampler always outputs 0 or 1 samples, subtracting one from the
      // ratio to so ratios like e.g. 1.001 get enough samples
      _src_rate_bf_spls = div_ceil (_tgt_rate_bf_spls, ratio - 1);
#else
      // The current implementation uses the output buffer as an input
      // deinterleaving buffer, so we set it to the desired block size. If this
      // requisite changes this #ifdef can be left on the other branch.
      _src_rate_bf_spls = max_block_size;
#endif
      // Having a remainder queue with a surplus of samples equal to the ratio.
      remainder_n_spls = ratio;
    }
    else {
      // downsampler first
      // A downsampler always outputs 0 or 1 samples, so the desired block size
      // will never be surpassed.
      _tgt_rate_bf_spls = max_block_size;
      uint ratio        = div_ceil (src_srate, tgt_srate);
      uint ratio_int    = src_srate / tgt_srate;
      auto ratio_frac   = _converter_out.ratio_frac();
      // in some stages the fractional resampler (max 1.999 ratio) will be
      // feeding 2 samples to the integer 2x-4x... upsampler. Other cycles one.
      // The periodicity has to be computed.
      uint n_spls2x     = ratio_frac.num - ratio_frac.den; // num: out, den: in
      uint in_cycle_len = ratio_frac.den;
      uint n_spls1x     = in_cycle_len - n_spls2x;

      uint full_cycles = max_block_size / in_cycle_len;
      uint cycle_rem   = max_block_size % in_cycle_len;

      _src_rate_bf_spls = full_cycles * ((n_spls2x * 2 + n_spls1x) * ratio_int);
      uint rem2x        = std::min (cycle_rem, n_spls2x);
      _src_rate_bf_spls += rem2x * 2 * ratio_int;
      cycle_rem -= rem2x;
      uint rem1x = std::min (cycle_rem, n_spls1x);
      _src_rate_bf_spls += rem1x * ratio_int;
      // Having a remainder queue with a surplus of samples equal to the ratio
      // rounded up.
      remainder_n_spls = ratio;
    }

    // ensure that we not allocate too big buffers on the stack
    uint dyn_tgt_rate_bf_spls = 0;
    uint dyn_src_rate_bf_spls = 0;

    if ((_tgt_rate_bf_spls + _src_rate_bf_spls) > max_stack_spls) {
      uint big   = std::max (_tgt_rate_bf_spls, _src_rate_bf_spls);
      uint small = std::min (_tgt_rate_bf_spls, _src_rate_bf_spls);

      if (big > max_stack_spls && small > max_stack_spls) {
        dyn_tgt_rate_bf_spls = _tgt_rate_bf_spls;
        dyn_src_rate_bf_spls = _src_rate_bf_spls;
      }
      else {
        assert (big > max_stack_spls);
        if (big == _tgt_rate_bf_spls) {
          dyn_tgt_rate_bf_spls = _tgt_rate_bf_spls;
        }
        else {
          dyn_src_rate_bf_spls = _src_rate_bf_spls;
        }
      }
    }
    // allocate and assign memory
    remainder_n_spls = pow2_round_ceil (remainder_n_spls);
    _raw_mem.clear();
    _raw_mem.resize (
      remainder_n_spls + dyn_tgt_rate_bf_spls + dyn_src_rate_bf_spls);
    auto dynmem = xspan {_raw_mem};

    _remainder.reset (dynmem.cut_head (remainder_n_spls));
    _tgt_rate_bf = dynmem.cut_head (dyn_tgt_rate_bf_spls).data();
    _src_rate_bf = dynmem.cut_head (dyn_src_rate_bf_spls).data();
  }
  //----------------------------------------------------------------------------
private:
  resampler<value_type, n_channels>       _converter_in;
  resampler<value_type, n_channels>       _converter_out;
  sample_type*                            _tgt_rate_bf;
  sample_type*                            _src_rate_bf;
  uint                                    _tgt_rate_bf_spls = 0;
  uint                                    _src_rate_bf_spls = 0;
  uint                                    _max_block_size   = 0;
  static_pow2_circular_queue<sample_type> _remainder;
  std::vector<sample_type>                _raw_mem;
};
//------------------------------------------------------------------------------
} // namespace artv
