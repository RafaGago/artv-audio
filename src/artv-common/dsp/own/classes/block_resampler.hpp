#pragma once

#include "artv-common/dsp/own/classes/circular_queue.hpp"
#include "artv-common/dsp/own/classes/fir.hpp"
#include "artv-common/dsp/own/classes/window.hpp"

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

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
  using sample_type = vec<value_type, n_channels>;
  //----------------------------------------------------------------------------
  // "block_processor_fn" has to have the signature below:
  //
  // void (crange<sample_type>)
  //
  // The block of samples to process is passed as a read and write contiguous
  // range.
  template <class U, class Functor>
  void process (
    crange<U*>       outs,
    crange<U const*> ins,
    uint             block_samples,
    Functor          block_processor_fn)
  {
    // Using the VLA extension (non-portable) to avoid touching a lot of memory
    // whenever it's possible. It may inhibit inlining but I guess it's not that
    // important at this level.
    sample_type stack_in_buff[_in_buff ? 0 : _in_buff_spls];
    sample_type stack_out_buff[_out_buff ? 0 : _out_buff_spls];

    auto blockbf
      = make_crange (_in_buff ? _in_buff : stack_in_buff, _in_buff_spls);
    auto convertbf
      = make_crange (_out_buff ? _out_buff : stack_out_buff, _out_buff_spls);

    uint spls_out = 0;

    for (uint spls_in = 0; spls_in < block_samples;) {
      uint spls_in_past = spls_in;
      uint blocksize    = 0;
      auto block        = blockbf;

      // resample the input, if "_converter_in" is upsampling a fractiona ratio
      // the block size might go above the desired size (variably between
      // calls).
      while (spls_in < block_samples && blocksize < _desired_block_size) {
        sample_type spl;
        for (uint c = 0; c < n_channels; ++c) {
          spl[c] = ins[c][spls_in];
        }
        ++spls_in;
        uint n_spls = _converter_in.tick (
          block.cast (value_type {}),
          make_crange ((value_type*) &spl, n_channels));
        block.cut_head (n_spls);
        blocksize += n_spls;
      }

      block = make_crange (blockbf.data(), blocksize);
      block_processor_fn (block); // It this doesn't inline, then CRTP

      // resample the block
      auto convert     = convertbf;
      uint n_converted = 0;
      while (block.size()) {
        auto in     = block.cut_head (1);
        uint n_spls = _converter_out.tick (
          convert.cast (value_type {}), in.cast (value_type {}));
        convert.cut_head (n_spls);
        n_converted += n_spls;
      }

      // try to clear remaining samples
      uint block_spls_rem = spls_in - spls_in_past;
      uint spls_out_end   = spls_out + block_spls_rem;
      for (; _remainder.size() && spls_out < spls_out_end; ++spls_out) {
        auto spl = _remainder.pop();
        for (uint c = 0; c < n_channels; ++c) {
          outs[c][spls_out] = spl[c];
        }
      }
      // generate as many outs as possible
      uint conv_out = 0;
      while (spls_out < spls_out_end && conv_out < n_converted) {
        for (uint c = 0; c < n_channels; ++c) {
          outs[c][spls_out] = convertbf[conv_out][c];
        }
        ++spls_out;
        ++conv_out;
      }
      // save the new remainder (if any)
      for (; conv_out < n_converted; ++conv_out) {
        _remainder.push (convertbf[conv_out]);
      }
    }
  }
  //----------------------------------------------------------------------------
  // "desired_block_size": The average block length. When the first stage is an
  //   upsampler it might be exceeded by "ceil(ratio) - 1" samples.
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
    uint  desired_block_size,
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

    _desired_block_size = desired_block_size;

    uint remainder_n_spls;

    const uint max_stack_spls
      = block_buffers_max_stack_bytes / sizeof (sample_type);

    // precompute intermediate buffer sizes
    if (tgt_srate > src_srate) {
      // upsampler first
      uint ratio = div_ceil (tgt_srate, src_srate);
      // The upsampler will output some times more samples than the desired
      // block size, hence we account for them.
      _in_buff_spls = round_ceil (desired_block_size + ratio, ratio);
      // A downsampler always outputs 0 or 1 samples, subtracting one from the
      // ratio to so ratios like e.g. 1.001 get enough samples
      _out_buff_spls = div_ceil (_in_buff_spls, ratio - 1);
      // Having a remainder queue with a surplus of samples equal to the ratio.
      remainder_n_spls = ratio;
    }
    else {
      // downsampler first
      uint ratio = div_ceil (src_srate, tgt_srate);
      // A downsampler always outputs 0 or 1 samples, so the desired block size
      // will never be surpassed.
      _in_buff_spls = desired_block_size;
      // The upsampler will never output more samples than the scaled up ratio
      _out_buff_spls = _in_buff_spls * ratio;
      // Having a remainder queue with a surplus of samples equal to the ratio.
      remainder_n_spls = ratio;
    }

    // ensure that we not allocate too big buffers on the stack
    uint dyn_in_buff_spls  = 0;
    uint dyn_out_buff_spls = 0;

    if ((_in_buff_spls + _out_buff_spls) > max_stack_spls) {
      uint big   = std::max (_in_buff_spls, _out_buff_spls);
      uint small = std::min (_in_buff_spls, _out_buff_spls);

      if (big > max_stack_spls && small > max_stack_spls) {
        dyn_in_buff_spls  = _in_buff_spls;
        dyn_out_buff_spls = _out_buff_spls;
      }
      else {
        assert (big > max_stack_spls);
        if (big == _in_buff_spls) {
          dyn_in_buff_spls = _in_buff_spls;
        }
        else {
          dyn_out_buff_spls = _out_buff_spls;
        }
      }
    }
    // allocate and assign memory
    remainder_n_spls = pow2_round_ceil (remainder_n_spls);
    _raw_mem.clear();
    _raw_mem.resize (remainder_n_spls + dyn_in_buff_spls + dyn_out_buff_spls);
    auto dynmem = make_crange (_raw_mem);

    _remainder.reset (dynmem.cut_head (remainder_n_spls));
    _in_buff  = dynmem.cut_head (dyn_in_buff_spls).data();
    _out_buff = dynmem.cut_head (dyn_out_buff_spls).data();
  }
  //----------------------------------------------------------------------------
private:
  resampler<value_type, n_channels>       _converter_in;
  resampler<value_type, n_channels>       _converter_out;
  sample_type*                            _in_buff;
  sample_type*                            _out_buff;
  uint                                    _in_buff_spls       = 0;
  uint                                    _out_buff_spls      = 0;
  uint                                    _desired_block_size = 0;
  static_pow2_circular_queue<sample_type> _remainder;
  std::vector<sample_type>                _raw_mem;
};
//------------------------------------------------------------------------------
} // namespace artv
