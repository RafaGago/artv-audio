#pragma once

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <type_traits>

namespace artv {
//------------------------------------------------------------------------------
template <size_t N> // N == channel count
struct channels_mute_solo {
public:
  //----------------------------------------------------------------------------
  void reset (u64 v)
  {
    auto m = ~get_mute_solo_mask (v);
    auto s = get_mute_solo_mask (v >> 1);
    s      = s ? s : std::numeric_limits<decltype (s)>::max();
    mask   = m & s;
  }
  //----------------------------------------------------------------------------
  u64 get_mask() const { return mask; }
  //----------------------------------------------------------------------------
  bool is_channel_active (uint ch) const { return !!((mask >> (ch * N)) & 1); }

private
  : //----------------------------------------------------------------------------
  static u64 get_mute_solo_mask (u64 v)
  {
    // v &= (u64) 0x5555555555555555; // 0101010101... just get mute or solo
    u64 r = 0;
    // go from e.g. 01000101 to 0001000000010001
    for (int i = 0; i < N; ++i) {
      r |= (v & bit<u64> (i * 2)) << (i * (N - 2));
    }
    // multiply by the mask, so 0001000000010000 -> 1111000011110000, which
    // can then be directly masked to inputs/outputs.
    return r * lsb_mask<u64> (N);
  }
  //----------------------------------------------------------------------------
  u64 mask;
};
//------------------------------------------------------------------------------
template <size_t N> // N == channel count
class order_and_buffering {
public:
  using bus_latency_arr = std::array<u16, N>;

  //----------------------------------------------------------------------------
  template <class T>
  void recompute (
    u64  ins, // bit stream of size N * N representing receives of mix chnls
    u64  send_outs, // bit stream of size N * N representing sends of mix chnls
    u64  mute_solo_v, // bit stream of size N * 2 (2 = mute + solo)
    u64  mixer2mixer_sends, // bit stream of size N.
    u64  inhibit_reordering, // bit stream of size N.
    uint n_groups_log2,
    bus_latency_arr const& mix_latencies,
    T                      channel_has_processing)
  {
    uint n_parallel_mix = N >> n_groups_log2;
    assert (n_parallel_mix);

    // If there are parallel groups, we want to disable serial connections
    // between channels of different groups for easier computations.
    auto m2ms_mask = bit<u64> (n_parallel_mix - 1); // e.g 1000 if p.mix = 4
    for (uint i = 0; i < n_groups_log2; ++i) {
      m2ms_mask |= (m2ms_mask << (i * n_parallel_mix));
    }
    // e.g 10001000 if n_parallel_mix == 4 and N == 8

    mixer2mixer_sends &= ~m2ms_mask;
    u64 mixer2mixer_receives = mixer2mixer_sends << 1;

    memset (&in, 0, sizeof in);
    memset (&receives, 0, sizeof receives);
    memset (&out, 0, sizeof out);
    memset (&swaps, 0, sizeof swaps);

    for (uint i = 0; i < N; ++i) {
      mix[i]   = i + 1; // buffer ids are 1-indexed
      order[i] = i;
    }

    std::array<bool, N> has_processing;
    for (uint i = 0; i < has_processing.size(); ++i) {
      has_processing[i] = channel_has_processing (i);
    }

    // apply mute and solo to the input routings
    mute_solo.reset (mute_solo_v);
    ins &= mute_solo.get_mask();
    send_outs &= mute_solo.get_mask();

    // we get the sends for each mixer to which output, to detect buffer deps
    // is easier to operate with receives on the outputs from the mixers
    u64 outs = from_send_to_rcvs (send_outs);

    // The FX processor runs in three passes over the channels:
    // -1) direct channel swaps
    // -2) summing inputs and processing gain, pan, etc.
    // -3) summing processed channels into outputs (DAW buffer).

    for (uint chnl_idx = 0; chnl_idx < N; ++chnl_idx) {
      chnl_bits chnl = get_chnl_bits (chnl_idx);

      // Channel swaps of SINGLE cross channel references. Has to be done
      // first. Swaps are not affected by mixer2mixer sends.
      bool bus_swap_opt_enabled = n_groups_log2 == 0;
      if (
        bus_swap_opt_enabled && is_pow2 (ins & chnl.bit_in_all)
        && ((ins & chnl.bit_in_self) == 0)) {
        // only one remote channel contains this chnl as input. Maybe
        // swap...
        auto swap_idx = first_bit_set (ins & chnl.bit_in_all);
        assert (swap_idx > 0);
        swap_idx -= 1; // "first_bit_set" uses 1 based indexing
        swap_idx /= N;

        chnl_bits swap = get_chnl_bits (swap_idx);

        if ((ins & swap.bit_in_all) == (swap.bit_in_all & chnl.mask)) {
          // only this channel has the remote channel as input. Good
          // to swap.
          swaps[chnl_idx] = swap_idx + 1; // buffers are 1 indexed
          // update "ins" so the swap is not seen on the next steps.
          ins &= ~chnl.bit_in_all;
          ins |= chnl.bit << chnl.offset;
          ins &= ~swap.bit_in_all;
          ins |= swap.bit << swap.offset;
        }
      }

      // Fill inputs and output buffer data. Important that it is done
      // after swapping.
      for (uint j = 0; j < N; ++j) {
        if ((ins >> chnl.offset) & (1 << j)) {
          in[chnl_idx][j] = j + 1;
        }
        if ((outs >> chnl.offset) & (1 << j)) {
          out[chnl_idx][j] = j + 1;
        }
      }

      // Try to find a channel ordering that tries to minimize internal
      // buffer requirements (touching more memory).
      uint order_start_pos
        = std::find (order.begin(), order.end(), chnl_idx) - order.begin();
      uint order_new_pos = order_start_pos;

      // the channel's DAW buffer is modified if there is input summing
      bool chnl_is_modified = (ins & chnl.mask) != chnl.bit_in_self;
      // the channel's DAW buffer is modified if there is processing
      chnl_is_modified |= has_processing[chnl_idx];

      for (uint curr = order_start_pos + 1; curr < N; ++curr) {
        if (
          (chnl.bit_in_self & mixer2mixer_sends) != 0
          || (chnl.bit_in_self & mixer2mixer_receives) != 0) {
          // this channel sends or receives from its neighboring mixer, so we
          // inhibit the reordering optimization. They need to be processed
          // sequentially left to right. Not done outside of the loop to not
          // create another indentation level.
          break;
        }
        if ((chnl.bit_in_self & inhibit_reordering) != 0) {
          // manually specified to not reorder this channel.
          break;
        }

        uint chnl_bit_in_current = chnl.bit << (order[curr] * N);

        if (chnl_is_modified && (ins & chnl_bit_in_current)) {
          // someone depends on this channel on the inputs. place afterwards.
          order_new_pos = curr;
          continue;
        }
        // On the output summing stage the DAW buffer for this channel
        // is only modified if there is summing with other channels or
        // when the channel itself is not present (memset 0).
        if ((outs & chnl.mask) == chnl.bit_in_self) {
          // no modification on the outputs.
          continue;
        }

        if (outs & chnl_bit_in_current) {
          // someone depends on this channel on the inputs. place afterwards.
          order_new_pos = curr;
          continue;
        }
      }

      if (order_start_pos == order_new_pos) {
        continue; // no changes in ordering
      }
      auto ob = order.begin();
      std::move (
        ob + order_start_pos + 1, ob + order_new_pos + 1, ob + order_start_pos);
      order[order_new_pos] = chnl_idx;
    }
    // redo ins and outs to match the order, so we can detect dependencies
    // with just the ">" operator and masking.
    ins  = update_bitfields_for_current_order (ins);
    outs = update_bitfields_for_current_order (outs);

    for (uint i = 0; i < N; ++i) {
      u64       chnl_idx       = order[i];
      chnl_bits chnl           = get_chnl_bits (chnl_idx);
      u64       ro_mask        = lsb_mask<u64> (N) << (i * N); // reordered mask
      u64       ro_bit_in_self = chnl.bit << (i * N); // reordered

      // Detect intermediate buffer needs.
      if (((ins & ~ro_mask) & chnl.bit_in_all) > ro_mask) {
        // Someone processed later depends on the passed input, we
        // might not be able to overwrite the buffer passed by the DAW.
        if ((ins & ro_mask) != ro_bit_in_self || has_processing[chnl_idx]) {
          // there is either summing or processing that will destroy
          // the passed DAW buffer's contents. Requires internal
          // buffer (negative idx)
          mix[chnl_idx] = -(chnl_idx + 1);
        }
      }
      if (((outs & ~ro_mask) & chnl.bit_in_all) > ro_mask) {
        // Someone processed later depends on this mix buffer, we
        // might not be able to use the buffer passed by the DAW as mix
        // buffer.
        if ((outs & ro_mask) != ro_bit_in_self) {
          // there is summing that will destroy the mix buffer's
          // contents. Requires internal buffer (negative idx)
          mix[chnl_idx] = -(chnl_idx + 1);
        }
      }
    }

    for (uint chnl = 0; chnl < N; ++chnl) {
      // update outs
      for (uint i = 0; i < N; ++i) {
        if (out[chnl][i] != 0) {
          out[chnl][i] = mix[i];
        }
      }
      // update neighboring mixer inputs
      if (((1 << chnl) & mixer2mixer_receives) != 0) {
        // receives from neighboring mixer
        assert (chnl > 0 && "The leftmost channel has no left neighbor mixer");
        receives[chnl] = mix[chnl - 1];
      }
    }

    // Grouping-aware channel reordering. When processing, "outs" will need to
    // be transformed to contain only the buses relevant to each group.
    uint bus_beg  = 0;
    uint wr_idx   = 0;
    auto order_cp = order;
    do {
      auto bus_end = bus_beg + n_parallel_mix;

      for (uint i = 0; i < order.size(); ++i) {
        auto bus = order_cp[i];
        if (bus >= bus_beg && bus < bus_end) {
          order[wr_idx] = bus;
          ++wr_idx;
        }
      }

      bus_beg = bus_end;
    } while (bus_beg < order.size());

    compute_latencies (
      mix_latencies,
      n_parallel_mix,
      ins,
      send_outs,
      mixer2mixer_sends,
      mixer2mixer_receives);
  }
  //----------------------------------------------------------------------------
  // just a variant adding default parameters, as this can't be done on the
  // original function because of the lamda
  void recompute (
    u64  ins, // bit stream of size N * N representing receives of mix chnls
    u64  send_outs, // bit stream of size N * N representing sends of mix chnls
    u64  mute_solo_v, // bit stream of size N * 2 (2 = mute + solo)
    u64  mixer2mixer_sends, // bit stream of size N.
    u64  inhibit_reordering, // bit stream of size N.
    uint n_groups_log2                   = 0,
    bus_latency_arr const& mix_latencies = bus_latency_arr {})
  {
    recompute (
      ins,
      send_outs,
      mute_solo_v,
      mixer2mixer_sends,
      inhibit_reordering,
      n_groups_log2,
      mix_latencies,
      [] (uint chnl) { return true; });
  }
  //----------------------------------------------------------------------------
  static_assert (N <= 127, "overflowing signed char");

  // input buffers to each mixer: 1 based indexes. It contains N+1 entries
  // because they can input from its neighboring mixer.
  //
  // Notice that each buffer ID maps to a bus.
  //
  // // input buffers to each output channel: 1 based indexes

  std::array<std::array<s8, N>, N> in; // ins from the DAW
  std::array<s8, N>                receives; // rcv from neighbouring mixer
  std::array<std::array<s8, N>, N> out;
  std::array<s8, N>                mix; // buffers: 1 based indexes
  std::array<s8, N>                swaps; // buffers: 1 based indexes
  std::array<u8, N>                order; // channels: 0 based indexes

  channels_mute_solo<N> mute_solo;

  bus_latency_arr post_input_mix_delay_samples;
  bus_latency_arr fx_delay_samples;
  bus_latency_arr pre_output_mix_delay_samples;
  uint            plugin_latency;

private:
  //----------------------------------------------------------------------------
  u64 update_bitfields_for_current_order (u64 v)
  {
    u64 ret = 0;
    for (uint i = 0; i < order.size(); ++i) {
      ret |= ((v >> (order[i] * N)) & lsb_mask<u64> (N)) << (i * N);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  struct chnl_bits {
    u64 bit;
    u64 mask;
    u64 offset;
    u64 bit_in_self;
    u64 bit_in_all;
    u64 bit_in_others;
  };
  //----------------------------------------------------------------------------
  static constexpr u64 get_bit_selection_mask (k_int<4>)
  {
    return 0x1111; // N == 4, channels are selected with this mask.
  }
  //----------------------------------------------------------------------------
  static constexpr u64 get_bit_selection_mask (k_int<8>)
  {
    return 0x0101010101010101; // N == 8, channels are selected with this mask.
  }
  //----------------------------------------------------------------------------
  static chnl_bits get_chnl_bits (u64 chnl)
  {
    chnl_bits r;
    r.bit         = bit<u64> (chnl);
    r.offset      = chnl * N;
    r.mask        = lsb_mask<u64> (N) << r.offset;
    r.bit_in_self = r.bit << r.offset;
    static_assert (N == 8 || N == 4, "bit hack below has to be redone");
    u64 selmask     = get_bit_selection_mask (k_int<N> {});
    r.bit_in_all    = r.bit * selmask;
    r.bit_in_others = r.bit_in_all & ~r.mask;
    return r;
  }
  //----------------------------------------------------------------------------
  u64 from_send_to_rcvs (u64 v)
  {
    static_assert ((N * N) <= (sizeof v * 8), "");
    u64 r       = 0;
    u64 selmask = get_bit_selection_mask (k_int<N> {});
    for (uint i = 0; i < N; ++i) {
      u64 selected = (v >> i) & selmask;
      // from e.g. 0001000100010001 to 1111
      for (uint j = 0; j < (N - 1); ++j) {
        selected |= selected >> (N - 1);
      }
      r |= (selected & lsb_mask<u64> (N)) << (N * i);
    }
    return r;
  }
  //----------------------------------------------------------------------------
  void compute_latencies (
    bus_latency_arr const& mix_latencies,
    uint                   n_parallel_mix,
    u64                    ins_as_bits,
    u64                    outs_as_bits,
    u64                    mixer2mixer_sends,
    u64                    mixer2mixer_receives)
  {
    // Summarizing the relevant parts only. Each channel/bus has on this order:
    //
    // 1. A mixer for its inputs
    // 2. A post input-mix delay line.
    // 3. A receive for the neighboring (to the left) channel send (if enabled).
    // 4. FX, with latency, oversampling, etc. The "fx_z_spls"
    //    parameter
    // 5. A send to the neighboring channel (to the right) (if enabled).
    //-6. A pre output-mix delay line.
    //-7  An output mixer.
    //
    // This function is about computing the delays on 2 and 6.

    this->post_input_mix_delay_samples = bus_latency_arr {};
    this->pre_output_mix_delay_samples = bus_latency_arr {};
    this->fx_delay_samples             = bus_latency_arr {};
    this->plugin_latency               = 0;

    // Cummulative latencies at point 5. See leading comment.
    bus_latency_arr post_fx_cummulative_dly = {};

    // notice that this classification and dead channel detection could be done
    // at the main function. By not detecting dead channels there I guess that
    // it's easier to avoid audio clicks when automating, as the audio was
    // still being processed.
    //
    // In the case of latency there is no room for miscalculation.
    std::array<u8, N> chnl_io = classify_io_detect_dead_ends (
      ins_as_bits, outs_as_bits, mixer2mixer_sends, mixer2mixer_receives);

    for (uint c = 0; c < N; ++c) {
      if (chnl_io[c] == dead_channel) {
        continue;
      }

      this->fx_delay_samples[c]  = mix_latencies[c];
      post_fx_cummulative_dly[c] = mix_latencies[c];

      u64  chnl_mask      = lsb_mask<u64> (N) << (c * N);
      bool has_daw_inputs = !!(ins_as_bits & chnl_mask);
      bool has_receives   = mixer2mixer_receives & bit<uint> (c);

      if (chnl_io[c] & has_recvs_bit) {
        // "mixer2mixer_receives" has always zeroes on the first element of a
        // group. It is safe to access the previous neighbor.
        uint rcv_delay_samples = post_fx_cummulative_dly[c - 1];
        post_fx_cummulative_dly[c] += rcv_delay_samples;
        // no imput compensation if there is no input.
        this->post_input_mix_delay_samples[c] = rcv_delay_samples;
        this->post_input_mix_delay_samples[c] *= has_daw_inputs;
      }
    }

    for (uint i = 0; i < (N / n_parallel_mix); ++i) {
      uint beg = (i * n_parallel_mix);
      uint end = beg + n_parallel_mix;

      uint max_group_latency = *std::max_element (
        &post_fx_cummulative_dly[beg], &post_fx_cummulative_dly[end]);

      for (uint c = beg; c < end; ++c) {
        if (chnl_io[c] != dead_channel) {
          this->pre_output_mix_delay_samples[c]
            = max_group_latency - post_fx_cummulative_dly[c];
        }
      }
      this->plugin_latency += max_group_latency;
    }
  }
  //----------------------------------------------------------------------------
  enum channel_classification {
    dead_channel  = 0,
    has_ins_bit   = 1 << 0,
    has_outs_bit  = 1 << 1,
    has_sends_bit = 1 << 2,
    has_recvs_bit = 1 << 3,
  };
  //----------------------------------------------------------------------------
  std::array<u8, N> classify_io_detect_dead_ends (
    u64 ins_as_bits,
    u64 outs_as_bits,
    u64 mixer2mixer_sends,
    u64 mixer2mixer_receives)
  {
    // this function requires "mixer2mixer_sends" and "mixer2mixer_receives" to
    // be corrected for grouping.

    std::array<u8, N> chnl_io {};

    // first channel catalogation
    for (uint c = 0; c < N; ++c) {
      u64  chnl_mask = lsb_mask<u64> (N) << (c * N);
      bool has_ins   = !!(ins_as_bits & chnl_mask);
      bool has_outs  = !!(outs_as_bits & chnl_mask);
      bool has_sends = mixer2mixer_sends & bit<uint> (c);
      bool has_recvs = mixer2mixer_receives & bit<uint> (c);

      chnl_io[c] = has_ins ? has_ins_bit : 0;
      chnl_io[c] |= has_outs ? has_outs_bit : 0;
      chnl_io[c] |= has_sends ? has_sends_bit : 0;
      chnl_io[c] |= has_recvs ? has_recvs_bit : 0;

      // dead by no io
      bool dead   = ((!has_ins && !has_recvs) || (!has_outs && !has_sends));
      uint prev_c = (c == 0) ? 0 : c - 1; // no underflow reads of prev chnl.
      // dead by dead receive
      dead |= (!has_ins && has_recvs) && (chnl_io[prev_c] == 0);

      chnl_io[c] *= dead ? 0 : 1;
    }
    // detecting dead channels because of dead sends requires the subsequent
    // channels to have been classified (so another iteration) and to iterate in
    // reverse to correctly propagate.
    static_assert (N >= 2, "");
    for (uint c = (N - 2); c < N; --c) {
      bool has_outs            = chnl_io[c] & has_outs_bit;
      bool has_sends           = chnl_io[c] & has_sends_bit;
      bool dead_right_neighbor = chnl_io[c + 1] == dead_channel;
      bool dead                = !has_outs && has_sends && dead_right_neighbor;

      chnl_io[c] *= dead ? 0 : 1;
    }
    return chnl_io;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
} // namespace artv
