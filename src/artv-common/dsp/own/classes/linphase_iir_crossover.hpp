#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/misc/delay_compensation_buffers.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/filters/composite/linear_iir_butterworth.hpp"

namespace artv {

// Linear phase
//------------------------------------------------------------------------------
template <uint N_Crossovers>
class linphase_iir_crossover {
public:
  static constexpr uint n_crossovers = N_Crossovers;
  static constexpr uint n_bands      = n_crossovers + 1;
  static constexpr uint max_order    = 8;

  static_assert (n_bands >= 2, "not a crossover");
  //----------------------------------------------------------------------------
  void reset (float samplerate, std::array<uint, n_crossovers> n_stages)
  {
    _mem.clear();

    xspan_memset (make_xspan (_order), 0);
    xspan_memset (make_xspan (_freq), 0);
    xspan_memset (make_xspan (_coeffs), 0);
    xspan_memset (make_xspan (_states), 0);
    xspan_memset (make_xspan (_in_delcomp), 0);
    xspan_memset (make_xspan (_out_delcomp), 0);

    _t_spl    = 1.f / samplerate;
    _n_stages = n_stages;

    // compute max requirements for crossovers
    std::array<uint, n_crossovers> stage_state; // states (samples)
    std::array<uint, n_crossovers> stage_in; // input/hp buffer size (samples)
    std::array<uint, n_crossovers> stage_out; // output buffer size (samples)

    uint total_in     = 0; // combined input/hp buffer sizes (samples)
    uint total_states = 0; // combined state sizes (samples)

    for (uint i = 0; i < n_crossovers; ++i) {
      uint n_states = lp_type::get_n_states (max_order, _n_stages[i]);
      uint latency  = lp_type::get_latency (max_order, _n_stages[i]);

      stage_state[i] = n_states;
      stage_in[i]    = latency;

      total_in += latency;
      total_states += n_states;
    }

    uint total_out    = 0; // combined output buffer size (samples)
    uint total_in_rem = total_in;

    for (uint i = 0; i < (n_crossovers - 1); ++i) {
      total_in_rem -= stage_in[i];
      stage_out[i] = total_in_rem;
      total_out += total_in_rem;
    }
    // not required on the last stage, it does exist to avoid conditionals.
    stage_out[n_crossovers - 1] = 0;

    // allocate all required memory.
    _mem.resize (block_mem::n_double_x2 + total_states + total_in + total_out);

    // assign the (worst-case) memory ranges
    double_x2* ptr = _mem.data();
    _block = reinterpret_cast<block_mem*> (ptr); // TODO: strict alias...
    ptr += block_mem::n_double_x2;

    for (uint i = 0; i < n_crossovers; ++i) {
      _states[i] = make_xspan (ptr, stage_state[i]);
      ptr += stage_state[i];

      _in_delcomp[i].reset (make_xspan (ptr, stage_in[i]), 0);
      ptr += stage_in[i];

      _out_delcomp[i].reset (make_xspan (ptr, stage_out[i]), 0);
      ptr += stage_out[i];
    }
  }
  //----------------------------------------------------------------------------
  uint get_latency() const
  {
    for (uint i = 0; i < n_crossovers; ++i) {
      if (_order[i] != 0) {
        // The first available crossover is the one having the information more
        // easily accessible.
        return _in_delcomp[i].delay() + _out_delcomp[i].delay();
      }
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  static uint get_n_stages (float freq, float t_spl, float snr_db)
  {
    return lp_type::get_n_stages (freq, t_spl, snr_db);
  }
  //----------------------------------------------------------------------------
  void set_crossover_point (
    uint  crossv_idx,
    float freq_l,
    float freq_r,
    uint  order) // order, 0 -> disabled, 2, 4, 8.
  {
    assert (crossv_idx < n_crossovers);
    uint idx = crossv_idx;

    if (
      _freq[idx][0] == freq_l && _freq[idx][1] == freq_r
      && order == _order[idx]) {
      return;
    }
    _freq[idx][0] = freq_l;
    _freq[idx][1] = freq_r;

    if (order != 0) {
      double_x2 f = {freq_l, freq_r};
      lp_type::reset_coeffs<double_x2> (_coeffs[idx], f, _t_spl, order);
    }

    if (order == _order[idx]) {
      return;
    }
    _order[idx] = order;
    // asuming an order change as able to cause a click.
    xspan_memset (_states[idx], 0);
    recompute_latencies (idx);
    memset (_block, 0, sizeof *_block);
  }
  //----------------------------------------------------------------------------
  std::array<double_x2, n_bands> tick (double_x2 in)
  {
    std::array<double_x2, n_bands> bands {};
    auto                           hp_out = in;

    for (uint c = 0; c < n_crossovers; ++c) {
      if (_order[c] == 0) {
        continue;
      }

      auto lp_out = lp_type::tick<double_x2> (
        _coeffs[c], _states[c], hp_out, _order[c], _n_stages[c], _sample_idx);

      auto delayed_in = _in_delcomp[c].exchange (hp_out);

      hp_out = delayed_in - lp_out; // Linear-phase: getting HP by subtraction.

      if (_out_delcomp[c].delay() > 0) {
        bands[c] = _out_delcomp[c].exchange (lp_out);
      }
      else {
        bands[c] = lp_out;
      }
    }
    ++_sample_idx;
    bands[n_crossovers] = hp_out;
    return bands;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (ins.size() >= 2);
    assert (ins[0] != nullptr);
    assert (ins[1] != nullptr);

    assert (outs.size() >= (n_bands * 2));
    for (uint b = 0; b < (n_bands * 2); ++b) {
      assert (outs[b] != nullptr);
    }

    // hp always contains the hp residue. The highpass out always sits on the
    // last channel independently on how many bands are activated.
    double_x2* hp = &_block->bands[n_crossovers][0];

    uint done = 0;
    while (done < samples) {
      uint blocksize = std::min (samples - done, block_mem::block_size);

      // to interleaved in
      for (uint i = 0; i < blocksize; ++i) {
        hp[i] = double_x2 {ins[0][done + i], ins[1][done + i]};
      }
      // process all bands
      for (uint c = 0; c < n_crossovers; ++c) {
        if (_order[c] == 0) {
          continue;
        }
        // process and place non latency compensated lp in band out
        memcpy (&_block->bands[c], hp, blocksize * sizeof hp[0]);
        lp_type::tick<double_x2> (
          _coeffs[c],
          _states[c],
          make_xspan (&_block->bands[c][0], blocksize),
          _order[c],
          _n_stages[c],
          _sample_idx);

        // compute hp for next band
        for (uint i = 0; i < blocksize; ++i) {
          auto lp         = _block->bands[c][i];
          auto delayed_hp = _in_delcomp[c].exchange (hp[i]);
          hp[i]           = delayed_hp - lp;
        }

        if (_out_delcomp[c].delay() == 0) {
          // no latency to compensate on this output
          continue;
        }
        for (uint i = 0; i < blocksize; ++i) {
          // latency compensate this band's output
          auto lp             = _block->bands[c][i];
          _block->bands[c][i] = _out_delcomp[c].exchange (lp);
        }
      }
      // block output deinterleave
      for (uint b = 0; b < n_bands; ++b) {
        uint l = b * 2;
        uint r = l + 1;
        for (uint i = 0; i < blocksize; ++i) {
          outs[l][done + i] = _block->bands[b][i][0];
          outs[r][done + i] = _block->bands[b][i][1];
        }
      }
      done += blocksize;
      _sample_idx += blocksize;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  struct block_mem {
    static constexpr uint block_size  = 64;
    static constexpr uint n_double_x2 = (block_size * n_bands);

    std::array<std::array<double_x2, block_size>, n_bands> bands;
  };
  //----------------------------------------------------------------------------
  void recompute_latencies (uint modified_idx)
  {
    auto idx       = modified_idx;
    auto new_order = _order[idx];
    if (new_order == 0) {
      _in_delcomp[idx].reset (0);
    }
    else {
      _in_delcomp[idx].reset (lp_type::get_latency (new_order, _n_stages[idx]));
    }
    // recompute input buffer total + memory clear
    uint total = 0;
    for (uint i = 0; i < n_crossovers; ++i) {
      auto del = _in_delcomp[i].delay();
      total += del;
      _in_delcomp[i].reset (del); // clear memory and position
    }
    // recompute output buffer totals + memory clear
    uint rem = total;
    for (uint i = 0; i < n_crossovers; ++i) {
      rem -= _in_delcomp[i].delay();
      _out_delcomp[i].reset (rem); // reset latency, clear memory and position
    }
  }
  //----------------------------------------------------------------------------
  struct buffering {
    uint in, out;
  };

  static constexpr uint n_chnls = 2;

  using lp_type = linear_iir_butterworth_lowpass_any_order;

  static constexpr uint n_crossv_coeffs = lp_type::get_n_coeffs (max_order);

  alignas (sse_bytes)
    std::array<std::array<double_x2, n_crossv_coeffs>, n_crossovers> _coeffs;

  std::array<xspan<double_x2>, n_crossovers>                    _states;
  std::array<delay_compensated_buffer<double_x2>, n_crossovers> _in_delcomp;
  std::array<delay_compensated_buffer<double_x2>, n_crossovers> _out_delcomp;

  std::vector<double_x2, overaligned_allocator<double_x2, sse_bytes>> _mem;
  block_mem*                                                          _block;

  std::array<uint, n_crossovers>                       _order {};
  std::array<uint, n_crossovers>                       _n_stages;
  std::array<std::array<float, n_chnls>, n_crossovers> _freq {};
  float                                                _t_spl;
  uint                                                 _sample_idx;
};
//------------------------------------------------------------------------------
} // namespace artv
