#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/composite/linkwitz_riley.hpp"

namespace artv {

// TODO(maybe): Remove template parameters to avoid possible code bloat if two
// crossovers of the same type live on the same program.
// TODO(maybe): When needed, if dry-wet processing is needed it is trivial to
//   process the dry signal with the allpass equivalent version of the
//   crossover. It is just running the dry signal through all allpasses.
//------------------------------------------------------------------------------
template <uint N_Crossovers, bool Controls_summing = true>
class crossover {
public:
  static constexpr uint n_crossovers     = N_Crossovers;
  static constexpr uint n_bands          = n_crossovers + 1;
  static constexpr bool controls_summing = Controls_summing;

  static_assert (n_bands >= 2, "not a crossover");
  //----------------------------------------------------------------------------
  void reset (float samplerate)
  {
    memset (this, 0, sizeof *this);
    _samplerate = samplerate;
  }
  //----------------------------------------------------------------------------
  void set_crossover_point (
    uint  crossv_idx,
    float freq_l,
    float freq_r,
    uint  order) // order, 0 -> disabled, 2, 4 8.
  {
    assert (crossv_idx < n_crossovers);
    uint idx = id_invert (crossv_idx);

    if (
      _freq[idx][0] == freq_l && _freq[idx][1] == freq_r
      && order == _order[idx]) {
      return;
    }
    _freq[idx][0] = freq_l;
    _freq[idx][1] = freq_r;

    if (order != 0) {
      double_x2 f = {freq_l, freq_r};
      lr::reset_coeffs (_coeffs[idx], f, _samplerate, order);
    }

    if (order == _order[idx]) {
      return;
    }
    _order[idx] = order;
    // reset states
    memset (&_states[idx], 0, sizeof _states[0]);

    // reset correction filters
    if (crossv_idx == 0) {
      return; // no correction filters for the first crossover (Highest Freq)
    }
    if constexpr (controls_summing) {
      auto correction_idx = crossv_idx - 1;
      memset (&_corr_states[correction_idx], 0, sizeof _corr_states[0]);
    }
    else {
      // Is there a better mathematical way than looping?
      uint correction_idx = 0;
      for (uint crv = 1; crv < (crossv_idx + 1); ++crv) {
        if (crv != crossv_idx) {
          correction_idx += crv;
          continue;
        }
        memset (&_corr_states[correction_idx], 0, crv * sizeof _corr_states[0]);
      }
    }
  }
  //----------------------------------------------------------------------------
  std::array<double_x2, n_bands> tick (double_x2 in)
  {
    // bands enable disabled.
    std::array<double_x2, n_bands> bands {};
    lr_crossover_out<double_x2>    crossv {};
    crossv.lp = in;

    // crossovers go top to bottom internally
    //
    //      c3   c2   c1   c0
    //  ------------------------
    // | b0 | b1 | b2 | b3 | b4 |
    //  ------------------------

    for (uint crv = 0, bnd = (n_bands - 1); crv < n_crossovers; ++crv, --bnd) {
      if (_order[crv] == 0) {
        continue;
      }
      crossv = lr::tick (_coeffs[crv], _states[crv], crossv.lp, _order[crv]);

      bands[bnd] = crossv.hp;
    }
    bands[0] = crossv.lp;
    return bands;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void tick (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    tick<T, !controls_summing> (outs, ins, samples);
  }
  //----------------------------------------------------------------------------
  std::array<double_x2, n_bands> phase_correct (
    std::array<double_x2, n_bands> bands)
  {
    // with three bands there are the same number of correctiong blocks wether
    // we can perform the summing optimizations or not.
    assert (!controls_summing || n_bands == 3);

    std::array<double_x2, n_bands> ret            = bands;
    uint                           correction_idx = 0;

    //        c3     c2     c1     c0
    //  ----------------------------------
    // |  b0  |  b1  |  b2  |  b3  |  b4  |
    //  ----------------------------------
    //                              ap-c1
    //                       ap-c2  ap-c2
    //                ap-c3  ap-c3  ap-c3

    for (uint crv = 1; crv < n_crossovers; ++crv) {
      for (uint i = 0; i < crv; ++i, ++correction_idx) {
        if (unlikely (_order[crv] == 0)) {
          continue;
        }
        uint band = n_crossovers - i;
        ret[band] = lr::apply_correction (
          _coeffs[crv], _corr_states[correction_idx], ret[band], _order[crv]);
      }
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  double_x2 sum (std::array<double_x2, n_bands> bands)
  {
    assert (controls_summing);

    double_x2 ret {};
    for (uint band = (n_bands - 1), crossv = 1; band > 1; --band, ++crossv) {
      uint correction_idx = crossv - 1;

      ret += bands[band];
      if (unlikely (_order[crossv] != 0)) {
        continue;
      }
      // apply correction to the sum of bands simultaneously
      ret = lr::apply_correction (
        _coeffs[crossv], _corr_states[correction_idx], ret, _order[crossv]);
    }
    ret += bands[0] + bands[1];
  }
  //----------------------------------------------------------------------------
  template <class T>
  void sum (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    static constexpr uint blocksize = 32;

    assert (outs.size() >= 2);
    assert (outs[0] != nullptr);
    assert (outs[1] != nullptr);

    assert (ins.size() >= (n_bands * 2));
    for (uint b = 0; b < (n_bands * 2); ++b) {
      assert (ins[b] != nullptr);
    }

    std::array<std::array<double_x2, n_bands>, blocksize> buffer;

    uint done = 0;
    while (done < samples) {
      uint subblock_size = std::min (samples - done, blocksize);

      for (uint b = 0; b < n_bands; ++b) {
        uint in = b * 2;
        for (uint i = 0; i < subblock_size; ++subblock_size) {
          buffer[i][b] = double_x2 {ins[in][done + i], ins[in + 1][done + i]};
        }
      }

      for (uint i = 0; i < subblock_size; ++subblock_size) {
        buffer[i][0] = tick (buffer[i]);
      }

      for (uint i = 0; i < subblock_size; ++subblock_size) {
        outs[0][done + i] = buffer[i][0][0];
        outs[1][done + i] = buffer[i][0][1];
      }
      done += subblock_size;
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T, bool do_phase_correct>
  void tick (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    static constexpr uint blocksize = 32;

    assert (ins.size() >= 2);
    assert (ins[0] != nullptr);
    assert (ins[1] != nullptr);

    assert (outs.size() >= (n_bands * 2));
    for (uint b = 0; b < (n_bands * 2); ++b) {
      assert (outs[b] != nullptr);
    }

    std::array<std::array<double_x2, n_bands>, blocksize> buffer;

    uint done = 0;
    while (done < samples) {
      uint subblock_size = std::min (samples - done, blocksize);

      for (uint i = 0; i < subblock_size; ++i) {
        buffer[i] = tick (double_x2 {ins[0][done + i], ins[1][done + i]});
      }
      if constexpr (do_phase_correct) {
        for (uint i = 0; i < subblock_size; ++i) {
          buffer[i] = phase_correct (buffer[i]);
        }
      }
      for (uint b = 0; b < n_bands; ++b) {
        for (uint i = 0; i < subblock_size; ++i) {
          uint out                = b * 2;
          outs[out][done + i]     = buffer[i][b][0];
          outs[out + 1][done + i] = buffer[i][b][1];
        }
      }
      done += subblock_size;
    }
  }
  //----------------------------------------------------------------------------
  using lr = linkwitz_riley_stage_any_order;
  //----------------------------------------------------------------------------
  uint id_invert (uint v)
  {
    // Crossover order is reversed. Splitting high to low freq but with the
    // memory layout increasing.
    assert (v < n_crossovers);
    return n_crossovers - 1 - v;
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_n_correctors()
  {
    return get_n_correctors (
      k_uint<0> {}, k_uint<0> {}, k_uint<n_crossovers> {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint Val, uint Stop>
  static constexpr auto get_n_correctors (
    k_uint<Idx>,
    k_uint<Val>,
    k_uint<Stop>)
  {
    return get_n_correctors (
      k_uint<Idx + 1> {}, k_uint<Val + Idx> {}, k_uint<Stop> {});
  }
  //----------------------------------------------------------------------------
  template <uint Result, uint Stop>
  static constexpr auto get_n_correctors (
    k_uint<Stop>,
    k_uint<Result>,
    k_uint<Stop>)
  {
    return k_uint<Result> {};
  }
  //----------------------------------------------------------------------------
  static_assert (
    decltype (get_n_correctors (
      k_uint<0> {},
      k_uint<0> {},
      k_uint<4> {}))::value
      == 6,
    "Broken get_n_correctors implementation...");
  //----------------------------------------------------------------------------
  static constexpr uint n_correctors = controls_summing
    ? n_crossovers - 1
    : decltype (get_n_correctors())::value;

  //----------------------------------------------------------------------------
  static constexpr uint n_chnls             = 2;
  static constexpr uint n_crossv_coeffs     = n_chnls * lr::n_coeffs;
  static constexpr uint n_crossv_states     = n_chnls * lr::n_states;
  static constexpr uint n_correction_states = n_chnls * lr::n_correction_states;

  alignas (double_x2)
    std::array<std::array<double, n_crossv_coeffs>, n_crossovers> _coeffs;
  alignas (double_x2)
    std::array<std::array<double, n_crossv_states>, n_crossovers> _states;
  alignas (double_x2) std::
    array<std::array<double, n_correction_states>, n_correctors> _corr_states;

  std::array<uint, n_crossovers>                 _order {};
  std::array<std::array<float, 2>, n_crossovers> _freq {};
  float                                          _samplerate;
};
//------------------------------------------------------------------------------
} // namespace artv
