#pragma once

#include <cmath>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/parts/filters/composite/butterworth.hpp"

namespace artv {

// An IIR crossover with perfect reconstruction (linear-phase), but missing
// phase coherence between bands, exact crossover point on the left band
// side, and a flat frequency response on the passband.
//
// It is perfect (TODO: verify) with one pole crossovers, but gets more "wonky"
// as the number of poles increases.
//
// It is using Butterworth just so at least one of the sides has a predictable
// cutoff -3dB point. It still reconstructs perfectly with any type of
// crossover.
//------------------------------------------------------------------------------
template <uint N_Crossovers>
class wonky_crossover {
public:
  static constexpr uint n_crossovers = N_Crossovers;
  static constexpr uint n_bands      = n_crossovers + 1;

  static_assert (n_bands >= 2, "not a crossover");
  //----------------------------------------------------------------------------
  void reset (float samplerate)
  {
    memset (this, 0, sizeof *this);
    _t_spl = 1.f / samplerate;
  }
  //----------------------------------------------------------------------------
  void set_crossover_point (
    uint  crossv_idx,
    float freq_l,
    float freq_r,
    uint  order) // order, 0 -> disabled, 2, 4 8.
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
      crossv::reset_coeffs<double_x2> (_coeffs[idx], f, _t_spl, order);
    }

    if (order == _order[idx]) {
      return;
    }
    _order[idx] = order;
    // reset states
    memset (&_states[idx], 0, sizeof _states[0]);
  }
  //----------------------------------------------------------------------------
  std::array<double_x2, n_bands> tick (double_x2 in)
  {
    // the lowest band, usually the one with the most energy on audio signals,
    // is the one that doesn't have any wonky response.
    std::array<double_x2, n_bands> bands {};
    double_x2                      residue = in;

    for (uint i = 0; i < n_crossovers; ++i) {
      if (_order[i] == 0) {
        continue;
      }
      bands[i]
        = crossv::tick<double_x2> (_coeffs[i], _states[i], residue, _order[i]);
      residue -= bands[i];
    }
    bands[n_crossovers] = residue;
    return bands;
  }
  //----------------------------------------------------------------------------
  template <class T>
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
private:
  //----------------------------------------------------------------------------
  using crossv                    = butterworth_lowpass_any_order;
  static constexpr uint max_order = 8;
  static_assert (max_order < crossv::max_order, "");

  static constexpr uint n_crossv_coeffs
    = crossv::n_coeffs_for_order (max_order);
  static constexpr uint n_crossv_states
    = crossv::n_states_for_order (max_order);

  alignas (double_x2)
    std::array<std::array<double_x2, n_crossv_coeffs>, n_crossovers> _coeffs;
  alignas (double_x2)
    std::array<std::array<double_x2, n_crossv_states>, n_crossovers> _states;

  std::array<uint, n_crossovers>                 _order {};
  std::array<std::array<float, 2>, n_crossovers> _freq {};
  float                                          _t_spl;
};
//------------------------------------------------------------------------------
} // namespace artv
