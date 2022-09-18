#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

#include <juce_dsp/juce_dsp.h>

namespace artv {

//------------------------------------------------------------------------------
// applies a ramp on a buffer but instead of interpolating on each sample, the
// ramp is incremented at SIMD size width intervals.
template <class T>
void simd_gain_ramp (xspan<T> in, T gain_beg, T gain_end)
{
  using simd                  = juce::dsp::SIMDRegister<T>;
  constexpr size_t simd_elems = simd::SIMDNumElements;

  // applies a gain ramp in steps of the current SSE type size.
  T inc  = ((T) simd_elems) * (gain_end - gain_beg) / (T) (in.size());
  T gain = gain_beg;

  auto simd_ramp = [=, &gain] (std::array<T*, 1> c, uint blocks) {
    T* end = c[0] + (blocks * simd_elems);

    while (c[0] < end) {
      auto v = simd::fromRawArray (c[0]);
      auto g = simd {gain};
      v *= g;
      gain += inc;
      c[0] += simd_elems;
    }
  };

  auto regular_ramp = [&gain] (std::array<T*, 1> c, uint remainder) {
    T* end = c[0] + remainder;
    while (c[0] < end) {
      *c[0] *= gain;
      ++c[0];
    }
  };

  block_divide (
    simd::SIMDRegisterSize,
    std::array<T*, 1> {in.data()},
    in.size(),
    simd_ramp,
    regular_ramp);
}
//------------------------------------------------------------------------------
// applies gain ramps and sums all of them on in[0]. It doesn't modify in[1].
template <class T>
void simd_gain_ramps_add (
  std::array<T*, 2> in,
  uint              in_size,
  std::array<T, 2>  gain_beg,
  std::array<T, 2>  gain_end)
{
  using simd                  = juce::dsp::SIMDRegister<T>;
  constexpr size_t simd_elems = simd::SIMDNumElements;

  std::array<T, 2> inc;
  for (uint i = 0; i < inc.size(); ++i) {
    inc[i] = ((T) simd_elems) * (gain_end[i] - gain_beg[i]) / (T) (in_size);
  }
  auto gain = gain_beg;

  auto simd_ramp = [=, &gain] (std::array<T*, 2> c, uint blocks) {
    T* end = c[0] + (blocks * simd_elems);

    while (c[0] < end) {
      auto v1 = simd::fromRawArray (c[0]);
      auto v2 = simd::fromRawArray (c[1]);
      auto g1 = simd {gain[0]};
      auto g2 = simd {gain[1]};
      v1 *= g1;
      v2 *= g2;
      v1 += v2;
      v1.copyToRawArray (c[0]);
      gain[0] += inc[0];
      gain[1] += inc[1];
      c[0] += simd_elems;
      c[1] += simd_elems;
    }
  };

  auto regular_ramp = [&gain] (std::array<T*, 2> c, uint remainder) {
    T* end = c[0] + remainder;
    while (c[0] < end) {
      *c[0] = (*c[0] * gain[0]) + (*c[1] * gain[1]);
      ++c[0];
      ++c[1];
    }
  };

  block_divide (simd::SIMDRegisterSize, in, in_size, simd_ramp, regular_ramp);
}
//------------------------------------------------------------------------------
} // namespace artv
