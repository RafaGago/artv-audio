#pragma once
#include <cassert>
#include <cstdint>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// applies a ramp on a buffer but instead of interpolating on each sample, the
// ramp is incremented at SIMD size width intervals.
template <class T>
void simd_gain_ramp (xspan<T> in, T gain_beg, T gain_end)
{
  static_assert (std::is_floating_point_v<T>);
  using vec_t           = vec16<T>;
  constexpr auto traits = vec_traits_t<vec_t> {}; // matching SSE
  assert (0 == ((std::uintptr_t) in.data() % 16)); // assert alignment

  // applies a gain ramp in steps of the current SSE type size.
  T inc  = ((T) traits.size) * (gain_end - gain_beg) / (T) (in.size());
  T gain = gain_beg;

  auto simd_ramp = [=, &gain] (std::array<T*, 1> c, uint blocks) {
    T* end = c[0] + (blocks * traits.size);

    while (c[0] < end) {
      auto v = vec_load<vec_t> (c[0]);
      v *= gain;
      gain += inc;
      vec_store (c[0], v);
      c[0] += traits.size;
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
    traits.bytes,
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
  static_assert (std::is_floating_point_v<T>);
  using vec_t           = vec16<T>;
  constexpr auto traits = vec_traits_t<vec_t> {}; // matching SSE
  assert (0 == ((std::uintptr_t) in.data() % 16)); // assert alignment

  std::array<T, 2> inc;
  for (uint i = 0; i < inc.size(); ++i) {
    inc[i] = ((T) traits.size) * (gain_end[i] - gain_beg[i]) / (T) (in_size);
  }
  auto gain = gain_beg;

  auto simd_ramp = [=, &gain] (std::array<T*, 2> c, uint blocks) {
    T* end = c[0] + (blocks * traits.size);

    while (c[0] < end) {
      auto v1 = vec_load<vec_t> (c[0]);
      auto v2 = vec_load<vec_t> (c[1]);
      v1 *= gain[0];
      v2 *= gain[1];
      v1 += v2;
      vec_store (c[0], v1);
      gain[0] += inc[0];
      gain[1] += inc[1];
      c[0] += traits.size;
      c[1] += traits.size;
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

  block_divide (traits.bytes, in, in_size, simd_ramp, regular_ramp);
}
//------------------------------------------------------------------------------
} // namespace artv
