#pragma once

#include <cassert>
#include <cstring>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {
//------------------------------------------------------------------------------
// bulk SIMD smoothing of many parameters of the same floating point type via a
// 1 pole lowpass. Values accessed by index.
template <class T, uint N, uint Simd_size = sse_bytes>
class value_smoother {
public:
  //----------------------------------------------------------------------------
  using value_type                = T;
  static constexpr uint simd_size = Simd_size;
  static constexpr uint size      = N;
  using vec_type                  = vec<T, simd_size / sizeof (T)>;
  static constexpr uint vec_size  = vec_traits_t<vec_type>::size;
  //----------------------------------------------------------------------------
  void reset (T t_spl)
  {
    memset (&_target, 0, sizeof _target);
    memset (&_current, 0, sizeof _current);
    using x1_t = vec<T, 1>;
    onepole_smoother::reset_coeffs (
      make_crange (_coeff).cast (x1_t {}), vec_set<x1_t> (1. / 0.1), t_spl);
  }
  //----------------------------------------------------------------------------
  void set_to_target() { _current = _target; }
  //----------------------------------------------------------------------------
  void set (T v, uint idx)
  {
    assert (idx < N);
    _target[idx / vec_size][idx % vec_size] = v; // vec_size is a power of 2
  }
  //----------------------------------------------------------------------------
  T get (uint idx) const
  {
    assert (idx < N);
    return _current[idx / vec_size][idx % vec_size]; // vec_size is a power of 2
  }
  //----------------------------------------------------------------------------
  void tick (uint samples = 1)
  {
    for (uint j = 0; j < _target.size(); ++j) {
      for (uint i = 0; i < samples; ++i) {
        _current[j] = onepole_smoother::tick<vec_type> (
          make_crange (_coeff), make_crange (_current[j]), _target[j]);
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  alignas (sse_bytes) simd_vec_array<T, N, sse_bytes> _target;
  alignas (sse_bytes) simd_vec_array<T, N, sse_bytes> _current;
  T _coeff;
};
//------------------------------------------------------------------------------
} // namespace artv
