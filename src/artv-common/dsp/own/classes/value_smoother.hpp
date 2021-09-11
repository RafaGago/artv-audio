#pragma once

#include <cassert>
#include <cstring>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
template <class T, uint N, uint Simd_size = sse_bytes>
class value_smoother {
public:
  //----------------------------------------------------------------------------
  using value_type                = T;
  static constexpr uint simd_size = Simd_size;
  static constexpr uint size      = N;
  using vec_type                  = vec<T, simd_size / sizeof (T)>;
  //----------------------------------------------------------------------------
  void reset (T samplerate)
  {
    memset (&_target, 0, sizeof _target);
    memset (&_current, 0, sizeof _current);
    onepole_smoother::reset_coeffs (
      make_crange (_coeff), make_vec_x1<T> (1. / 0.1), samplerate);
  }
  //----------------------------------------------------------------------------
  void set_to_target() { _current = _target; }
  //----------------------------------------------------------------------------
  void set (T v, uint idx)
  {
    assert (idx < N);
    _target[idx] = v;
  }
  //----------------------------------------------------------------------------
  T get (uint idx) const
  {
    assert (idx < N);
    return _current[idx];
  }
  //----------------------------------------------------------------------------
  void tick (uint samples = 1)
  {
    static constexpr uint step = vec_traits_t<vec_type>::size;
    for (uint i = 0; i < samples; ++i) {
      for (uint j = 0; j < _target.size(); j += step) {
        vec_type out = onepole_smoother::tick (
          make_crange (_coeff),
          make_crange (&_current[j], step),
          vec_load<vec_type> (&_target[j]),
          single_coeff_set_tag {});
        vec_store (&_current[j], out);
      }
    }
  }
  //----------------------------------------------------------------------------
  // for when smoothing both coefficients and values
  auto const& get_all() const { return _current; }
  //----------------------------------------------------------------------------
private:
  alignas (sse_bytes) simd_array<T, N, sse_bytes> _target;
  alignas (sse_bytes) simd_array<T, N, sse_bytes> _current;
  T _coeff;
};
//------------------------------------------------------------------------------
} // namespace artv
