#pragma once

#include <array>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/misc/interpolation.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint N, class V, enable_if_vec_of_float_point_t<V>* = nullptr>
class householder_matrix {
public:
  //----------------------------------------------------------------------------
  static constexpr uint size = N;
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) // mem.size() / N == a power of 2
  {
    auto ptr   = mem.data();
    uint msize = mem.size() / N;

    assert (msize && ((msize * N) == mem.size()));

    for (uint i = 0; i < N; ++i, ptr += msize) {
      _del[i].reset (make_crange (ptr, msize));
    }
  }
  //----------------------------------------------------------------------------
  void set (std::array<uint, N> t, V feedback)
  {
    set_times (t);
    set_feedback (feedback);
  }
  //----------------------------------------------------------------------------
  void set_times (std::array<uint, N> t)
  {
    for (uint v : t) {
      assert (v < _del[0].size());
    }
    _t = t;
  }
  //----------------------------------------------------------------------------
  void set_feedback (V feedback) { _fb = feedback; }
  //----------------------------------------------------------------------------
  V tick (V in)
  {
    using T = vec_value_type_t<V>;

    std::array<V, N> v;
    V                factor = vec_set<V> ((T) 0);

    for (uint i = 0; i < N; ++i) {
      v[i] = _del[i].get (_t[i]);
      factor += v[i];
    }
    factor *= vec_set<V> ((T) -0.5);
    for (uint i = 0; i < N; ++i) {
      _del[i].push (in + (v[i] + factor) * _fb);
    }
    return factor;
  }
  //----------------------------------------------------------------------------
private:
  std::array<pow2_circular_buffer<V>, N> _del {};
  std::array<uint, N>                    _t {};
  V                                      _fb {};
};
//------------------------------------------------------------------------------<
} // namespace artv
