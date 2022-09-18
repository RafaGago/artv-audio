#pragma once

#include <cassert>
#include <limits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/classes/ducker.hpp"
namespace artv {
// ducker, for e.g. usage on reverbs/echoes
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
class add_ducker {
public:
  //----------------------------------------------------------------------------
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v)
  {
    if (v == _ducking_threshold) {
      return;
    }
    _ducking_threshold = v;
    _ducker.set_threshold (vec_set<V> (v));
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -40.f, 12.f, 12.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    if (v == _ducking_speed) {
      return;
    }
    _ducking_speed = v;
    v *= 0.01f;
    _ducker.set_speed (vec_set<V> (v * v), (vec_value_type_t<V>) _t_spl);
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  template <class T, class F>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples, F child_f)
  {
    assert (outs.size() >= 2);
    assert (ins.size() >= 2);

    using VT = vec_value_type_t<V>;

    constexpr uint                block_max_size = 32;
    std::array<V, block_max_size> ducker_gain;
    std::array<T*, 2>             out {{outs[0], outs[1]}};
    std::array<T const*, 2>       in {{ins[0], ins[1]}};

    for (uint b = 0; b < samples; b += block_max_size) {
      uint blocksize = std::min (samples - b, block_max_size);

      for (uint i = 0; i < blocksize; ++i) {
        ducker_gain[i] = _ducker.tick (V {(VT) in[0][i], (VT) in[1][i]});
      }

      child_f (xspan {out}, xspan {in}, blocksize);

      for (uint i = 0; i < blocksize; ++i) {
        out[0][i] *= ducker_gain[i][0];
        out[1][i] *= ducker_gain[i][1];
      }

      out[0] += block_max_size;
      out[1] += block_max_size;
      in[0] += block_max_size;
      in[1] += block_max_size;
    }
  }
  //----------------------------------------------------------------------------
  void reset (float samplerate)
  {
    _t_spl             = 1. / samplerate;
    _ducking_threshold = 2000.f;
    _ducking_speed     = 2000.f;
    set (
      ducking_threshold_tag {},
      get_parameter (ducking_threshold_tag {}).defaultv);
    set (ducking_speed_tag {}, get_parameter (ducking_speed_tag {}).defaultv);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  ducker<V> _ducker;
  float     _ducking_threshold {};
  float     _ducking_speed {};
  float     _t_spl {};
};
//------------------------------------------------------------------------------
} // namespace artv
