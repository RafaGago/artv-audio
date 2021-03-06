#pragma once

#include <cmath>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
// Naive but CPU-light pitch shifter, found on Shabronic's FDN Reverb JSFX.
// Based on "sin" crosfading between two buffers. It has a latency of half the
// buffer size. Only for buffers with power of 2 sizes.
//------------------------------------------------------------------------------
template <
  class V,
  bool Free_running,
  enable_if_vec_of_float_point_t<V>* = nullptr>
class pitch_shift_sin {
public:
  //----------------------------------------------------------------------------
  struct reader {
    fixed_point<20, 22, false> pos {}; // u64
    fixed_point<10, 22, false> delta {}; // u32
  };
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) // has a latency of "mem.size() / 2"
  {
    assert (mem.size() <= (decltype (reader::pos)::int_max) + 1);

    _mem.reset (mem);
    _size_recip = 1. / (float) mem.size();
  }
  //----------------------------------------------------------------------------
  // many shift factors can be read from the same buffer.
  void set_reader (reader& r, float amt_semitones, bool resync = true)
  {
    r.delta = exp (amt_semitones * (1. / 12.) * M_LN2);
    if (resync) {
      r.pos = _mem.abs_pos();
    }
  }
  //----------------------------------------------------------------------------
  V read (reader& tk)
  {
    using T = vec_value_type_t<V>;

    uint wpos = _mem.abs_pos();

    if (Free_running) {
      tk.pos += tk.delta;
    }
    else {
      tk.pos = wpos;
      tk.pos.mul (tk.delta);
    }

    auto rint  = tk.pos.integer();
    auto rfrac = tk.pos.fraction();

    V p1  = _mem.get_abs (rint);
    V p2  = _mem.get_abs (rint + 1);
    V ps1 = linear_interp::tick ({p1, p2}, vec_set<V> (rfrac));

    p1    = _mem.get_abs (rint + _mem.size() / 2);
    p2    = _mem.get_abs (rint + 1 + _mem.size() / 2);
    V ps2 = linear_interp::tick ({p1, p2}, vec_set<V> (rfrac));

    // sin equal-power crossfade
    uint del     = rint - wpos;
    T    crossf1 = (T) _mem.wrap_abs_pos (del) * _size_recip * M_PI;
    del += _mem.size() / 2;
    T crossf2 = (T) _mem.wrap_abs_pos (del) * _size_recip * M_PI;
    return sin (crossf1) * ps1 + sin (crossf2) * ps2;
  }
  //----------------------------------------------------------------------------
  void push (V in) { _mem.push (in); }
  //----------------------------------------------------------------------------
private:
  pow2_circular_buffer<V> _mem;
  float                   _size_recip;
};
//------------------------------------------------------------------------------
} // namespace artv
