#pragma once

#include <cmath>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
// Naive but CPU-light pitch shifter, found on Shabronic's FDN Reverb JSFX.
// Based on "sin" crosfading between two buffers. It has a latency of half the
// buffer size. Only for buffers with power of 2 sizes.
//------------------------------------------------------------------------------
template <class V, bool Free_running, enable_if_floatpt_vec_t<V>* = nullptr>
class pitch_shift_sin {
public:
  //----------------------------------------------------------------------------
  struct reader {
    fixpt<0, 20, 22, 0> pos {}; // u64
    fixpt<0, 10, 22, 0> delta {}; // u32
  };
  //----------------------------------------------------------------------------
  void reset (xspan<V> mem) // has a latency of "mem.size() / 2"
  {
    assert (mem.size() <= (decltype (reader::pos)::max_int()) + 1);

    _mem.reset (mem);
    _size_recip = 1. / (float) mem.size();
  }
  //----------------------------------------------------------------------------
  // many shift factors can be read from the same buffer.
  void set_reader (reader& r, float amt_semitones, bool resync = true)
  {
    r.delta.load_float (exp (amt_semitones * (1. / 12.) * M_LN2));
    if (resync) {
      r.pos.load_int (_mem.abs_pos());
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
      tk.pos.load_int (wpos);
      tk.pos *= tk.delta;
    }

    auto rint  = tk.pos.as_int();
    auto rfrac = tk.pos.float_fractional();

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
