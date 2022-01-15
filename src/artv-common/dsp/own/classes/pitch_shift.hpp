#pragma once

#include <cmath>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/misc/interpolation.hpp"

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
  void reset (crange<V> mem) // has a latency of "mem.size() / 2"
  {
    _mem.reset (mem);
    _size_recip = 1. / (float) mem.size();
    _shift      = 0;
    _rint       = 0;
    _rfrac      = 0.;
  }
  //----------------------------------------------------------------------------
  void set (float amt_semitones)
  {
    _shift = exp (amt_semitones * (1. / 12.) * M_LN2);
    _rint  = _mem.abs_pos();
    _rfrac = 0.;
  }
  //----------------------------------------------------------------------------
  V tick (V in)
  {
    using T = vec_value_type_t<V>;
    _mem.push (in);
    uint wpos = _mem.abs_pos();

    T rpos = (Free_running) ? (T) _rint + _rfrac + _shift : (T) wpos * _shift;
    _rint  = (uint) rpos;
    _rfrac = rpos - (T) _rint;
    // avoid accumulating precission loss by not letting "_rint" to become big,
    // this is for the free running case.
    _rint = _mem.wrap_abs_pos (_rint);

    V p1  = _mem.get_abs (_rint);
    V p2  = _mem.get_abs (_rint + 1);
    V ps1 = linear_interp::get ({p1, p2}, vec_set<V> (_rfrac));

    p1    = _mem.get_abs (_rint + _mem.size() / 2);
    p2    = _mem.get_abs (_rint + 1 + _mem.size() / 2);
    V ps2 = linear_interp::get ({p1, p2}, vec_set<V> (_rfrac));

    uint del       = _rint - wpos;
    T    del1_norm = (T) _mem.wrap_abs_pos (del) * _size_recip;
    del += _mem.size() / 2;
    T del2_norm = (T) _mem.wrap_abs_pos (del) * _size_recip;
    // sin window
    return vec_set<V> (sin (del1_norm * M_PI)) * ps1
      + vec_set<V> (sin (del2_norm * M_PI)) * ps2;
  }
  //----------------------------------------------------------------------------
private:
  pow2_circular_buffer<V> _mem;
  float                   _shift;
  float                   _size_recip;
  uint                    _rint;
  float                   _rfrac;
};
//------------------------------------------------------------------------------
} // namespace artv
