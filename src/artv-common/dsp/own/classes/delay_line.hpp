#pragma once

#include <array>
#include <cstring>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/misc/interpolation.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <class V>
class pow2_circular_buffer {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem)
  {
    assert (is_pow2 (mem.size()));
    _pos  = 0;
    _mem  = mem.data();
    _mask = mem.size() - 1;
  };
  //----------------------------------------------------------------------------
  V get (uint idx) const { return get_abs (_pos - idx); }
  V operator[] (uint idx) const { return get (idx); }
  //----------------------------------------------------------------------------
  void push (const V& value)
  {
    _mem[_pos & _mask] = value;
    ++_pos;
  }
  //----------------------------------------------------------------------------
  uint size() const { return _mask + 1; }
  //----------------------------------------------------------------------------
  // For pitch shifters, getting a position that can have a factor multiplied
  // and then accessed, so instead of returning "pos" from "0" to "size -1" it
  // is returned from "size" to "(2 * size) -1"
  uint abs_pos() const { return wrap_abs_pos (_pos) + 1 + _mask; }
  //----------------------------------------------------------------------------
  // For pitch shifters, accessing by the absolute position of the buffer.
  // Normally accessing a position relative to the write head is what is
  // desired, those are "get" and "operator[]".
  V get_abs (uint abs_pos) const { return _mem[wrap_abs_pos (abs_pos)]; }
  //----------------------------------------------------------------------------
  // also for pitch shifters
  uint wrap_abs_pos (uint abs_pos) const { return abs_pos & _mask; }
  //----------------------------------------------------------------------------
private:
  V*   _mem;
  uint _mask;
  uint _pos;
};
//------------------------------------------------------------------------------
// just a single allocation multichannel delay line with random access.
// Capacity is always rounded up to a power of two to avoid divisions when
// wrapping around.
template <class T>
class delay_line {
public:
  //----------------------------------------------------------------------------
  void reset (uint channels, uint max_delay)
  {
    reset_with_extra_tail_samples (channels, max_delay, 0);
  }
  //----------------------------------------------------------------------------
  void get (crange<T> channels_out, uint idx) // idx 0 is previously pushed
  {
    assert (channels_out.size() >= _channels);
    assert (idx < max_delay_samples());

    uint pos = ((_wpos + idx) & _capacity_mask) * _channels;
    memcpy (channels_out.data(), &_samples[pos], sizeof (T) * _channels);
  }
  //----------------------------------------------------------------------------
  T get (uint channel, uint idx)
  {
    assert (channel < _channels);
    assert (idx < max_delay_samples());

    return _samples[(((_wpos + idx) & _capacity_mask) * _channels) + channel];
  }
  //----------------------------------------------------------------------------
  void push (crange<const T> channels_in)
  {
    assert (channels_in.size() >= _channels);
    --_wpos;
    uint pos = (_wpos & _capacity_mask) * _channels;
    memcpy (&_samples[pos], channels_in.data(), sizeof (T) * _channels);
  }
  //----------------------------------------------------------------------------
  uint max_delay_samples() const { return _capacity_mask + 1; }
  //----------------------------------------------------------------------------
  uint n_channels() const { return _channels; }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  T* reset_with_extra_tail_samples (
    uint channels,
    uint max_delay,
    uint tail_samples)
  {
    max_delay = pow2_round_ceil (max_delay);
    _samples.clear();
    _samples.resize ((max_delay + tail_samples) * channels); // 0 fill
    _wpos          = 0;
    _capacity_mask = max_delay - 1;
    _channels      = channels;
    return &_samples[_samples.size() - channels];
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  std::vector<T> _samples;
  uint           _wpos;
  uint           _capacity_mask;
  uint           _channels;
};
//------------------------------------------------------------------------------
// interpolated and modulated delay line, requires an stateless interpolator
// that supports random access. The modulator is left external, so it can be
// shared by different delay lines.
template <
  class V,
  class Interp                       = linear_interp,
  enable_if_vec_of_float_point_t<V>* = nullptr>
class interp_mod_delay_line {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem)
  {
    assert (mem.size() > interp_margin);
    _mem.reset (mem);
    _max_delay = (float) (_mem.size() - 1 - interp_margin);
  }
  //----------------------------------------------------------------------------
  void set (float time, float mod_freq, float mod_depth, float samplerate)
  {
    set_time (time);
    set_mod_depth (mod_depth);
  }
  //----------------------------------------------------------------------------
  void set_time (float samples)
  {
    _delay = std::min ((float) _max_delay, samples);
  }
  //----------------------------------------------------------------------------
  // the delay time will oscillate on an excursion of ± "samples" centered
  // around the delay time.
  void set_mod_depth (float samples) { _depth = samples; }
  //----------------------------------------------------------------------------
  V tick (V in, vec_value_type_t<V> mod_in = (vec_value_type_t<V>) 0)
  {
    assert (mod_in >= (vec_value_type_t<V>) -1);
    assert (mod_in <= (vec_value_type_t<V>) 1);

    using T = vec_value_type_t<V>;

    float fpdel = _delay + mod_in * _depth;
    fpdel       = std::clamp (fpdel, 0.f, _max_delay);
    auto del    = (uint) fpdel;
    auto frac   = fpdel - (float) del;

    std::array<V, Interp::n_points> samples;
    for (uint i = 0; i < samples.size(); ++i) {
      samples[i] = _mem[del + i];
    }
    V ret = Interp::get (samples, vec_set<V> (frac));
    _mem.push (in);
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  static constexpr uint interp_margin = Interp::n_points - 1;

  pow2_circular_buffer<V> _mem;
  float                   _delay;
  float                   _depth;
  float                   _max_delay;
};
//------------------------------------------------------------------------------
// Basic reverb building block. Strictly not a delay line.
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
class schroeder_allpass {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem) { _mem.reset (mem); }
  //----------------------------------------------------------------------------
  void set (uint delay_samples, V gain)
  {
    set_time (delay_samples);
    set_gain (gain);
  }
  //----------------------------------------------------------------------------
  void set_time (uint samples)
  {
    assert (samples < _mem.size());
    _delay = samples;
  }
  //----------------------------------------------------------------------------
  void set_gain (V gain) { _gain = gain; }
  //----------------------------------------------------------------------------
  V tick (V in)
  {
    V y1 = _mem[_delay];
    V y  = in + y1 * -_gain;
    _mem.push (y);
    return y * _gain + y1;
  }
  //----------------------------------------------------------------------------
private:
  pow2_circular_buffer<V> _mem;
  uint                    _delay {};
  V                       _gain {};
  V                       _z1 {};
};
//------------------------------------------------------------------------------
// A delay line that can be fractional but not modulated fast. Untested, wasn't
// needed after all.
template <class T>
class allpass_interpolated_delay_line : public delay_line<T> {
public:
  //----------------------------------------------------------------------------
  // shadowing.
  //----------------------------------------------------------------------------
  void reset (uint channels, uint capacity_log2)
  {
    _allpass_states
      = this->reset_with_extra_tail_samples (channels, capacity_log2, 1);
  }
  //----------------------------------------------------------------------------
  T get (uint channel, float idx)
  {
    uint  i_idx = (float) idx;
    float fract = idx - i_idx;
    return get (channel, i_idx, get_fractional_coeff (fract));
  }
  //----------------------------------------------------------------------------
  constexpr double get_fractional_coeff (double fraction) // 0 to 1
  {
    return (1.0 - fraction) / (1.0 + fraction);
  }
  //----------------------------------------------------------------------------
  // just to be able to save divisions in cases the fraction is constant. See
  // "get (uint channel, float idx)"
  T get (uint channel, uint idx, double fractional_coeff)
  {
    assert (channel < this->n_channels());
    assert (idx < (this->max_delay_samples() - 1));

    std::array<T, 2> points = {get (channel, idx), get (channel, idx + 1)};
    // this is built based on 1 push 1 retrieval
    return warped_allpass_interpolate (
      fractional_coeff, _allpass_states[channel], points);
  }
  //----------------------------------------------------------------------------
private:
  // magic from RS-MET's library. Seems to be
  // https://ccrma.stanford.edu/~jos/pasp/First_Order_Allpass_Interpolation.html
  double warped_allpass_interpolate (
    double           coeff,
    float&           prev,
    std::array<T, 2> points)
  {
    double ret = points[0] + (coeff * points[1]) - (0.999 * coeff * prev);
    prev       = ret;
    return ret;
  }
  //----------------------------------------------------------------------------
  T* _allpass_states;
};
//------------------------------------------------------------------------------
} // namespace artv
