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
// - non-owned memory
// - single channel
// - random access.
// - power of 2 capacity.
// - not interleaved
// - highest indexes = older samples
template <class V>
class pow2_circular_buffer {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem)
  {
    assert (is_pow2 (mem.size()));
    crange_memset (mem, 0);
    _pos  = 0;
    _mem  = mem.data();
    _mask = mem.size() - 1;
  };
  //----------------------------------------------------------------------------
  // idx from newest to oldest, so idx=0 -> z, idx=1 -> z-1, etc.
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
// - not managed memory
// - multichannel
// - random access.
// - power of 2 capacity.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
//
// Note: It doesn't reuse "pow2_circular_buffer" to be smaller, as "mask" and
// "pos" is not required for each channel.
//------------------------------------------------------------------------------

namespace detail {
//------------------------------------------------------------------------------
template <class T>
struct static_delay_line_base {
public:
  //----------------------------------------------------------------------------
  using value_type = T;
  //----------------------------------------------------------------------------
  constexpr void reset (crange<T> mem, uint n_channels)
  {
    uint size = mem.size() / n_channels;
    assert ((mem.size() % n_channels) == 0);
    assert (size != 0 && is_pow2 (size));

    crange_memset (mem, 0);
    _z          = mem.data();
    _mask       = size - 1;
    _pos        = 0;
    _n_channels = n_channels;
  }
  //----------------------------------------------------------------------------
  constexpr uint size() const { return _mask + 1; }
  //----------------------------------------------------------------------------
  constexpr uint n_channels() const { return _n_channels; }
  //----------------------------------------------------------------------------
  T*   _z;
  uint _mask;
  uint _pos;
  uint _n_channels;
};
//------------------------------------------------------------------------------
} // namespace detail

template <class T, bool Interleaved>
class static_delay_line;

// interleaved overload
template <class T>
class static_delay_line<T, true> : private detail::static_delay_line_base<T> {
private:
  using base = detail::static_delay_line_base<T>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr value_type get (time_type sample, uint chnl) const
  {
    auto& t = *this;

    assert (sample <= t._mask && "Unintended wraparound?");
    assert (chnl < t._n_channels);

    uint pos = (t._pos - sample) & t._mask;
    return t._z[pos * t._n_channels + chnl];
  }
  //----------------------------------------------------------------------------
  constexpr void push (const crange<value_type> row)
  {
    auto& t = *this;

    assert (row.size() >= t._n_channels);
    ++t._pos;
    uint pos = t._pos & t._mask;
    memcpy (&t._z[pos * t._n_channels], row.data(), sizeof (T) * t._n_channels);
  }
  //----------------------------------------------------------------------------
};

// non-interleaved overload
template <class T>
class static_delay_line<T, false> : private detail::static_delay_line_base<T> {
private:
  using base = detail::static_delay_line_base<T>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr T get (time_type sample, uint chnl) const
  {
    auto& t = *this;

    assert (sample <= t._mask && "Unintended wraparound?");
    assert (chnl < t._n_channels);

    uint pos   = (t._pos - sample) & t._mask;
    uint chmem = (t._mask + 1) * chnl;
    return t._z[chmem + pos];
  }
  //----------------------------------------------------------------------------
  constexpr void push (const crange<value_type> row)
  {
    auto& t = *this;

    assert (row.size() >= t._n_channels);
    ++t._pos;
    uint pos = t._pos & t._mask;
    for (uint i = 0, chmem = 0; i < t._n_channels; ++i, chmem += t._mask + 1) {
      t._z[chmem + pos] = row[i];
    }
  }
  //----------------------------------------------------------------------------
};
namespace detail {
template <class Delay_line_base, class Interpolation>
class interpolated_delay_line : private Delay_line_base {
private:
  using base = Delay_line_base;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::push;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using interp     = Interpolation;
  using time_type  = std::conditional_t<interp::n_points == 1, uint, float>;
  //----------------------------------------------------------------------------
  // integer index access
  value_type get_raw (uint sample, uint channel)
  {
    return Delay_line_base::get (sample, channel);
  }
  //----------------------------------------------------------------------------
  // interpolated overload
  value_type get (time_type sample, uint channel)
  {
    if constexpr (interp::n_points == 1) {
      return get_raw (sample, channel);
    }
    else {
      std::array<value_type, interp::n_points> y;

      auto sample_uint = (uint) sample;
      auto frac        = sample - (time_type) sample_uint;

      for (uint i = 0; i < interp::n_points; ++i) {
        y[i] = get_raw (sample_uint + i, channel);
      }
      if constexpr (is_vec_v<value_type>) {
        return interp::get (y, vec_set<value_type> (frac));
      }
      else {
        return interp::get (y, frac);
      }
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

} // namespace detail

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - random access.
// - power of 2 capacity.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
//------------------------------------------------------------------------------
template <class T, class Interp = linear_interp, bool Interleaved = false>
using interpolated_delay_line
  = detail::interpolated_delay_line<static_delay_line<T, Interleaved>, Interp>;

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - random access.
// - power of 2 capacity.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
// - convenience functions for taking modulation
//------------------------------------------------------------------------------
template <class T, class Interp = linear_interp, bool Interleaved = false>
class modulable_delay_line
  : private interpolated_delay_line<T, Interp, Interleaved> {
private:
  using base = interpolated_delay_line<T, Interp, Interleaved>;

public:
  //----------------------------------------------------------------------------
  using base::get;
  using base::get_raw;
  using base::n_channels;
  using base::push;
  using base::size;
  using time_type  = typename base::time_type;
  using interp     = typename base::interp;
  using value_type = typename base::value_type;
  //----------------------------------------------------------------------------
  template <class... Ts>
  constexpr void reset (Ts&&... args)
  {
    base::reset (std::forward<Ts> (args)...);
    _max_delay = size() - 1 - interp::n_points - 1;
  }
  //----------------------------------------------------------------------------
  // interpolated overload with modulation, just for convenience
  value_type get (time_type sample, float mod_amt, float max_mod, uint channel)
  {
    // in case we are using a non interpolated delay line, no-op otherwise
    auto z = (float) sample;
    z      = z + mod_amt * max_mod;
    z      = std::clamp (z, 0.f, _max_delay);
    return get ((time_type) z, channel);
  }
  //----------------------------------------------------------------------------
private:
  float _max_delay;
};
//------------------------------------------------------------------------------
// - single dynamic allocation
// - multichannel
// - random access.
// - power of 2 capacity.
// - not interleaved
// - same sized channels
// - highest indexes = older samples
template <class T, bool Interleaved>
class dynamic_delay_line {
public:
  //----------------------------------------------------------------------------
  using value_type = T;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  ~dynamic_delay_line() { do_delete(); }
  //----------------------------------------------------------------------------
  void reset (uint n_channels, time_type max_delay)
  {
    max_delay = pow2_round_ceil (max_delay);
    do_delete();
    uint elems = max_delay * n_channels;
    _mem       = new value_type[elems];
    _z.reset (make_crange (_mem, elems), n_channels);
  }
  //----------------------------------------------------------------------------
  T get (time_type sample, uint channel) { return _z.get (sample, channel); }
  //----------------------------------------------------------------------------
  void push (const crange<T> x) { _z.push (x); }
  //----------------------------------------------------------------------------
  uint size() const { return _z.size(); }
  //----------------------------------------------------------------------------
  uint n_channels() const { return _z.n_channels(); }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void do_delete()
  {
    if (_mem) {
      delete[] _mem;
      _mem = nullptr;
    }
  }
  //----------------------------------------------------------------------------
  static_delay_line<T, Interleaved> _z;
  T*                                _mem = nullptr;
};
//------------------------------------------------------------------------------
// Basic reverb building block. Strictly not a delay line. Not sure if it
// belongs here. It probably needs to be moved.
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
} // namespace artv
