#pragma once

#include <array>
#include <cstring>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"

#include "artv-common/dsp/own/parts/interpolation/stateful.hpp"
#include "artv-common/dsp/own/parts/interpolation/stateless.hpp"

namespace artv {

//------------------------------------------------------------------------------
// - non-owned memory
// - single channel
// - random access.
// - power of 2 capacity.
// - not interleaved
// - highest indexes = older samples (LIFO)
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
    ++_pos;
    _mem[_pos & _mask] = value;
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
// - non-owned memory
// - single channel
// - random access.
// - not interleaved
// - highest indexes = older samples (LIFO)
template <class V>
class circular_buffer {
public:
  //----------------------------------------------------------------------------
  void reset (crange<V> mem)
  {
    assert (mem.size());
    crange_memset (mem, 0);
    _pos  = 0;
    _mem  = mem.data();
    _size = mem.size();
  };
  //----------------------------------------------------------------------------
  // idx from newest to oldest, so idx=0 -> z, idx=1 -> z-1, etc.
  V get (uint idx) const
  {
    assert (idx < _size);
    idx = _pos - idx;
    idx += idx > _size ? _size : 0;
    return _mem[idx];
  }
  //----------------------------------------------------------------------------
  V operator[] (uint idx) const { return get (idx); }
  //----------------------------------------------------------------------------
  void push (const V& value)
  {
    ++_pos;
    _pos       = _pos < _size ? _pos : 0;
    _mem[_pos] = value;
  }
  //----------------------------------------------------------------------------
  uint size() const { return _size; }
  //----------------------------------------------------------------------------
private:
  V*   _mem;
  uint _pos;
  uint _size;
};
//------------------------------------------------------------------------------
namespace detail {

//------------------------------------------------------------------------------
template <class T, bool Use_pow2_sizes_sizes>
struct static_delay_line_base;

//------------------------------------------------------------------------------
template <class T>
struct static_delay_line_base<T, true> {
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
template <class T>
struct static_delay_line_base<T, false> {
public:
  //----------------------------------------------------------------------------
  using value_type = T;
  //----------------------------------------------------------------------------
  constexpr void reset (crange<T> mem, uint n_channels)
  {
    uint size = mem.size() / n_channels;
    assert ((mem.size() % n_channels) == 0);
    assert (size != 0);

    crange_memset (mem, 0);
    _z          = mem.data();
    _size       = size;
    _pos        = 0;
    _n_channels = n_channels;
  }
  //----------------------------------------------------------------------------
  constexpr uint size() const { return _size; }
  //----------------------------------------------------------------------------
  constexpr uint n_channels() const { return _n_channels; }
  //----------------------------------------------------------------------------
  T*   _z;
  uint _size;
  uint _pos;
  uint _n_channels;
};

//------------------------------------------------------------------------------
} // namespace detail

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - random access.
// - Optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
//
// Note: It doesn't reuse "pow2_circular_buffer" to be smaller, as "mask" and
// "pos" is not required for each channel.
//------------------------------------------------------------------------------
template <class T, bool Interleaved, bool Use_pow2_sizes_sizes>
class static_delay_line;

// interleaved and pow2 overload
template <class T>
class static_delay_line<T, true, true>
  : private detail::static_delay_line_base<T, true> {
private:
  using base = detail::static_delay_line_base<T, true>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr value_type get (time_type delay_spls, uint chnl) const
  {
    auto& t = *this;

    assert (delay_spls <= t._mask && "Unintended wraparound?");
    assert (chnl < t._n_channels);

    uint pos = (t._pos - delay_spls) & t._mask;
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

// interleaved non pow2 overload
template <class T>
class static_delay_line<T, true, false>
  : private detail::static_delay_line_base<T, false> {
private:
  using base = detail::static_delay_line_base<T, false>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr value_type get (time_type delay_spls, uint chnl) const
  {
    auto& t = *this;

    assert (delay_spls <= t._size && "Out of range");
    assert (chnl < t._n_channels);

    uint pos = t._pos - delay_spls;
    pos += pos > t._size ? t._size : 0;

    return t._z[pos * t._n_channels + chnl];
  }
  //----------------------------------------------------------------------------
  constexpr void push (const crange<value_type> row)
  {
    auto& t = *this;

    assert (row.size() >= t._n_channels);
    ++t._pos;
    t._pos = t._pos < t._size ? t._pos : 0;
    memcpy (
      &t._z[t._pos * t._n_channels], row.data(), sizeof (T) * t._n_channels);
  }
  //----------------------------------------------------------------------------
};

// non-interleaved pow 2 overload
template <class T>
class static_delay_line<T, false, true>
  : private detail::static_delay_line_base<T, true> {
private:
  using base = detail::static_delay_line_base<T, true>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr T get (time_type delay_spls, uint chnl) const
  {
    auto& t = *this;

    assert (delay_spls <= t._mask && "Unintended wraparound?");
    assert (chnl < t._n_channels);

    uint pos   = (t._pos - delay_spls) & t._mask;
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

// non-interleaved non-pow 2 overload
template <class T>
class static_delay_line<T, false, false>
  : private detail::static_delay_line_base<T, false> {
private:
  using base = detail::static_delay_line_base<T, false>;

public:
  //----------------------------------------------------------------------------
  using base::n_channels;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  using time_type  = uint;
  //----------------------------------------------------------------------------
  constexpr T get (time_type delay_spls, uint chnl) const
  {
    auto& t = *this;

    assert (delay_spls <= t._size && "Out of range");
    assert (chnl < t._n_channels);

    uint pos = t._pos - delay_spls;
    pos += pos > t._size ? t._size : 0;
    uint chmem = t._size * chnl;
    return t._z[chmem + pos];
  }
  //----------------------------------------------------------------------------
  constexpr void push (const crange<value_type> row)
  {
    auto& t = *this;

    assert (row.size() >= t._n_channels);
    ++t._pos;
    t._pos = t._pos < t._size ? t._pos : 0;

    for (uint i = 0, chmem = 0; i < t._n_channels; ++i, chmem += t._size) {
      t._z[chmem + t._pos] = row[i];
    }
  }
  //----------------------------------------------------------------------------
};

namespace detail {
//------------------------------------------------------------------------------
// "Delay_line_base" = the underlying buffer type, e.g. "static_delay_line".
template <class Delay_line_base>
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
  //----------------------------------------------------------------------------
  // integer index access
  value_type get_raw (uint sample, uint channel)
  {
    return base::get (sample, channel);
  }
  //----------------------------------------------------------------------------
  template <uint N_points, uint X_offset, class InterpFunctor>
  value_type get (
    uint            delay_spls_int,
    float           delay_spls_frac,
    uint            channel,
    InterpFunctor&& interp)
  {
    std::array<value_type, N_points> y;
    // probably not making it an assert but a clamp, let's see?
    assert (
      std::clamp<float> (
        delay_spls_frac + (float) delay_spls_int,
        X_offset,
        size() - N_points + X_offset)
      == (delay_spls_frac + (float) delay_spls_int));

    // e.g. Catmull-Rom interpolates between the 2 central points of 4.
    delay_spls_int -= X_offset;

    for (uint i = 0; i < y.size(); ++i) {
      y[i] = get_raw (delay_spls_int + i, channel);
    }
    return interp (y, make_vec<value_type> (delay_spls_frac));
  }
  //----------------------------------------------------------------------------
  template <uint N_points, uint X_offset, class InterpFunctor>
  value_type get (float delay_spls, uint channel, InterpFunctor&& interp)
  {
    auto spls = (uint) delay_spls;
    auto frac = delay_spls - (float) spls;
    return get<N_points, X_offset> (spls, frac, channel, std::move (interp));
  }
  //----------------------------------------------------------------------------
  // Stateless_interp = a class on
  // "artv-common/dsp/own/parts/interpolation/stateless.hpp"
  template <class Stateless_interp>
  value_type get (float delay_spls, uint channel)
  {
    using interp = Stateless_interp;
    return get<interp::n_points, interp::x_offset> (
      delay_spls, channel, [] (auto y, auto x) { return interp::tick (y, x); });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// "Delay_line_base" = the underlying buffer type, e.g. "static_delay_line".
// "Interpolation" one of the classes on
// "artv-common/dsp/own/parts/interpolation/stateful.hpp"
template <class Delay_line_base, class Interpolation>
class statefully_interpolated_delay_line
  : private interpolated_delay_line<Delay_line_base> {
private:
  using base = interpolated_delay_line<Delay_line_base>;

public:
  //----------------------------------------------------------------------------
  using base::get_raw;
  using base::n_channels;
  using base::push;
  using base::size;
  using value_type = typename base::value_type;
  using interp     = Interpolation;
  static constexpr uint n_interp_states
    = interp::n_states + interp::n_coeffs + interp::n_coeffs_int;
  //----------------------------------------------------------------------------
  void reset (crange<value_type> mem, uint n_channels)
  {
    _interp     = mem.cut_head (n_interp_states * n_channels);
    _delay_spls = mem.cut_head (1 * n_channels);
    crange_memset (_delay_spls, 0);
    base::reset (mem, n_channels);
  }
  //----------------------------------------------------------------------------
  static constexpr uint interp_overhead_elems (uint n_channels)
  {
    // + 1 for the delays
    return n_channels * (n_interp_states + 1);
  }
  //----------------------------------------------------------------------------
  // Stateful interpolators need this function to be called once setup and:
  // -for the "get" function to always have the same "delay_spls" value.
  // -for the "get" function to be called once for each "push" call.
  void reset_interpolator (float delay_spls, uint channel)
  {
    reset_interpolator (
      delay_spls, delay_spls - (float) ((uint) delay_spls), channel, true);
  }
  //----------------------------------------------------------------------------
  value_type get (uint channel)
  {
    assert (channel < n_channels());

    using V = make_vector_t<value_type>;

    return base::template get<interp::n_points, interp::x_offset> (
      get_delay_spls (channel), channel, [channel, this] (auto y, auto x) {
        return interp::tick (
          this->get_interp_coeffs (channel).template cast<const V>(),
          this->get_interp_states (channel).template cast<V>(),
          y,
          x);
      });
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  void reset_interpolator (
    float delay_spls,
    float delay_spls_frac,
    uint  channel,
    bool  clear_state)
  {
    assert (channel < n_channels());
    assert (delay_spls_frac == delay_spls - (float) ((uint) delay_spls));

    auto co = get_interp_coeffs (channel);
    auto st = get_interp_states (channel);

    using V = make_vector_t<value_type>;

    interp::reset_coeffs (co.cast (V {}), vec_set<V> (delay_spls_frac));
    if (clear_state) {
      interp::reset_states (st.cast (V {}));
    }
    if constexpr (is_vec_v<value_type>) {
      _delay_spls[channel] = make_vec<V> (delay_spls);
    }
    else {
      _delay_spls[channel] = delay_spls;
    }
  }
  //----------------------------------------------------------------------------
  crange<value_type> get_interp_coeffs (uint channel)
  {
    assert (channel < n_channels());
    auto chptr = &_interp[channel * (interp::n_states + interp::n_coeffs)];
    return {&chptr[0], interp::n_coeffs};
  }
  //----------------------------------------------------------------------------
  crange<value_type> get_interp_states (uint channel)
  {
    assert (channel < n_channels());
    auto chptr = &_interp[channel * (interp::n_states + interp::n_coeffs)];
    return {&chptr[interp::n_coeffs], interp::n_states};
  }
  //----------------------------------------------------------------------------
  float get_delay_spls (uint channel)
  {
    assert (channel < n_channels());
    if constexpr (is_vec_v<value_type>) {
      return (float) _delay_spls[channel][0];
    }
    else {
      return (float) _delay_spls[channel];
    }
  }
  //----------------------------------------------------------------------------
  // this one ignores "_delay_spls"
  value_type get_interpolated_unchecked (
    uint  spls_uint,
    float spls_frac,
    uint  channel)
  {
    assert (channel < n_channels());

    using V = make_vector_t<value_type>;

    return base::template get<interp::n_points, interp::x_offset> (
      spls_uint, spls_frac, channel, [channel, this] (auto y, auto x) {
        return interp::tick (
          this->get_interp_coeffs (channel).template cast<const V>(),
          this->get_interp_states (channel).template cast<V>(),
          y,
          x);
      });
  }
  //----------------------------------------------------------------------------
private:
  crange<value_type> _delay_spls;
  crange<value_type> _interp;
};
//------------------------------------------------------------------------------
} // namespace detail

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - random access interpolator.
// - optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
//------------------------------------------------------------------------------
template <class T, bool Interleaved = false, bool Use_pow2_sizes = true>
using interpolated_delay_line = detail::interpolated_delay_line<
  static_delay_line<T, Interleaved, Use_pow2_sizes>>;

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - stateful interpolator.
// - optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
//------------------------------------------------------------------------------
template <
  class T,
  class Interp        = thiran_interp<1>,
  bool Interleaved    = false,
  bool Use_pow2_sizes = true>
using statefully_interpolated_delay_line
  = detail::statefully_interpolated_delay_line<
    static_delay_line<T, Interleaved, Use_pow2_sizes>,
    Interp>;

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel
// - random access interpolator.
// - optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
// - convenience functions for taking modulation
//------------------------------------------------------------------------------
template <class T, bool Interleaved = false, bool Use_pow2_sizes = true>
class modulable_delay_line
  : private interpolated_delay_line<T, Interleaved, Use_pow2_sizes> {
private:
  using base = interpolated_delay_line<T, Interleaved, Use_pow2_sizes>;

public:
  //----------------------------------------------------------------------------
  using base::get;
  using base::get_raw;
  using base::n_channels;
  using base::push;
  using base::reset;
  using base::size;
  using value_type = typename base::value_type;
  //----------------------------------------------------------------------------
  // interpolated overload with modulation, just for convenience
  // Stateless_interp = a class on
  // "artv-common/dsp/own/parts/interpolation/stateless.hpp"
  template <class Stateless_interp>
  value_type get (float delay_spls, float mod_amt, float max_mod, uint channel)
  {
    // in case we are using a non interpolated delay line, no-op otherwise
    auto z         = (float) delay_spls;
    auto max_delay = size() - 1 - Stateless_interp::n_points - 1;
    z              = z + mod_amt * max_mod;
    z              = std::clamp<float> (z, 0.f, max_delay);
    return get<Stateless_interp> ((float) z, channel);
  }
  //----------------------------------------------------------------------------
};

namespace detail {
// Allpass interpolator taking the fractional part as a coefficient, as seen on
// Freeverb3
struct raw_allpass_interp {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points = 2;
  static constexpr uint x_offset = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V fractional)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {
    st[y1] = vec_set<V> ((vec_value_type_t<V>) 0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>         co,
    crange<V>               st,
    std::array<V, n_points> y_points,
    V                       x)
  {
    auto& z = y_points; // to avoid nomenclature clash
    V     y = z[0] + x * (z[1] - st[y1]);
    st[y1]  = y;
    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// 1-st order thiran, but with no previous inputs, as the delay line will pass
// them.
struct thiran_interp_1 {
  //----------------------------------------------------------------------------
  enum coeffs { a, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points = 2;
  static constexpr uint x_offset = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V fractional)
  {
    using T = vec_value_type_t<V>;
    // D parameter between 0.418 and 1.418
    V d = fractional + (T) 0.418;
#if 0
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    // TODO: Is this really faster than 1 divs?
    // TODO: It seems to cause noise/be broken.

    V v = ((T) 1. - d) * (T) (1. / 2.);

    V v_p2 = v * v;
    V v_p4 = v_p2 * v_p2;
    V v_p8 = v_p4 * v_p4;

    co[a] = v;
    co[a] *= ((T) 1 + v);
    co[a] *= ((T) 1 + v_p2);
    co[a] *= ((T) 1 + v_p4);
    co[a] *= ((T) 1 + v_p8);

#else
    co[a]  = (1 - d) / (1 + d);
#endif
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>         co,
    crange<V>               st,
    std::array<V, n_points> y_points,
    V                       x [[maybe_unused]])
  {
    auto& z = y_points;
    V     y = z[0] * co[a] + z[1] - co[a] * st[y1];
    st[y1]  = y;
    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// 2-st order thiran, but with no previous inputs, as the delay line will pass
// them -> Direct form 1 because of this.
struct thiran_interp_2_df1 {
  //----------------------------------------------------------------------------
  enum coeffs { a1, a2, n_coeffs };
  enum coeffs_int { n_coeffs_int };
  enum state { y1, y2, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points = 3;
  static constexpr uint x_offset = 0;
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V fractional)
  {
    using T = vec_value_type_t<V>;
    // D parameter between 1.5 and 2.5
    V d = fractional + (T) 1.5;
#if 0
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    // TODO: Is this really faster than 2 divs?
    // TODO: It seems to cause noise/be broken.

    V v1 = ((T) 3. - d) * (T) (1. / 4.);
    V v2 = ((T) 2. - d) * (T) (1. / 4.);

    V v1_p2 = v1 * v1;
    V v1_p4 = v1_p2 * v1_p2;
    V v1_p8 = v1_p4 * v1_p4;

    co[a1] = (T) -2 * ((T) 0.25 - v1);
    co[a1] *= ((T) 1 + v1);
    co[a1] *= ((T) 1 + v1_p2);
    co[a1] *= ((T) 1 + v1_p4);
    co[a1] *= ((T) 1 + v1_p8);

    V v2_p2 = v2 * v2;
    V v2_p4 = v2_p2 * v2_p2;
    V v2_p8 = v2_p4 * v2_p4;

    co[a2] = co[a1] * (T) -0.5 * ((T) 0.25 - v2);
    co[a2] *= ((T) 1 + v2);
    co[a2] *= ((T) 1 + v2_p2);
    co[a2] *= ((T) 1 + v2_p4);
    co[a2] *= ((T) 1 + v2_p8);

#else
    co[a1] = -(d - 2) / (d + 1);
    co[a2] = ((d - 1) * (d - 2)) / ((d + 1) * (d + 2));
#endif
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V>         co,
    crange<V>               st,
    std::array<V, n_points> y_points,
    V                       x [[maybe_unused]])
  {
    auto& z = y_points;

    V y = co[a2] * z[0] + co[a1] * z[1] + z[2];
    y -= co[a1] * st[y1] + co[a2] * st[y2];
    st[y2] = st[y1];
    st[y1] = y;
    return y;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// Slowly modulable allpass. The main difference with
// "statefully_interpolated_delay_line" is that it is setting the interpolator
// states to linear interpolation and optionally running an extra  warmup. This
// is specialized for the allpass classes above, hence it is on the ::detail
// namespace.

template <
  class T,
  uint Allpass_warmup_spls,
  bool Interleaved,
  bool Use_pow2_sizes,
  class Interp>
class modulable_allpass_delay_line
  : private artv::statefully_interpolated_delay_line<
      T,
      Interp,
      Interleaved,
      Use_pow2_sizes> {
  using base = artv::
    statefully_interpolated_delay_line<T, Interp, Interleaved, Use_pow2_sizes>;

public:
  //----------------------------------------------------------------------------
  using base::get_raw;
  using base::interp_overhead_elems;
  using base::n_channels;
  using base::push;
  using base::reset;
  using base::size;
  using raw_interp = typename base::interp;
  struct interp {
    static constexpr uint n_points = raw_interp::n_points + Allpass_warmup_spls;
    static constexpr uint x_offset = raw_interp::x_offset;
  };
  using value_type                      = typename base::value_type;
  static constexpr uint n_interp_states = base::n_interp_states;
  static constexpr uint n_warmup_spls   = Allpass_warmup_spls;
  //----------------------------------------------------------------------------
  static constexpr uint interp_overhead_elems (uint n_channels)
  {
    return base::interp_overhead_elems (n_channels);
  }
  //----------------------------------------------------------------------------
  void set_resync_delta (float spls) { _resync_delta = spls; }
  //----------------------------------------------------------------------------
  value_type get (float delay_spls, uint channel)
  {
    assert (channel < n_channels());

    auto diff     = delay_spls - this->get_delay_spls (channel);
    auto diff_abs = abs (diff);
    if (unlikely (diff_abs == 0.f)) {
      // no modulation, not the main use case for this...
      return base::get (channel);
    }

    uint  spls = (uint) delay_spls;
    float frac = delay_spls - (float) spls;
    this->reset_interpolator (delay_spls, frac, channel, false);

    if (diff_abs >= _resync_delta) {
      // hack the previous states to linear interpolation as a starting point
      std::array<value_type, raw_interp::n_points> z;

      bool delay_is_increasing = diff > 0.f;
      if (delay_is_increasing) {
        // move in increasing direction: away from the write pointer (idx 0)
        for (uint i = z.size() - 1; i < z.size(); --i) {
          z[i] = get_raw (spls - n_warmup_spls - i, channel);
        }
      }
      else {
        // move in decreasing direction: towards the write pointer (idx 0)
        for (uint i = z.size() - 1; i < z.size(); --i) {
          z[i] = get_raw (spls + n_warmup_spls + i, channel);
        }
      }
      auto st = this->get_interp_states (channel);
      for (uint i = 0; i < z.size() - 1; ++i) {
        if constexpr (is_vec_v<value_type>) {
          using builtin = vec_value_type_t<value_type>;
          st[i]         = linear_interp::tick (
            make_array (z[i], z[i + 1]), vec_set<value_type> ((builtin) frac));
        }
        else {
          st[i] = linear_interp::tick (make_array (z[i], z[i + 1]), frac);
        }
      }
      // run some warmup samples on the allpass (if any)
      if (delay_is_increasing) {
        // move in increasing direction: away from the write pointer (idx 0)
        for (uint i = 0; i < n_warmup_spls; ++i) {
          base::get_interpolated_unchecked (
            spls - n_warmup_spls + i, frac, channel);
        }
      }
      else {
        // move in decreasing direction: towards the write pointer (idx 0)
        for (uint i = 0; i < n_warmup_spls; ++i) {
          base::get_interpolated_unchecked (
            spls + n_warmup_spls - i, frac, channel);
        }
      }
    }
    return base::get_interpolated_unchecked (spls, frac, channel);
  }
  //----------------------------------------------------------------------------
private:
  float _resync_delta {};
};
} // namespace detail

// Adding (slow) modulation to allpass interpolated delay lines
template <
  class T,
  uint Allpass_warmup_n_spls = 0,
  bool Interleaved           = false,
  bool Use_pow2_sizes        = true>
using modulable_raw_allpass_delay_line = detail::modulable_allpass_delay_line<
  T,
  Allpass_warmup_n_spls,
  Interleaved,
  Use_pow2_sizes,
  detail::raw_allpass_interp>;

template <
  class T,
  uint Allpass_warmup_n_spls = 0,
  bool Interleaved           = false,
  bool Use_pow2_sizes        = true>
using modulable_thiran1_delay_line = detail::modulable_allpass_delay_line<
  T,
  Allpass_warmup_n_spls,
  Interleaved,
  Use_pow2_sizes,
  detail::thiran_interp_1>;

template <
  class T,
  uint Allpass_warmup_n_spls = 0,
  bool Interleaved           = false,
  bool Use_pow2_sizes        = true>
using modulable_thiran2_delay_line = detail::modulable_allpass_delay_line<
  T,
  Allpass_warmup_n_spls,
  Interleaved,
  Use_pow2_sizes,
  detail::thiran_interp_2_df1>;
//------------------------------------------------------------------------------
// - single dynamic allocation
// - multichannel
// - random access.
// - power of 2 capacity.
// - not interleaved
// - same sized channels
// - highest indexes = older samples
template <class T, bool Interleaved, bool Use_pow2_sizes = true>
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
    if constexpr (Use_pow2_sizes) {
      max_delay = pow2_round_ceil (max_delay);
    }
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
  static_delay_line<T, Interleaved, Use_pow2_sizes> _z;
  T*                                                _mem = nullptr;
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
    V y  = in - y1 * _gain;
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
