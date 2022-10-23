#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

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
  void reset (xspan<V> mem)
  {
    assert (is_pow2 (mem.size()));
    _mem  = mem.data();
    _mask = mem.size() - 1;
    clear();
  };
  //----------------------------------------------------------------------------
  void clear()
  {
    assert (_mem);
    _pos = 0;
    memset (_mem, 0, sizeof *_mem * size());
  }
  //----------------------------------------------------------------------------
  // idx from newest to oldest, so idx=0 -> z, idx=1 -> z-1, etc.
  V get (uint idx) const { return get_abs (_pos - idx); }
  V operator[] (uint idx) const { return get (idx); }
  //----------------------------------------------------------------------------
  void push (V const& value)
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
  void reset (xspan<V> mem)
  {
    assert (mem.size());
    _mem  = mem.data();
    _size = mem.size();
    clear();
  };
  //----------------------------------------------------------------------------
  void clear()
  {
    assert (_mem);
    _pos = 0;
    memset (_mem, 0, sizeof *_mem * size());
  }
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
  void push (V const& value)
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
  constexpr void reset (xspan<T> mem, uint n_channels)
  {
    uint size = mem.size() / n_channels;
    assert ((mem.size() % n_channels) == 0);
    assert (size != 0 && is_pow2 (size));

    _z          = mem.data();
    _mask       = size - 1;
    _n_channels = n_channels;
    clear();
  }
  //----------------------------------------------------------------------------
  constexpr void clear()
  {
    assert (_z);
    _pos = 0;
    memset (_z, 0, size() * _n_channels * sizeof *_z);
  }
  //----------------------------------------------------------------------------
  constexpr uint size() const { return _mask + 1; }
  //----------------------------------------------------------------------------
  constexpr uint n_channels() const { return _n_channels; }
  //----------------------------------------------------------------------------
  // allows querying how many elements of the "value_type" type have to be
  // allocated and passed to "reset" for a given delay line size and number of
  // channels.
  //----------------------------------------------------------------------------
  static constexpr uint n_required_elems (uint n_spls, uint n_channels)
  {
    return (n_spls * n_channels);
  }
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
  constexpr void reset (xspan<T> mem, uint n_channels)
  {
    uint size = mem.size() / n_channels;
    assert ((mem.size() % n_channels) == 0);
    assert (size != 0);

    _z          = mem.data();
    _n_channels = n_channels;
    _size       = size;
    clear();
  }
  //----------------------------------------------------------------------------
  constexpr void clear()
  {
    assert (_z);
    _pos = 0;
    memset (_z, 0, size() * _n_channels * sizeof *_z);
  }
  //----------------------------------------------------------------------------
  constexpr uint size() const { return _size; }
  //----------------------------------------------------------------------------
  constexpr uint n_channels() const { return _n_channels; }
  //----------------------------------------------------------------------------
  // allows querying how many elements of the "value_type" type have to be
  // allocated and passed to "reset" for a given delay line size and number of
  // channels.
  //----------------------------------------------------------------------------
  static constexpr uint n_required_elems (uint n_spls, uint n_channels)
  {
    return (n_spls * n_channels);
  }
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
// - multichannel, but all of them of the same maximum size.
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
  using base::clear;
  using base::n_channels;
  using base::n_required_elems;
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
  constexpr void push (xspan<value_type> const row)
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
  using base::clear;
  using base::n_channels;
  using base::n_required_elems;
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
  constexpr void push (xspan<value_type> const row)
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
  using base::clear;
  using base::n_channels;
  using base::n_required_elems;
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
  constexpr void push (xspan<value_type> const row)
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
  using base::clear;
  using base::n_channels;
  using base::n_required_elems;
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
  constexpr void push (xspan<value_type> const row)
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
// "Interp "a class on
// "artv-common/dsp/own/parts/interpolation/stateless.hpp"
// "Delay_line_base" = the underlying buffer type, e.g. "static_delay_line".
// "Interp" = interpolator type
// "External_global_interp_coeffs"= used for reusing interpolator coefficients,
// e.g. for the case of sharing sinc interpolator tables.
template <
  class Delay_line_base,
  class Interp,
  bool External_global_interp_coeffs = false>
class interpolated_delay_line : private Delay_line_base {
private:
  using base = Delay_line_base;

public:
  //----------------------------------------------------------------------------
  using base::clear;
  using base::n_channels;
  using base::push;
  using base::size;

  using value_type                           = typename base::value_type;
  using interp                               = Interp;
  static constexpr bool interp_global_co_ext = External_global_interp_coeffs;
  static_assert (!interp_global_co_ext || interp::coeffs_are_global);

  // As of now the additional interpolation boundaries can't be abstacted
  // away because it conflicts the requirement of the memory been allocated
  // externally, especially when it has to be a power of 2. It would be nice
  // to find a solution and being able to always retrieve e.g. position 0.
  //
  // In practice on reverbs and delays it is not a big deal but it might be
  // for other applications.

  //----------------------------------------------------------------------------
  static constexpr uint min_delay_spls() { return interp::x_offset; };
  //----------------------------------------------------------------------------
  uint max_delay_spls()
  {
    return size() - (interp::n_points + interp::x_offset);
  };
  //----------------------------------------------------------------------------
  static constexpr uint min_size_spls() { return interp::n_points; };
  //----------------------------------------------------------------------------
  // allows querying how many elements of the "value_type" type have to be
  // allocated and passed to "reset" for a given delay line size and number of
  // channels.
  static constexpr uint n_required_elems (uint size, uint n_channels)
  {
    uint n            = 0;
    uint n_coeff_sets = interp::coeffs_are_global ? 1 : n_channels;

    if constexpr (storage::size_co_elem == 1 && storage::size_st_elem == 1) {
      // no padding coeff/state-wise between channels but one padding at the
      // end might apply.
      n += div_ceil (
        interp::n_coeffs * n_coeff_sets + interp::n_states * n_channels,
        storage::vec_size);
    }
    else {
      n += (storage::size_co_padded / storage::vec_size) * n_coeff_sets;
      n += (storage::size_st_padded / storage::vec_size) * n_channels;
    }
    n += base::n_required_elems (size, n_channels);
    return n;
  }
  //----------------------------------------------------------------------------
  // distributes the memory, but if the interpolators require initialization
  // they will have to be reset separately.
  void reset (xspan<value_type> mem, uint n_channels)
  {
    if constexpr (interp::n_coeffs > 0 || interp::n_states > 0) {
      uint n_coeff_sets = interp::coeffs_are_global ? 1 : n_channels;
      uint n_cut        = 0;
      if constexpr (storage::size_co_elem == 1 && storage::size_st_elem == 1) {
        // no padding coeff/state-wise between channels but one padding at the
        // end might apply.
        n_cut += round_ceil (
          interp::n_coeffs * n_coeff_sets + interp::n_states * n_channels,
          storage::vec_size);
      }
      else {
        n_cut += storage::size_co_padded * n_coeff_sets;
        n_cut += storage::size_st_padded * n_channels;
      }
      _mem = mem.template cast<builtin_type>().cut_head (n_cut);
      xspan_memset (_mem, 0);
      mem = mem.template cast<builtin_type>()
              .advanced (n_cut)
              .template cast<value_type>();
    }
    assert ((mem.size() / n_channels) >= min_size_spls());
    base::reset (mem, n_channels);
  }
  //----------------------------------------------------------------------------
  template <
    class... Args,
    std::enable_if_t<
      (sizeof...(Args) == sizeof...(Args)) && !interp_global_co_ext>* = nullptr>
  void reset_interpolator (uint channel, bool reset_state, Args&&... interp_arg)
  {
    assert (channel < n_channels());
    if (channel == 0 || !interp::coeffs_are_global) {
      // only reset the interpolator by passing channel 0 when the coefficients
      // are global. This is to avoid unnecessary coefficient recomputations
      // while still allow resetting the states for each individual channel.
      interp::reset_coeffs (
        get_interp_coeffs (channel), std::forward<Args> (interp_arg)...);
    }
    if (reset_state) {
      interp::reset_states (get_interp_states (channel));
    }
  }
  //----------------------------------------------------------------------------
  template <
    class T,
    std::enable_if_t<std::is_same_v<T, T> && interp_global_co_ext>* = nullptr>
  void reset_interpolator (
    uint           channel,
    bool           reset_state,
    xspan<T const> external_coeffs)
  {
    static_assert (std::is_same_v<T, coeffs_type>);
    if (channel == 0) {
      // only reset the interpolator by passing channel 0 when the coefficients
      // are global. This is to keep simmetry with "reset_interpolator" when
      // "interp_global_co_ext" is false.
      static_assert (interp::coeffs_are_global);
      _ext_coeffs = external_coeffs;
    }
    if (reset_state) {
      interp::reset_states (get_interp_states (channel));
    }
  }
  //----------------------------------------------------------------------------
  // integer index access
  value_type get_raw (uint sample, uint channel)
  {
    return base::get (sample, channel);
  }
  //----------------------------------------------------------------------------
  value_type get (float delay_spls, uint channel)
  {
    auto spls_uint = (uint) delay_spls;
    return get (spls_uint, delay_spls - spls_uint, channel);
  }
  //----------------------------------------------------------------------------
  value_type get (uint delay_spls_int, float delay_spls_frac, uint channel)
  {
    assert (channel < n_channels());
    auto coeffs = get_interp_coeffs (channel);
    auto states = get_interp_states (channel);

    std::array<value_type, interp::n_points> y;
    // would clamping instead of an assert be better?
    assert (
      std::clamp<float> (
        delay_spls_frac + (float) delay_spls_int,
        min_delay_spls(),
        max_delay_spls())
      == (delay_spls_frac + (float) delay_spls_int));

    // e.g. Catmull-Rom interpolates between the 2 central points of 4. Sinc
    // interpolators also have offset
    delay_spls_int -= interp::x_offset;

    for (uint i = 0; i < y.size(); ++i) {
      y[i] = get_raw (delay_spls_int + i, channel);
    }
    auto ret = interp::template tick<value_type> (
      coeffs, states, y, make_vec<value_type> (delay_spls_frac));
    return ret;
  }
  //----------------------------------------------------------------------------
  void get (xspan<value_type> out, xspan<float> n_spls)
  {
    assert (out.size() >= n_channels());
    assert (n_spls.size() >= n_channels());
    for (uint i = 0; i < n_channels(); ++i) {
      out[i] = get (n_spls[i], i);
    }
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  using coeffs_type = std::conditional_t<
    (is_vec_v<value_type> && (interp::n_coeffs > 0) && !interp::coeffs_are_vec),
    vec_value_type_t<value_type>,
    value_type>;

  using states_type = std::conditional_t<
    (is_vec_v<value_type> && (interp::n_states > 0) && !interp::states_are_vec),
    vec_value_type_t<value_type>,
    value_type>;
  //----------------------------------------------------------------------------
  auto get_interp_coeffs (uint channel)
  {
    if constexpr (interp::n_coeffs == 0) {
      return xspan<coeffs_type> {};
    }
    else if constexpr (interp_global_co_ext) {
      return _ext_coeffs;
    }
    else {
      auto ret = _mem;
      if (!interp::coeffs_are_global) {
        ret.cut_head (channel * storage::channel_size);
      }
      ret = ret.get_head (storage::size_co_padded);
      return ret.template cast<coeffs_type>();
    }
  }
  //----------------------------------------------------------------------------
  xspan<states_type> get_interp_states (uint channel)
  {
    if constexpr (interp::n_states == 0) {
      return {};
    }
    else {
      auto ret = _mem;
      if (interp::coeffs_are_global) {
        ret.cut_head (storage::size_co_padded);
        ret.cut_head (channel * storage::channel_size);
        return ret.template cast<states_type>();
      }
      else {
        ret.cut_head (channel * storage::channel_size);
        ret.cut_head (storage::size_co_padded);
        return ret.template cast<states_type>();
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  struct storage {
    // size of the vector type (if any)
    static constexpr uint vec_size
      = is_vec_v<value_type> ? vec_traits_t<value_type>::size : 1;

    // size in number of floats of one coefficient/state element
    static constexpr uint size_co_elem = (interp::coeffs_are_vec ? vec_size : 1)
      * (interp_global_co_ext ? 0 : 1);
    static constexpr uint size_st_elem = interp::states_are_vec ? vec_size : 1;

    // unpadded size in number of floats of all coefficient/state elements
    static constexpr uint size_co = interp::n_coeffs * size_co_elem;
    static constexpr uint size_st = interp::n_states * size_st_elem;

    // padded size in number of floats of all coefficient/state elements
    static constexpr uint size_co_padded
      = (size_co_elem == size_st_elem || size_co_elem != 1)
      ? size_co
      : round_ceil (size_co, size_st_elem);

    static constexpr uint size_st_padded
      = (size_co_elem == size_st_elem || size_st_elem != 1 || size_co_elem == 0)
      ? size_st
      : round_ceil (size_st, size_co_elem);

    // the size of each channel
    static constexpr uint channel_size = interp::coeffs_are_global
      ? size_st_padded
      : (size_co_padded + size_st_padded);
  };

  using builtin_type = std::conditional_t<
    is_vec_v<value_type>,
    vec_value_type_t<value_type>,
    value_type>;

  using mem_storage_type = std::conditional_t<
    interp::n_coeffs == 0 && interp::n_states == 0,
    null_type,
    xspan<builtin_type>>;

  using ext_interp_coeffs_type = std::
    conditional_t<interp_global_co_ext, xspan<coeffs_type const>, null_type>;

  mem_storage_type       _mem {};
  ext_interp_coeffs_type _ext_coeffs {};
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------

// A class mainly done to run thiran/allpass interpolators unmodulated. It
// stores the fractional delay and makes obvious that the delay line is not
// random access, as the get method doesn't have a sample count parameter.
//
// "Delay_line_base" = the underlying buffer type, e.g. "static_delay_line".
// "Interpolation" one of the classes on
// "artv-common/dsp/own/parts/interpolation/stateful.hpp"

template <class Delay_line_base, class Interpolation>
class statefully_interpolated_delay_line
  : private interpolated_delay_line<Delay_line_base, Interpolation> {
private:
  using base = interpolated_delay_line<Delay_line_base, Interpolation>;

public:
  // TODO Clear
  //----------------------------------------------------------------------------
  using base::get_raw;
  using base::n_channels;
  using base::push;
  using base::size;

  using value_type = typename base::value_type;
  using interp     = Interpolation;
  //----------------------------------------------------------------------------
  static constexpr uint min_delay_spls() { return base::min_delay_spls(); };
  constexpr uint        max_delay_spls() { return base::max_delay_spls(); };
  //----------------------------------------------------------------------------
  // this is for thiran interpolators, which take the previous inputs from the
  // delay line itself, so the minimum size has to account for the states
  static constexpr uint min_size_spls()
  {
    return base::min_size_spls() + interp::n_states;
  };
  //----------------------------------------------------------------------------
  static constexpr uint n_required_elems (uint n_spls, uint n_channels)
  {
    uint n = 0;
    if constexpr (is_vec_v<value_type>) {
      n += div_ceil (n_channels, vec_traits_t<value_type>::size);
    }
    else {
      n += interp::n_channels;
    }
    n += base::n_required_elems (n_spls, n_channels);
    return n;
  }
  //----------------------------------------------------------------------------
  void clear()
  {
    for (uint i = 0; i < n_channels(); ++i) {
      reset_interpolator (0.f, 0.f, i, true);
    }
    base::clear();
  }
  //----------------------------------------------------------------------------
  void reset (xspan<value_type> mem, uint n_channels)
  {
    if constexpr (is_vec_v<value_type>) {
      _delay_spls = mem.template cast<delay_type>().cut_head (n_channels);
      mem.cut_head (div_ceil (n_channels, vec_traits_t<value_type>::size));
    }
    else {
      _delay_spls = mem.cut_head (n_channels);
    }
    xspan_memset (_delay_spls, 0);
    base::reset (mem, n_channels);
  }
  //----------------------------------------------------------------------------
  // Stateful interpolators need this function to be called once setup and:
  // -for the "get" function to always have the same "delay_spls" value.
  // -for the "get" function to be called once for each "push" call. A sample
  // has to be taken out for each read (case example Thiran interp).
  void reset_interpolator (float delay_spls, uint channel)
  {
    reset_interpolator (
      delay_spls, delay_spls - (float) ((uint) delay_spls), channel, true);
  }
  //----------------------------------------------------------------------------
  value_type get (uint channel)
  {
    assert (channel < n_channels());
    return base::get (get_delay_spls (channel), channel);
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

    using V = make_vector_t<value_type>;
    base::reset_interpolator (
      channel, clear_state, vec_set<V> (delay_spls_frac));
    _delay_spls[channel] = delay_spls;
  }
  //----------------------------------------------------------------------------
  float get_delay_spls (uint channel)
  {
    assert (channel < n_channels());
    return (float) _delay_spls[channel];
  }
  //----------------------------------------------------------------------------
  // this one ignores "_delay_spls"
  value_type get_interpolated_unchecked (
    uint  spls_uint,
    float spls_frac,
    uint  channel)
  {
    assert (channel < n_channels());
    return base::get (spls_uint, spls_frac, channel);
  }
  //----------------------------------------------------------------------------
  using base::get_interp_coeffs;
  using base::get_interp_states;
  //----------------------------------------------------------------------------
private:
  using delay_type = std::conditional_t<
    is_vec_v<value_type>,
    vec_value_type_t<value_type>,
    value_type>;

  xspan<delay_type> _delay_spls;
};
//------------------------------------------------------------------------------
} // namespace detail

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel, but all of them of the same maximum size.
// - random access interpolator.
// - optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
//------------------------------------------------------------------------------
template <
  class T,
  class Interp                       = linear_interp,
  bool Interleaved                   = false,
  bool Use_pow2_sizes                = true,
  bool External_global_interp_coeffs = false>
using interpolated_delay_line = detail::interpolated_delay_line<
  static_delay_line<T, Interleaved, Use_pow2_sizes>,
  Interp,
  External_global_interp_coeffs>;

//------------------------------------------------------------------------------
// - not managed memory
// - multichannel, but all of them of the same maximum size.
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
// - multichannel, but all of them of the same maximum size.
// - random access interpolator.
// - optional power of 2 capacity optimization.
// - channels might or not be interleaved
// - same sized channels
// - highest indexes = older samples
// - fractional delay line support
// - convenience functions for taking modulation
//------------------------------------------------------------------------------
template <
  class T,
  class Interp        = linear_interp,
  bool Interleaved    = false,
  bool Use_pow2_sizes = true>
class modulable_delay_line
  : private interpolated_delay_line<T, Interp, Interleaved, Use_pow2_sizes> {
private:
  using base = interpolated_delay_line<T, Interp, Interleaved, Use_pow2_sizes>;

public:
  //----------------------------------------------------------------------------
  using base::clear;
  using base::get;
  using base::get_raw;
  using base::n_channels;
  using base::n_required_elems;
  using base::push;
  using base::reset;
  using base::size;

  using interp     = Interp;
  using value_type = typename base::value_type;
  //----------------------------------------------------------------------------
  // interpolated overload with modulation, just for convenience
  // Stateless_interp = a class on
  // "artv-common/dsp/own/parts/interpolation/stateless.hpp"
  value_type get (float delay_spls, float mod_amt, float max_mod, uint channel)
  {
    // in case we are using a non interpolated delay line, no-op otherwise
    auto z         = (float) delay_spls;
    auto max_delay = size() - 1 - interp::n_points - 1;
    z              = z + mod_amt * max_mod;
    z              = std::clamp<float> (z, 0.f, max_delay);
    return base::get ((float) z, channel);
  }
  //----------------------------------------------------------------------------
};

namespace detail {
// Allpass interpolator taking the fractional part as a coefficient, as seen on
// Freeverb3
struct raw_allpass_interp {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs };
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points          = 2;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = false; // 1 set per channel
  static constexpr bool states_are_vec    = true;
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V fractional)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {
    st[y1] = vec_set<V> ((vec_value_type_t<V>) 0);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (
    xspan<V const>          co,
    xspan<V>                st,
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
  enum state { y1, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points          = 2;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = false; // 1 set per channel
  static constexpr bool states_are_vec    = true;
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V freq, vec_value_type_t<V> srate)
  {
    using T = vec_value_type_t<V>;
    auto d  = vec_tan (M_PI * freq / srate);
    co[a]   = ((T) 1 - d) / ((T) 1 + d);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V fractional)
  {
    using T = vec_value_type_t<V>;
    // D parameter between 0.418 and 1.418
    V d = fractional + (T) 0.418;
#if 1
    co[a] = ((T) 1 - d) / ((T) 1 + d);
#else
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    V v = ((T) 1. - d) * (T) (1. / 2.);

    V v_p2 = v * v;
    V v_p4 = v_p2 * v_p2;
    V v_p8 = v_p4 * v_p4;

    co[a] = v;
    co[a] *= ((T) 1 + v);
    co[a] *= ((T) 1 + v_p2);
    co[a] *= ((T) 1 + v_p4);
    co[a] *= ((T) 1 + v_p8);
#endif
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (
    xspan<V const>          co,
    xspan<V>                st,
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
  enum state { y1, y2, n_states };
  //----------------------------------------------------------------------------
  static constexpr uint n_points          = 3;
  static constexpr uint x_offset          = 0;
  static constexpr bool coeffs_are_vec    = true;
  static constexpr bool coeffs_are_global = false; // 1 set per channel
  static constexpr bool states_are_vec    = true;
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V freq, vec_value_type_t<V> srate)
  {
    using T = vec_value_type_t<V>;
    auto d  = vec_tan (M_PI * freq / srate);
    co[a1]  = (1 - d) / (1 + d);
    co[a2]  = ((d - 1) * (d - 2)) / ((d + 1) * (d + 2));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_coeffs (xspan<V> co, V fractional)
  {
    using T = vec_value_type_t<V>;
    V d     = fractional + (T) 1.403;
#if 1
    co[a1] = -(d - 2) / (d + 1);
    co[a2] = ((d - 1) * (d - 2)) / ((d + 1) * (d + 2));
#else
    // See https://dafx09.como.polimi.it/proceedings/papers/paper_72.pdf
    // chapter 6.
    // These measured worse than the divisions on a Ryzen7-5800x running 8x16
    // delay lines.
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
#endif
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static void reset_states (xspan<V> st)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_floatpt_vec_t<V>* = nullptr>
  static V tick (
    xspan<V const>          co,
    xspan<V>                st,
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
// "statefully_interpolated_delay_line" is that it allows to modulate the delay
// time.

template <
  class T,
  uint Allpass_warmup_spls, // warmup n of samples when the delay changes
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
  using base::max_delay_spls;
  using base::min_delay_spls;
  using base::n_channels;
  using base::n_required_elems;
  using base::push;
  using base::reset;
  using base::size;
  using raw_interp = typename base::interp;
  struct interp {
    static constexpr uint n_points = raw_interp::n_points + Allpass_warmup_spls;
    static constexpr uint x_offset = raw_interp::x_offset;
  };
  using value_type                    = typename base::value_type;
  static constexpr uint n_warmup_spls = Allpass_warmup_spls;
  //----------------------------------------------------------------------------
  // this is for thiran interpolators, which take the previous inputs from the
  // delay line itself, so the minimum size has to account for the states
  //----------------------------------------------------------------------------
  static constexpr uint min_size_spls()
  {
    return base::min_size_spls() + n_warmup_spls;
  };
  //----------------------------------------------------------------------------
  // Not so sure of the resync feature, on delays with low modulation it has
  // been good enough to set this value high and only run with recomputation of
  // coefficients.
  void set_resync_delta (float spls) { _resync_delta = abs (spls); }
  //----------------------------------------------------------------------------
  // The delay time sample difference against last call's delay time that
  // triggers a coefficient reset on the interpolator. Useful e.g. when the
  // sample delay is externally smoothed by a lowpass never hitting the
  // target value. Very small values are required to not be audible, at maximum
  // probably around "0.00001f".
  void set_interp_delta (float spls) { _epsilon = abs (spls); }
  //----------------------------------------------------------------------------
  // when calling blockwise a block of samples is fetched, processed, stored
  // and then a block of the same size is inserted.
  //
  // As the interpolators may reset each time "delay_spls" changes, the
  // "delay_spls_offset" allows distinguishing between modulation and block
  // access and allows the function to act accordingly.
  value_type get (float delay_spls, uint channel, uint delay_spls_offset = 0)
  {
    assert (channel < n_channels());

    auto del_prev    = this->get_delay_spls (channel);
    auto epsilon     = delay_spls - del_prev;
    auto epsilon_abs = abs (epsilon);

    if (unlikely (epsilon_abs <= _epsilon)) {
      // Not the main use case for this
      uint  spls = (uint) delay_spls;
      float frac = delay_spls - (float) spls;
      spls -= delay_spls_offset;
      return base::get_interpolated_unchecked (spls, frac, channel);
    }
    uint  spls = (uint) delay_spls;
    float frac = delay_spls - (float) spls;
    this->reset_interpolator (delay_spls, frac, channel, false);
    spls -= delay_spls_offset;
    assert (spls >= min_delay_spls());

    if (epsilon_abs >= _resync_delta) {
      // hack the previous states to linear interpolation as a starting point
      std::array<value_type, raw_interp::n_points> z;

      bool delay_is_increasing = epsilon > 0.f;
      // TODO: is the direction thing correct?
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
  constexpr void get (xspan<value_type> dst, xspan<float> del_spls)
  {
    assert (dst.size() >= n_channels());
    assert (del_spls.size() >= n_channels());
    for (uint i = 0; i < n_channels(); ++i) {
      dst[i] = get (del_spls[i], i);
    }
  }
  //----------------------------------------------------------------------------
private:
  float _resync_delta {50.f};
  float _epsilon {1e-30};
};
} // namespace detail

// Adding (slow) modulation to allpass interpolated delay lines.
// These Thirans need good DC blocking if used on feedback loops.
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
// - multichannel, but all of them of the same maximum size.
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
    _z.reset (xspan {_mem, elems}, n_channels);
  }
  //----------------------------------------------------------------------------
  T get (time_type sample, uint channel) { return _z.get (sample, channel); }
  //----------------------------------------------------------------------------
  void push (xspan<T> const x) { _z.push (x); }
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
template <class T, uint N, bool is_pow2>
class multisized_delay_line;

template <class T, uint N>
class multisized_delay_line<T, N, true> {
public:
  using value_type                          = T;
  static constexpr uint n_channels          = N;
  static constexpr bool requires_pow2_sizes = true;
  //----------------------------------------------------------------------------
  static_assert (std::is_floating_point_v<T>);
  //----------------------------------------------------------------------------
  template <class U>
  xspan<value_type> reset (xspan<value_type> mem, xspan<U const> sizes)
  {
    static_assert (std::is_unsigned_v<U>);

    assert (sizes.size() >= n_channels);
    for (uint i = 0; i < _delays.size(); ++i) {
      assert (is_pow2 (sizes[i]));
      _delays[i] = delay_data {mem.cut_head (sizes[i]).data(), sizes[i] - 1};
    }
    return mem; // return the remainder
  }
  //----------------------------------------------------------------------------
  constexpr void push (xspan<value_type> const row)
  {
    assert (row.size() >= n_channels);
    ++_pos;
    for (uint i = 0; i < n_channels; ++i) {
      auto& del                = _delays[i];
      del.ptr[del.mask & _pos] = row[i];
    }
  }
  //----------------------------------------------------------------------------
  constexpr value_type get (uint del_spls, uint chnl)
  {
    assert (chnl < n_channels);
    auto& del = _delays[chnl];
    assert (del_spls <= del.mask); // unintended wraparound?
    return del.ptr[del.mask & (_pos - del_spls)];
  }
  //----------------------------------------------------------------------------
  constexpr void get (xspan<value_type> dst, xspan<uint> del_spls)
  {
    assert (dst.size() >= n_channels);
    assert (del_spls.size() >= n_channels);
    for (uint i = 0; i < n_channels; ++i) {
      dst[i] = get (del_spls[i], i);
    }
  }
  //----------------------------------------------------------------------------
private:
  struct delay_data {
    value_type* ptr;
    uint        mask;
  };
  //----------------------------------------------------------------------------
  std::array<delay_data, N> _delays {};
  uint                      _pos {};
};
//------------------------------------------------------------------------------
template <class T, uint N>
class multisized_delay_line<T, N, false> {
public:
  using value_type                          = T;
  static constexpr uint n_channels          = N;
  static constexpr bool requires_pow2_sizes = false;
  //----------------------------------------------------------------------------
  static_assert (std::is_floating_point_v<T>);
  //----------------------------------------------------------------------------
  template <class U>
  xspan<value_type> reset (xspan<value_type> mem, xspan<U const> sizes)
  {
    static_assert (std::is_unsigned_v<U>);

    assert (sizes.size() >= n_channels);
    for (uint i = 0; i < _delays.size(); ++i) {
      _delays[i] = delay_data {mem.cut_head (sizes[i]).data(), sizes[i], 0};
    }
    return mem; // return the remainder
  }
  //----------------------------------------------------------------------------
  constexpr void push (xspan<value_type> const row)
  {
    assert (row.size() >= n_channels);
    for (uint i = 0; i < n_channels; ++i) {
      auto& del = _delays[i];
      ++del.pos;
      del.pos          = del.pos < del.size ? del.pos : 0;
      del.ptr[del.pos] = row[i];
    }
  }
  //----------------------------------------------------------------------------
  constexpr value_type get (uint del_spls, uint chnl)
  {
    assert (chnl < n_channels);
    auto& del = _delays[chnl];
    assert (del_spls < del.size); // unintended wraparound?
    del_spls = del.pos - del_spls;
    del_spls += del_spls > del.size ? del.size : 0;
    return del.ptr[del_spls];
  }
  //----------------------------------------------------------------------------
  constexpr void get (xspan<value_type> dst, xspan<uint> del_spls)
  {
    assert (dst.size() >= n_channels);
    assert (del_spls.size() >= n_channels);
    for (uint i = 0; i < n_channels; ++i) {
      dst[i] = get (del_spls[i], i);
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_required_elems (xspan<uint const> sizes)
  {
    return std::accumulate (sizes.begin(), sizes.end(), 0);
  }
  //----------------------------------------------------------------------------
private:
  struct delay_data {
    value_type* ptr;
    uint        size;
    uint        pos;
  };
  //----------------------------------------------------------------------------
  std::array<delay_data, N> _delays {};
};
//------------------------------------------------------------------------------
template <class T, uint N, class Stateless_interp, bool is_pow2>
class interpolated_multisized_delay_line
  : private multisized_delay_line<T, N, is_pow2> {
private:
  using base = multisized_delay_line<T, N, is_pow2>;
  //----------------------------------------------------------------------------
public:
  using value_type                 = typename base::value_type;
  static constexpr auto n_channels = base::n_channels;
  using interp                     = Stateless_interp;

  using base::push;
  using base::reset;
  //----------------------------------------------------------------------------
  // Stateless_interp = a class on
  // "artv-common/dsp/own/parts/interpolation/stateless.hpp"
  std::array<value_type, n_channels> get (
    std::array<float, n_channels> delay_spls)
  {
    return get<interp::n_points, interp::x_offset> (
      delay_spls, [] (auto y, auto x) { return interp::tick (y, x); });
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  template <uint N_points, uint X_offset, class InterpFunctor>
  std::array<value_type, n_channels> get (
    std::array<float, n_channels> delay_spls,
    InterpFunctor&&               interp)
  {
    // calculate delays
    std::array<uint, n_channels>  delay_spls_int;
    std::array<float, n_channels> delay_spls_frac;

    for (uint i = 0; i < n_channels; ++i) {
      delay_spls_int[i]  = (uint) delay_spls[i];
      delay_spls_frac[i] = delay_spls[i] - (float) delay_spls_int[i];
      delay_spls_int[i] -= X_offset;
    }
    // fetch all samples to interpolate
    array2d<value_type, N_points, n_channels> spls;
    for (uint c = 0; c < n_channels; ++c) {
      for (uint i = 0; i < N_points; ++i) {
        spls[c][i] = base::get (delay_spls_int[c] + i, c);
      }
    }
    // interpolate
    std::array<value_type, n_channels> ret;
    for (uint c = 0; c < n_channels; ++c) {
      auto frac = make_vec ((value_type) delay_spls_frac[c]);
      std::array<decltype (frac), N_points> p;
      for (uint i = 0; i < p.size(); ++i) {
        p[i] = make_vec (spls[c][i]);
      }
      ret[c] = interp (p, frac)[0]; // vectors of size 1
    }
    return ret;
  }
};
//------------------------------------------------------------------------------
// TODO: Thiran and others with the multichannel interface.
//------------------------------------------------------------------------------
} // namespace artv
