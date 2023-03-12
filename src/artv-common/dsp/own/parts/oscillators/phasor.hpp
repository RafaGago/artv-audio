#pragma once

#include <array>
#include <cassert>
#include <limits>

#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

namespace phase_tag {
// these all map degrees from 0 to 360. So they start at 0 and end up at 0 (2pi)
struct degrees {}; // 0 to 360
struct normalized {}; // 0 to 1
struct normalized_bipolar {}; // -1 to 1

// these all map degrees from -180 to 180. So they start and end at pi (3pi)
struct shifted_dregrees {}; // range = -180 to 180
struct shifted {}; // 0 to 1
struct shifted_bipolar {}; // -1 to 1
} // namespace phase_tag

// REMINDER: this is not constexpr because vector types aren't
//------------------------------------------------------------------------------
template <uint N = 1>
class phase {
public:
  static constexpr uint n_channels = N;

  using degrees            = phase_tag::degrees;
  using normalized         = phase_tag::normalized;
  using normalized_bipolar = phase_tag::normalized_bipolar;

  using shifted         = phase_tag::shifted;
  using shifted_bipolar = phase_tag::shifted_bipolar;

  using phase_uint = vec<u32, N>;
  using phase_int  = vec<s32, N>;

  using float_vec = vec<float, N>;

  using scalar_uint = typename vec_traits_t<phase_uint>::value_type;
  using scalar_sint = typename vec_traits_t<phase_int>::value_type;

  template <uint S, uint I, uint F>
  using fixpt_type      = fixpt_s<S, I, F, 0, vec_fp_types_trait<N>>;
  using fixpt_type_uint = fixpt_type<0, 0, sizeof (scalar_uint) * 8>;
  using fixpt_type_int  = fixpt_type<1, 0, (sizeof (scalar_uint) * 8) - 1>;
  //----------------------------------------------------------------------------
  phase() : _v {0} {}
  //----------------------------------------------------------------------------
  explicit phase (phase_uint raw) : _v {raw} {}
  //----------------------------------------------------------------------------
  explicit phase (phase_int raw) : _v {(phase_uint) raw} {}
  //----------------------------------------------------------------------------
  phase (degrees, float_vec deg) // 0 to 360
    : _v {from_degrees (deg)}
  {}
  //----------------------------------------------------------------------------
  phase (normalized, float_vec norm) // 0 to 1 mapping from o to 360
    : _v {from_normalized (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (normalized_bipolar, float_vec norm) // -1 to 1 mapping from o to 360
    : _v {from_normalized_bipolar (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (shifted, float_vec norm) // 0 to 1 mapping from -180 to 180
    : _v {from_shifted (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (shifted_bipolar, float_vec norm) // -1 to 1 mapping from -180 to 180
    : _v {from_bipolar_shifteded (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (degrees, std::array<float, N> deg) // 0 to 360
    : _v {from_degrees (vec_from_array (deg))}
  {}
  //----------------------------------------------------------------------------
  phase (normalized, std::array<float, N> norm) // 0 to 1 mapping from o to 360
    : _v {from_normalized (vec_from_array (norm))}
  {}
  //----------------------------------------------------------------------------
  phase (normalized_bipolar, std::array<float, N> norm) // -1 to 1
    : _v {from_normalized_bipolar (vec_from_array (norm))}
  {}
  //----------------------------------------------------------------------------
  phase (shifted, std::array<float, N> norm) // 0 to 1 mapping from -180 to 180
    : _v {from_shifted (vec_from_array (norm))}
  {}
  //----------------------------------------------------------------------------
  phase (
    shifted_bipolar,
    std::array<float, N> norm) // -1 to 1 mapping from -180 to 180
    : _v {from_shifted_bipolar (vec_from_array (norm))}
  {}
  //----------------------------------------------------------------------------
  template <class... Ts>
  phase (degrees, Ts&&... deg) // 0 to 360
    : _v {from_degrees (float_vec {std::forward<Ts> (deg)...})}
  {
    static_assert (sizeof...(Ts) == N);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  phase (normalized, Ts&&... norm) // 0 to 1
    : _v {from_normalized (float_vec {std::forward<Ts> (norm)...})}
  {
    static_assert (sizeof...(Ts) == N);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  phase (normalized_bipolar, Ts&&... norm) // -1 to 1
    : _v {from_normalized_bipolar (float_vec {std::forward<Ts> (norm)...})}
  {
    static_assert (sizeof...(Ts) == N);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  phase (shifted, Ts&&... norm) // 0 to 1 mapping from -180 to 180
    : _v {from_shifted (float_vec {std::forward<Ts> (norm)...})}
  {
    static_assert (sizeof...(Ts) == N);
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  phase (shifted_bipolar, Ts&&... norm) // -1 to 1 mapping from -180 to 180
    : _v {from_shifted_bipolar (float_vec {std::forward<Ts> (norm)...})}
  {
    static_assert (sizeof...(Ts) == N);
  }
  //----------------------------------------------------------------------------
  phase (phase const& other)            = default;
  phase& operator= (phase const& other) = default;
  ~phase()                              = default;
  //----------------------------------------------------------------------------
  phase& operator+= (phase const& other)
  {
    _v += other._v;
    return *this;
  }
  friend phase operator+ (phase lhs, phase rhs)
  {
    lhs += rhs;
    return lhs;
  }
  phase& operator+= (phase_uint other)
  {
    _v += other;
    return *this;
  }
  friend phase operator+ (phase lhs, phase_uint rhs)
  {
    lhs += rhs;
    return lhs;
  }
  //----------------------------------------------------------------------------
  phase& operator-= (phase const& other)
  {
    _v -= other._v;
    return *this;
  }
  friend phase operator- (phase lhs, phase rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  phase& operator-= (phase_uint other)
  {
    _v -= other;
    return *this;
  }
  friend phase operator- (phase lhs, phase_uint rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  //----------------------------------------------------------------------------
  phase& operator*= (phase const& other)
  {
    _v *= other._v;
    return *this;
  }
  friend phase operator* (phase lhs, phase rhs)
  {
    lhs *= rhs;
    return lhs;
  }
  phase& operator*= (phase_uint v)
  {
    _v *= v;
    return *this;
  }
  friend phase operator* (phase lhs, phase_uint v)
  {
    lhs *= v;
    return lhs;
  }
  //----------------------------------------------------------------------------
  static phase_uint to_uint (phase_uint x) { return x; }
  phase_uint        to_uint() const { return to_uint (_v); }
  //----------------------------------------------------------------------------
  // discontinuty at half the range 0 to max, -max to 0.
  static phase_int to_int (phase_uint x) { return (phase_int) x; }
  phase_int        to_int() const // discontinuty at half the range
  {
    return to_int (_v);
  }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> 0 to 360
  static float_vec to_degrees (phase_uint x)
  {
    return vec_cast<float> (x) * (1.f / deg_step);
  }
  float_vec to_degrees() const { to_degrees (_v); }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> 0 to 1
  static float_vec to_normalized (phase_uint x)
  {
    return vec_cast<float> (x) * (1.f / range_step);
  }
  float_vec to_normalized() const { return to_normalized (_v); }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> -1 to 1
  static float_vec to_normalized_bipolar (phase_uint x)
  {
    return to_shifted_bipolar (x - half_range_uint);
  }
  float_vec to_normalized_bipolar() const { return to_normalized_bipolar (_v); }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> -180 to 180
  static float_vec to_shifted_degrees (phase_uint x)
  {
    return vec_cast<float> (vec_cast<scalar_sint> (x)) * (1.f / half_deg_step);
  }
  float_vec to_shifted_degrees() const { return to_shifted_degrees (_v); }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> [0.5 to 1] [0 to 0.5]
  static float_vec to_shifted (phase_uint x)
  {
    return to_normalized (x + half_range_uint);
  }
  float_vec to_shifted() const { return to_shifted (_v); }
  //----------------------------------------------------------------------------
  // 0 to 2pi -> [0 to 1] [-1 to 0]
  static float_vec to_shifted_bipolar (phase_uint x)
  {
    return vec_cast<float> (vec_cast<scalar_sint> (x))
      * (1.f / half_range_step);
  }
  float_vec to_shifted_bipolar() const { return to_shifted_bipolar (_v); }
  //----------------------------------------------------------------------------
  // o to 360
  fixpt_type_int to_fixpt() { return fixpt_type_int::from (_v); }
  //----------------------------------------------------------------------------
  // o to 180, -180 to 0
  fixpt_type_uint to_fixpt_unsigned() { return fixpt_type_uint::from (_v); }
  //----------------------------------------------------------------------------
  scalar_uint get_raw (uint channel) const
  {
    assert (channel < N);
    return _v[channel];
  }
  //----------------------------------------------------------------------------
  void set_raw (scalar_uint v, uint channel)
  {
    assert (channel < N);
    _v[channel] = v;
  }
  //----------------------------------------------------------------------------
  phase_uint get_raw() const { return _v; }
  //----------------------------------------------------------------------------
  void set_raw (phase_uint v) const { return _v = v; }
  //----------------------------------------------------------------------------
  phase<N> set_at_uniform_spacing() { return set_at_uniform_spacing_from (0); }
  //----------------------------------------------------------------------------
  phase<N> set_at_uniform_spacing_from (uint reference_index)
  {
    constexpr scalar_uint step = std::numeric_limits<scalar_uint>::max() / N;
    assert (reference_index < N);
    auto old = _v;
    for (uint i = 0; i < N; ++i) {
      _v[i] = i * step;
    }
    _v += old[reference_index];
    return *this;
  }
  //----------------------------------------------------------------------------
private:
  static constexpr auto deg_step
    = (float) (((double) std::numeric_limits<scalar_uint>::max()) / 360.);

  static constexpr auto half_deg_step
    = (float) (((double) std::numeric_limits<scalar_uint>::max()) / 180.);

  static constexpr auto range_step
    = (float) std::numeric_limits<scalar_uint>::max();

  static constexpr auto half_range_uint
    = std::numeric_limits<scalar_uint>::max();

  static constexpr auto half_range_step = (float) (half_range_uint / 2);
  //----------------------------------------------------------------------------
  phase_uint from_degrees (float_vec v) // 0 to 360 -> 0 to 2pi
  {
    return vec_cast<scalar_uint> (v * deg_step);
  }
  //----------------------------------------------------------------------------
  phase_uint from_normalized (float_vec v) // 0 to 1 -> 0 to 2pi
  {
    return vec_cast<scalar_uint> (v * range_step);
  }
  //----------------------------------------------------------------------------
  phase_uint from_normalized_bipolar (float_vec v) // -1 to 1 -> 0 to 2pi
  {
    return from_shifted_bipolar (v) - half_range_uint;
  }
  //----------------------------------------------------------------------------
  phase_uint from_shifted_degrees (float_vec v) // -180 to 180 -> 1 to 3pi
  {
    return vec_cast<scalar_uint> (vec_cast<scalar_sint> (v * half_deg_step));
  }
  //----------------------------------------------------------------------------
  phase_uint from_shifted (float_vec v) // 0 to 1 -> 1 to 3pi
  {
    return from_normalized (v) + half_range_uint;
  }
  //----------------------------------------------------------------------------
  phase_uint from_shifted_bipolar (float_vec v) // -1 to 1 -> 1 to 3pi
  {
    return vec_cast<scalar_uint> (vec_cast<scalar_sint> (v * half_range_step));
  }
  //----------------------------------------------------------------------------
  phase_uint _v; // 0 to 2pi
};
//------------------------------------------------------------------------------
template <uint N = 1>
class int_phasor {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = N;
  //----------------------------------------------------------------------------
  void reset()
  {
    _phase     = phase<N> {vec_set<N> (0u)};
    _increment = phase<N> {vec_set<N> (0u)};
  }
  //----------------------------------------------------------------------------
  void set_increment (phase<N> p) { _increment = p; }
  //----------------------------------------------------------------------------
  void set_freq (vec<float, N> f, float t_spl)
  {
    constexpr auto pmax
      = (float) std::numeric_limits<typename phase<N>::scalar_uint>::max();
    set_increment (phase<N> {vec_cast<u32> (pmax * t_spl * f)});
  }
  //----------------------------------------------------------------------------
  void     set_phase (phase<N> p) { _phase = p; }
  phase<N> get_phase() const { return _phase; }
  phase<N> get_increment() const { return _increment; }
  //----------------------------------------------------------------------------
  phase<N> tick()
  {
    _phase += _increment;
    return _phase;
  }

  phase<N> tick (uint n)
  {
    using ph_uint = typename phase<N>::phase_uint;
    _phase += phase<N> (_increment) * vec_set<ph_uint> (n);
    return _phase;
  }
  //----------------------------------------------------------------------------
  struct tick_ext_ret {
    phase<N> ph;
    typename vec_traits_t<typename phase<N>::phase_uint>::same_size_int_type
      new_cycle;
  };
  // As tick, but sends a flag indicating a new cycle
  tick_ext_ret tick_ext (uint n = 1)
  {
    for (uint i = 0; i < N; ++i) {
      assert (
        (((double) _increment.to_uint()[i]) * n) < ((
          double) std::numeric_limits<typename phase<N>::scalar_uint>::max())
        && "Phase increment too big for new cycle detection to work reliably");
    }

    tick_ext_ret ret;
    auto         phaseprev = _phase.to_uint();
    ret.ph                 = tick (n);
    ret.new_cycle          = phaseprev > _phase.to_uint();
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  phase<N> _increment;
  phase<N> _phase;
};
//------------------------------------------------------------------------------
} // namespace artv
