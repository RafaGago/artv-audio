#pragma once

#include <array>
#include <cassert>
#include <limits>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

namespace phase_tag {
// these all map degrees from 0 to 360. So they start at 0 and end up at 0 (2pi)
struct degrees {}; // 0 to 360
struct normalized {}; // 0 to 1
struct normalized_bipolar {}; // -1 to 1

// these all map degrees from -180 to 180. So they start and end at pi (3pi)
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

  using raw_uint = typename vec_traits_t<phase_uint>::value_type;
  using raw_int  = typename vec_traits_t<phase_int>::value_type;
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
  static phase_int to_int (phase_uint x) // -INT32_MAX to +INT32_MAX
  {
    return make_symetric_int ((phase_int) x);
  }
  phase_int to_int() const // -INT32_MAX to +INT32_MAX
  {
    return to_int (_v);
  }
  //----------------------------------------------------------------------------
  // 0 to 360
  static float_vec to_degrees (phase_uint x)
  {
    to_normalized (x) * vec_set<N> (360.f);
  }
  float_vec to_degrees() const { to_degrees (_v); }
  //----------------------------------------------------------------------------
  // 0 to 1,
  static float_vec to_normalized (phase_uint x)
  {
    constexpr auto pstep
      = 1. / ((double) std::numeric_limits<phase_uint>::max());
    vec<double, N> v = vec_cast<double> (x);
    return vec_cast<float> (v * pstep);
  }
  float_vec to_normalized() const { return to_normalized (_v); }
  //----------------------------------------------------------------------------
  // -1 to 1, mapping from 0 to 360deg
  static float_vec to_normalized_bipolar (phase_uint x)
  {
    constexpr auto pstep
      = 1. / ((double) std::numeric_limits<raw_uint>::max() * 2.);
    vec<double, N> v = vec_cast<double> (x);
    return vec_cast<float> (v * pstep) - 1.f;
  }
  float_vec to_normalized_bipolar() const { return to_normalized_bipolar (_v); }
  //----------------------------------------------------------------------------
  // 0 to 1, mapping from -180 to 180deg. 0 degrees at 0.5.
  static float_vec to_shifted (phase_uint x)
  {
    constexpr auto half = std::numeric_limits<raw_uint>::max() / 2;
    return to_normalized (x + half);
  }
  float_vec to_shifted() const { return to_shifted (_v); }
  //----------------------------------------------------------------------------
  // 0 to 1, mapping from -180 to 180deg. 0 degrees at 0.5.
  static float_vec to_shifted_bipolar (phase_uint x)
  {
    constexpr auto half = std::numeric_limits<raw_uint>::max() / 2;
    return to_normalized_bipolar (x + half);
  }
  float_vec to_shifted_bipolar() const { return to_shifted_bipolar (_v); }
  //----------------------------------------------------------------------------
  static phase_int make_symetric_int (phase_int v)
  {
    // avoid UB, this value has no positive equivalent on a signed int, but we
    // still want phases to have a symmetric scaling on the outside. This
    // introduces an extremely small non accumulative error.
    v += (v == std::numeric_limits<phase_int>::min());
    return v;
  }
  //----------------------------------------------------------------------------
  raw_uint get_raw (uint channel) const
  {
    assert (channel < N);
    return _v[channel];
  }
  //----------------------------------------------------------------------------
  void set_raw (raw_uint v, uint channel)
  {
    assert (channel < N);
    _v[channel] = v;
  }
  //----------------------------------------------------------------------------
private:
  // C++ **** for being able to make constexpr constructors on an unevaluated
  // context (TODO: this wasn't using vector types previously. It seems that
  // vector types break constexprness, if at some point on the future they
  // don't, then add constexprness back).
  phase_uint from_degrees (float_vec deg)
  {
    constexpr auto pstep
      = ((double) std::numeric_limits<raw_uint>::max()) / 360.;
    vec<double, N> deg_dbl = vec_cast<double> (deg);
    return vec_cast<raw_uint> (deg_dbl * pstep);
  }

  phase_uint from_normalized (float_vec norm)
  {
    constexpr auto pstep = (double) std::numeric_limits<raw_uint>::max();
    return vec_cast<raw_uint> (vec_cast<double> (norm) * pstep);
  }

  phase_uint from_normalized_bipolar (float_vec bnorm)
  {
    constexpr auto pstep
      = ((double) (std::numeric_limits<raw_uint>::max() / 2));
    return vec_cast<raw_uint> (vec_cast<double> (bnorm + 1.) * pstep);
  }

  phase_uint from_shifted (float_vec norm)
  {
    constexpr auto half = std::numeric_limits<raw_uint>::max() / 2;
    phase_uint     n    = from_normalized (norm);
    return n - half;
  }

  phase_uint from_shifted_bipolar (float_vec bnorm)
  {
    constexpr auto half = std::numeric_limits<raw_uint>::max() / 2;
    phase_uint     n    = from_normalized_bipolar (bnorm);
    return n - half;
  }
  //----------------------------------------------------------------------------
  phase_uint _v;
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
      = (float) std::numeric_limits<typename phase<N>::raw_uint>::max();
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
        (((double) _increment.to_uint()[i]) * n)
          < ((double) std::numeric_limits<typename phase<N>::raw_uint>::max())
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
