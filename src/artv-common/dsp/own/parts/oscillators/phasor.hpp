#pragma once

#include <cassert>
#include <limits>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

// REMINDER: this is not constexpr because vector types aren't
//------------------------------------------------------------------------------
template <uint N = 1>
class phase {
public:
  struct degrees {};
  struct normalized {};
  struct bipolar_normalized {};

  static constexpr uint n_channels = N;

  using phase_uint = vec<u32, N>;
  using phase_int  = vec<s32, N>;

  using raw_uint = typename vec_traits_t<phase_uint>::value_type;
  using raw_int  = typename vec_traits_t<phase_int>::value_type;
  //----------------------------------------------------------------------------
  phase() : _v {0} {}
  //----------------------------------------------------------------------------
  explicit phase (phase_uint raw) : _v {raw} {}
  //----------------------------------------------------------------------------
  explicit phase (phase_int raw) : _v {(phase_uint) raw} {}
  //----------------------------------------------------------------------------
  phase (vec<float, N> deg, degrees) // 0 to 360
    : _v {from_degrees (deg)}
  {}
  //----------------------------------------------------------------------------
  phase (vec<float, N> norm, normalized) // 0 to 1
    : _v {from_normalized (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (vec<float, N> norm, bipolar_normalized) // -1 to 1
    : _v {from_bipolar_normalized (norm)}
  {}
  //----------------------------------------------------------------------------
  phase (phase const& other) = default;
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
  phase_uint to_uint() const { return _v; }
  //----------------------------------------------------------------------------
  phase_int to_int() const // -INT32_MAX to +INT32_MAX
  {
    return make_symetric_int ((phase_int) _v);
  }
  //----------------------------------------------------------------------------
  // 0 to 1,
  vec<float, N> to_normalized() const
  {
    constexpr auto pstep
      = 1. / ((double) std::numeric_limits<phase_uint>::max());
    vec<double, N> v = vec_cast<double> (_v);
    return vec_cast<float> (v * pstep);
  }
  //----------------------------------------------------------------------------
  // 0 to 360
  vec<float, N> to_degrees() const { to_normalized() * vec_set<N> (360.f); }
  //----------------------------------------------------------------------------
  // -1 to 1, mapping from 0 to 360deg
  vec<float, N> to_normalized_bipolar() const
  {
    constexpr auto pstep = 1. / ((double) std::numeric_limits<raw_int>::max());
    vec<double, N> v     = vec_cast<double> (to_int());
    return vec_cast<float> (v * pstep);
  }
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
  phase_uint from_degrees (vec<float, N> deg)
  {
    constexpr auto pstep
      = ((double) std::numeric_limits<raw_uint>::max()) / 360.;
    vec<double, N> deg_dbl = vec_cast<double> (deg);
    return vec_cast<raw_uint> (deg_dbl * pstep);
  }

  phase_uint from_normalized (vec<float, N> norm)
  {

    constexpr auto pstep = (double) std::numeric_limits<raw_uint>::max();
    return vec_cast<raw_uint> (vec_cast<double> (norm) * pstep);
  }

  phase_uint from_bipolar_normalized (vec<float, N> bnorm)
  {
    constexpr auto pstep = ((double) std::numeric_limits<raw_int>::max());
    return vec_cast<raw_int> (vec_cast<double> (bnorm) * pstep);
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
  void set_freq (vec<float, N> f, float samplerate)
  {
    constexpr auto pmax
      = (float) std::numeric_limits<typename phase<N>::raw_uint>::max();
    set_increment (phase<N> {vec_cast<u32> (pmax / samplerate * f)});
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
