#pragma once

#include <cassert>
#include <limits>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

//------------------------------------------------------------------------------
class phase {
public:
  struct degrees {};
  struct normalized {};
  struct bipolar_normalized {};

  using phase_uint = u32;
  using phase_int  = s32;
  //----------------------------------------------------------------------------
  constexpr phase() : _v {0} {}
  //----------------------------------------------------------------------------
  constexpr explicit phase (phase_uint raw) : _v {raw} {}
  //----------------------------------------------------------------------------
  constexpr explicit phase (phase_int raw) : _v {(phase_uint) raw} {}
  //----------------------------------------------------------------------------
  constexpr phase (float deg, degrees) // 0 to 360
    : _v {from_degrees (deg)}
  {}
  //----------------------------------------------------------------------------
  constexpr phase (float norm, normalized) // 0 to 1
    : _v {from_normalized (norm)}
  {}
  //----------------------------------------------------------------------------
  constexpr phase (float norm, bipolar_normalized) // -1 to 1
    : _v {from_bipolar_normalized (norm)}
  {}
  //----------------------------------------------------------------------------
  constexpr phase (phase const& other) = default;
  constexpr phase& operator= (phase const& other) = default;
  ~phase()                                        = default;
  //----------------------------------------------------------------------------
  constexpr phase& operator+= (phase const& other)
  {
    _v += other._v;
    return *this;
  }
  friend constexpr phase operator+ (phase lhs, phase rhs)
  {
    lhs += rhs;
    return lhs;
  }
  //----------------------------------------------------------------------------
  constexpr phase& operator-= (phase const& other)
  {
    _v -= other._v;
    return *this;
  }
  friend constexpr phase operator- (phase lhs, phase rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  //----------------------------------------------------------------------------
  constexpr phase& operator*= (phase_uint v)
  {
    _v *= v;
    return *this;
  }
  friend constexpr phase operator* (phase lhs, phase_uint v)
  {
    lhs *= v;
    return lhs;
  }
  //----------------------------------------------------------------------------
  constexpr phase_uint to_uint() const { return _v; }
  //----------------------------------------------------------------------------
  constexpr phase_int to_int() const // -INT32_MAX to +INT32_MAX
  {
    return make_symetric_int ((phase_int) _v);
  }
  //----------------------------------------------------------------------------
  constexpr float to_degrees() const // 0 to 360
  {
    constexpr auto pstep
      = 360. / ((double) std::numeric_limits<phase_uint>::max());
    return (float) ((double) _v) * pstep;
  }
  //----------------------------------------------------------------------------
  constexpr float to_normalized() const // 0 to 1, mapping from 0 to 360deg
  {
    constexpr auto pstep
      = 1. / ((double) std::numeric_limits<phase_uint>::max());
    return (float) ((double) _v) * pstep;
  }
  //----------------------------------------------------------------------------
  constexpr float to_normalized_bipolar()
    const // -1 to 1, mapping from 0 to 360deg
  {
    constexpr auto pstep
      = 1. / ((double) std::numeric_limits<phase_int>::max());
    return (float) ((double) to_int()) * pstep;
  }
  //----------------------------------------------------------------------------
  constexpr static phase_int make_symetric_int (phase_int v)
  {
    // avoid UB, this value has no positive equivalent on a signed int, but we
    // still want phases to have a symmetric scaling on the outside. This
    // introduces an extremely small non accumulative error.
    v += (v == std::numeric_limits<phase_int>::min());
    return v;
  }
  //----------------------------------------------------------------------------
private:
  // C++ **** for being able to make constexpr constructors on an unevaluated
  // context.
  constexpr phase_uint from_degrees (float deg)
  {
    constexpr auto pstep
      = ((double) std::numeric_limits<phase_uint>::max()) / 360.;
    return (phase_uint) (deg * pstep);
  }

  constexpr phase_uint from_normalized (float norm)
  {
    constexpr auto pstep = (double) std::numeric_limits<phase_uint>::max();
    return (phase_uint) (norm * pstep);
  }

  constexpr phase_uint from_bipolar_normalized (float bnorm)
  {
    constexpr auto pstep = ((double) std::numeric_limits<phase_int>::max());
    return (phase_int) (bnorm * pstep);
  }
  //----------------------------------------------------------------------------
  phase_uint _v;
};

//------------------------------------------------------------------------------
class int_phasor {
public:
  //----------------------------------------------------------------------------
  void reset()
  {
    _phase     = phase {0u};
    _increment = phase {0u};
  }
  //----------------------------------------------------------------------------
  void set_freq (float f, float samplerate)
  {
    constexpr auto pmax = (float) std::numeric_limits<phase::phase_uint>::max();
    _increment          = phase {(u32) (pmax / samplerate * f)};
  }
  //----------------------------------------------------------------------------
  void  set_phase (phase p) { _phase = p; }
  phase get_phase() const { return _phase; }
  //----------------------------------------------------------------------------
  phase tick()
  {
    _phase += _increment;
    return _phase;
  }

  phase tick (uint n)
  {
    _phase += phase (_increment) * n;
    return _phase;
  }
  //----------------------------------------------------------------------------
  struct tick_ext_ret {
    phase ph;
    bool  new_cycle;
  };

  tick_ext_ret tick_ext (uint n = 1)
  {
    assert (
      (((double) _increment.to_uint()) * n)
        < ((double) std::numeric_limits<phase::phase_uint>::max())
      && "Phase increment too big for new cycle detection to work reliably");

    tick_ext_ret ret;
    auto         phaseprev = _phase.to_uint();
    ret.ph                 = tick (n);
    ret.new_cycle          = phaseprev > _phase.to_uint();
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  phase _increment;
  phase _phase;
};
//------------------------------------------------------------------------------
} // namespace artv
