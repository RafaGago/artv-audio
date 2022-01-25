#pragma once

#include <cmath>
#include <limits>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/phasor.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

//------------------------------------------------------------------------------
class lfo {
public:
  //----------------------------------------------------------------------------
  void reset()
  {
    _phasor.reset();
    _smooth_coeff = 0.f;
    _smooth_state = 0.f;
  }
  //----------------------------------------------------------------------------
  void set_freq (float f, float samplerate)
  {
    _phasor.set_freq (f, samplerate);
    onepole_smoother::reset_coeffs (
      make_crange (_smooth_coeff).cast (vec<float, 1> {}),
      make_vec (f * 2),
      samplerate);
  }
  //----------------------------------------------------------------------------
  void set_phase (phase p) { _phasor.set_phase (p); }
  //----------------------------------------------------------------------------
  phase get_phase() const { return _phasor.get_phase(); }
  //----------------------------------------------------------------------------
  float tick_sine (uint n = 1)
  {
    // parabola from 0 to -int32_min. Low precision, but enough for an LFO.
    // solving for points. basically this:
    // https://web.archive.org/web/20180919204945/http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/
    // but ranged from the phasor integer values to save the PI scaling
    //
    // Solving for points
    // x = 0, y = 0
    // x = -1073741824, y = -1
    // x = -2147483648, y = 0
    //
    // a = 1.862645149230957e-09
    // b = 8.673617379884035e-19
    auto phase     = _phasor.tick (n).to_int(); // -INT_MAX to INT_MAX
    bool negative  = phase < 0;
    auto calcphase = negative ? phase : -phase;
    auto p         = (double) calcphase;
    auto dsin      = 1.862645149230957e-09 * p + 8.673617379884035e-19 * p * p;
    return (float) negative ? dsin : -dsin;
  }
  //----------------------------------------------------------------------------
  float tick_saw (uint n = 1)
  {
    return _phasor.tick (n).to_normalized_bipolar();
  }
  //----------------------------------------------------------------------------
  float tick_square (uint n = 1)
  {
    return _phasor.tick (n).to_int() < 0 ? -1.f : 1.f;
  }
  //----------------------------------------------------------------------------
  float tick_triangle (uint n = 1)
  {
    // triangle by abs of phase (saw) + scaling. Adding a quarter cycle so its
    // starts at 0
    constexpr auto quarter_cycle = phase {90., phase::degrees {}};
    constexpr auto quarter_range = quarter_cycle.to_uint();

    auto ph = _phasor.tick (n);
    ph += quarter_cycle;
    auto ph_int = ph.to_int();
    ph_int      = (ph_int < 0) ? -ph_int : ph_int;
    // make it bipolar
    ph_int -= quarter_range; // -INT32_MAX/2 to INT32_MAX/2

    return ((double) ph_int)
      * (2. / std::numeric_limits<phase::phase_int>::max());
  }
  //----------------------------------------------------------------------------
  // A trapezoid done by clipping a triangle and rescaling to -1 to 1.
  float tick_trapezoid (float clip_level, uint n = 1)
  {
    // triangle by abs of phase + scaling
    float tri = tick_triangle (n);
    tri       = std::min (clip_level, tri);
    tri       = std::max (-clip_level, tri);
    return tri * (1. / clip_level);
  }
  //----------------------------------------------------------------------------
  float tick_sample_hold (uint n = 1)
  {
    // triangle by abs of phase + scaling
    if (_phasor.tick_ext (n).new_cycle) {
      _noise_value = _whitenoise (1.f);
    }
    return _noise_value;
  }
  //----------------------------------------------------------------------------
  float tick_filtered_noise (uint n = 1)
  {
    // triangle by abs of phase + scaling
    if (_phasor.tick_ext (n).new_cycle) {
      _noise_value = _whitenoise (1.f);
    }
    float ret;
    for (uint i = 0; i < n; ++i) {
      ret = onepole_smoother::tick (
        make_crange (_smooth_coeff),
        make_crange (_smooth_state).cast (vec<float, 1> {}),
        make_vec (_noise_value))[0];
    }
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  int_phasor            _phasor;
  float                 _noise_value;
  white_noise_generator _whitenoise;
  static_assert (onepole_smoother::n_coeffs == 1, "");
  static_assert (onepole_smoother::n_coeffs_int == 0, "");
  float _smooth_coeff;
  static_assert (onepole_smoother::n_states == 1, "");
  float _smooth_state;
  //----------------------------------------------------------------------------
};

} // namespace artv
