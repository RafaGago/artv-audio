#pragma once

#include <cmath>
#include <limits>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/oscillators/phasor.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

//------------------------------------------------------------------------------
template <uint N>
class lfo {
public:
  //----------------------------------------------------------------------------
  static constexpr uint n_channels = N;
  using value_type                 = vec<float, N>;
  using phase_type                 = phase<N>;
  //----------------------------------------------------------------------------
  void reset()
  {
    _phasor.reset();
    _smooth_coeff = vec_set<N> (0.f);
    _smooth_state = vec_set<N> (0.f);
  }
  //----------------------------------------------------------------------------
  void set_freq (value_type f, float t_spl)
  {
    _phasor.set_freq (f, t_spl);
    onepole_smoother::reset_coeffs (make_xspan (_smooth_coeff), f * 2.f, t_spl);
  }
  //----------------------------------------------------------------------------
  void set_phase (phase_type p) { _phasor.set_phase (p); }
  //----------------------------------------------------------------------------
  phase_type get_phase() const { return _phasor.get_phase(); }
  //----------------------------------------------------------------------------
  value_type tick_sine (uint n = 1)
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
    auto negative  = phase < 0;
    auto calcphase = negative ? phase : -phase;
    auto p         = vec_cast<double> (calcphase);
    auto dsin      = 1.862645149230957e-09 * p + 8.673617379884035e-19 * p * p;
    auto dsinf     = vec_cast<float> (dsin);
    return (negative ? dsinf : -dsinf);
  }
  //----------------------------------------------------------------------------
  value_type tick_saw (uint n = 1)
  {
    return _phasor.tick (n).to_normalized_bipolar();
  }
  //----------------------------------------------------------------------------
  value_type tick_square (uint n = 1)
  {
    return _phasor.tick (n).to_int() < 0 ? -1.f : 1.f;
  }
  //----------------------------------------------------------------------------
  value_type tick_triangle (uint n = 1)
  {
    // triangle by abs of phase (saw) + scaling. Adding a quarter cycle so its
    // starts at 0
    const auto quarter_cycle
      = phase<N> {typename phase<N>::degrees {}, vec_set<N> (90.f)};
    const auto quarter_range = quarter_cycle.to_uint();

    auto ph = _phasor.tick (n);
    ph += quarter_cycle;
    auto ph_int = ph.to_int();
    ph_int      = (ph_int < 0) ? -ph_int : ph_int;
    // make it bipolar
    ph_int -= quarter_range; // -INT32_MAX/2 to INT32_MAX/2

    static constexpr double step
      = 2. / (double) std::numeric_limits<typename phase<N>::raw_int>::max();

    return vec_cast<float> (vec_cast<double> (ph_int) * step);
  }
  //----------------------------------------------------------------------------
  // A trapezoid done by clipping a triangle and rescaling to -1 to 1.
  value_type tick_trapezoid (value_type clip_level, uint n = 1)
  {
    // triangle by abs of phase + scaling
    auto tri = tick_triangle (n);
    tri      = vec_min (clip_level, tri);
    tri      = vec_max (-clip_level, tri);
    return tri * (1. / clip_level);
  }
  //----------------------------------------------------------------------------
  value_type tick_sample_hold (uint n = 1)
  {
    auto new_cycle = _phasor.tick_ext (n).new_cycle;
    for (uint i = 0; i < N; ++i) {
      if (new_cycle[i]) {
        _specific_data[i] = _whitenoise (1.f);
      }
    }
    return _specific_data;
  }
  //----------------------------------------------------------------------------
  value_type tick_filt_sample_and_hold (uint n = 1)
  {
    value_type ret {};
    for (uint i = 0; i < n; ++i) {
      ret = tick_sample_hold (1);
      ret = onepole_smoother::tick<decltype (ret)> (
        make_xspan (_smooth_coeff), make_xspan (_smooth_state), ret);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  value_type tick_pm (phase<N> delta, uint n = 1)
  {
    auto old  = _phasor.get_increment();
    auto newv = old + delta;
    _phasor.set_increment (newv);
    auto ret = _phasor.tick_sine (n);
    _phasor.set_increment (old);
    return ret;
  }
  //----------------------------------------------------------------------------
private:
  int_phasor<N>         _phasor;
  value_type            _specific_data;
  white_noise_generator _whitenoise;
  static_assert (onepole_smoother::n_coeffs == 1, "");
  static_assert (onepole_smoother::n_coeffs_int == 0, "");
  value_type _smooth_coeff;
  static_assert (onepole_smoother::n_states == 1, "");
  value_type _smooth_state;
  //----------------------------------------------------------------------------
};

} // namespace artv
