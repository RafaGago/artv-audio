#pragma once

#include <cassert>

#include <gcem.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/moving_average.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"

namespace artv {

//------------------------------------------------------------------------------
class transient_gate {
private:
  using vec_type     = double_x2;
  using vec_cmp_type = vec<s64, 2>;

public:
  //----------------------------------------------------------------------------
  enum detector_mode {
    detector_stereo,
    detector_mid,
    detector_l,
    detector_r,
    n_channel_modes
  };
  //----------------------------------------------------------------------------
  void set_channel_mode (uint m)
  {
    assert (m < n_channel_modes);
    _channel_mode = (detector_mode) m;
  }
  //----------------------------------------------------------------------------
  void set_detector_threshold (float db)
  {
    auto hyst     = _hi_threshold - _lo_threshold;
    _hi_threshold = db;
    set_detector_hysteresis (hyst);
  }
  //----------------------------------------------------------------------------
  void set_detector_hysteresis (float db)
  {
    _lo_threshold = _hi_threshold - db;
  }
  //----------------------------------------------------------------------------
  void set_detector_speed (float ratio)
  {
    assert (ratio >= 0. && ratio <= 1.);
    static constexpr double beta  = gcem::log (80);
    static constexpr double alpha = gcem::log (180.) - beta;

    auto speed = exp (alpha * ratio + beta);

    _followers.reset_coeffs<flw_slow> (
      vec_set<vec_type> (speed * 0.001),
      vec_set<vec_type> (0.2 * 0.001),
      _sample_rate);
  }
  //----------------------------------------------------------------------------
  void set_decay (float ratio)
  {
    assert (ratio >= 0. && ratio <= 1.);
    static constexpr double beta  = gcem::log (40);
    static constexpr double alpha = gcem::log (320.) - beta;

    auto decay = exp (alpha * ratio + beta);
    // the sample rate is divided by 2 because this was prototyped on JSFX using
    // another envelope code. TODO: fix on the constants above.
    _followers.reset_coeffs<flw_gate> (
      vec_set<vec_type> (0.001 * 0.001),
      vec_set<vec_type> (decay * 0.001),
      _sample_rate);
  }
  //----------------------------------------------------------------------------
  void set_decay_curve (uint v)
  {
    assert (v <= 5);
    _decay_curve = v;
  }
  //----------------------------------------------------------------------------
  void set_detector_curve (uint v)
  {
    assert (v <= 5);
    _detector_curve = v;
  }
  //----------------------------------------------------------------------------
  void reset (uint sample_rate)
  {
    _sample_rate    = sample_rate;
    _hi_threshold   = 0.;
    _lo_threshold   = 0.;
    _channel_mode   = detector_stereo;
    _decay_curve    = 2;
    _detector_curve = 2;
    _followers.reset_states_cascade();
    _in_filter.reset_states();
    // the sample rate is divided by 2 because this was prototyped on JSFX using
    // another envelope code. TODO: fix on the constants above.
    _followers.reset_coeffs<flw_fast> (
      vec_set<vec_type> (1. * 0.001),
      vec_set<vec_type> (1. * 0.001),
      _sample_rate);
  };
  //----------------------------------------------------------------------------
  vec_type tick (vec_type in)
  {
    // filter
    vec_type lp = _in_filter.tick (in);
    // not strictly a correct hipass without an input delay, as the boxcar has
    // latency. This is still enough for the purposes of a detector.
    vec_type hp = in - lp;

    // detector channel/mode select
    switch (_channel_mode) {
    case detector_stereo:
      hp = vec_abs (hp);
      break;
    case detector_mid:
      hp[0] = (std::abs (hp[0]) + std::abs (hp[1])) * 0.5;
      hp[1] = hp[0];
      break;
    case detector_l:
      hp[1] = hp[0];
      hp    = vec_abs (hp);
      break;
    case detector_r:
      hp[0] = hp[1];
      hp    = vec_abs (hp);
      break;
    default:
      assert (false);
    }

    hp = apply_curve (hp, _detector_curve); // tons of aliasing on the detector?

    // dB space detection.
    vec_type gain = 20. * vec_log10 (vec_max (vec_set<vec_type> (0.001), hp));
    vec_type fast = _followers.tick<flw_fast> (gain);
    vec_type slow = _followers.tick<flw_slow> (gain);
    vec_type diff = fast - slow;

    auto over  = diff > vec_set<vec_type> (_hi_threshold);
    auto under = diff < vec_set<vec_type> (_lo_threshold);
    _detect_on = _detect_on && !under ? !vec_set<vec_cmp_type> (0) : over;
    vec_type gate
      = _detect_on ? vec_set<vec_type> (1.) : vec_set<vec_type> (0.);

    // Smooth the transient gate signal with a third follower
    gate = _followers.tick<flw_gate> (gate);
    gate = apply_curve (gate, _decay_curve);

    return in * gate;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  vec_type apply_curve (vec_type in, uint mode)
  {
    switch (mode) {
    case 0:
      in = vec_cbrt (in);
      break;
    case 1:
      in = vec_sqrt (in);
      break;
    case 2:
      break;
    case 3:
      in *= in;
      break;
    case 4:
      in *= in * in;
      break;
    default:
      assert (false);
    }
    return in;
  }
  //----------------------------------------------------------------------------
  enum followers { flw_gate, flw_fast, flw_slow, n_flw };
  part_to_class<vec_type, slew_limiter, n_flw> _followers;
  part_to_class<vec_type, moving_average<4>>   _in_filter;

  double        _sample_rate;
  double        _hi_threshold;
  double        _lo_threshold;
  detector_mode _channel_mode;
  uint          _decay_curve;
  uint          _detector_curve;
  vec_cmp_type  _detect_on {};
};
//------------------------------------------------------------------------------
}; // namespace artv
