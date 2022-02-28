#pragma once

#include <algorithm>
#include <cassert>

#include <gcem.hpp>

#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

#include "artv-common/dsp/own/parts/filters/onepole.hpp"
#include "artv-common/dsp/own/parts/misc/slew_limiter.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"

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
    n_det_modes
  };
  //----------------------------------------------------------------------------
  void set_detector_hipass (float hz)
  {
    _in_filter.reset_coeffs (vec_set<vec_type> (hz), _sample_rate);
  }
  //----------------------------------------------------------------------------
  void set_detector_recovery (float ratio)
  {
    assert (ratio >= 0. && ratio <= 1.);

    static constexpr double minv = gcem::log (1. / 10.);
    static constexpr double maxv = gcem::log (1.) - minv;

    double release_sec = exp (maxv * ratio + minv);

    _followers.reset_coeffs<flw_slow> (
      vec_set<vec_type> (attack_sec),
      vec_set<vec_type> (release_sec),
      _sample_rate);

    _followers.reset_coeffs<flw_fast> (
      vec_set<vec_type> (attack_sec * 0.1),
      vec_set<vec_type> (release_sec),
      _sample_rate);
  }
  //----------------------------------------------------------------------------
  void set_detector_channels (uint m)
  {
    assert (m < n_det_modes);
    _det_mode = (detector_mode) m;
  }
  //----------------------------------------------------------------------------
  void set_curve_decay (float ratio)
  {
    assert (ratio >= 0. && ratio <= 1.);

    static constexpr double minv = gcem::log (1. / 1000.);
    static constexpr double maxv = gcem::log (2.5) - minv;

    double release_sec = exp (maxv * ratio + minv);

    _followers.reset_coeffs<flw_gate> (
      vec_set<vec_type> (1. / 2000.),
      vec_set<vec_type> (release_sec),
      _sample_rate);
  }
  //----------------------------------------------------------------------------
  void set_curve_decay_shape (uint v)
  {
    assert (v <= 5);
    _curve_shape = v;
  }
  //----------------------------------------------------------------------------
  void set_detector_shape (uint v)
  {
    assert (v <= 5);
    _det_shape = v;
  }
  //----------------------------------------------------------------------------
  void reset (uint sample_rate)
  {
    _sample_rate = sample_rate;
    _det_mode    = detector_stereo;
    _det_shape   = 0;
    _curve_shape = 0;

    _followers.reset_states_cascade();
    _in_filter.reset_states();
    _gate_filter.reset_states();
    _gate_filter.reset_coeffs (vec_set<vec_type> (250.), _sample_rate);
  };
  //----------------------------------------------------------------------------
  vec_type tick (vec_type in)
  {
    vec_type spl = _in_filter.tick (in);

    switch (_det_mode) {
    case detector_stereo:
      spl = vec_abs (spl);
      break;
    case detector_mid:
      spl[0] = (std::abs (spl[0]) + std::abs (spl[1])) * 0.5;
      spl[1] = spl[0];
      break;
    case detector_l:
      spl[1] = spl[0];
      spl    = vec_abs (spl);
      break;
    case detector_r:
      spl[0] = spl[1];
      spl    = vec_abs (spl);
      break;
    default:
      assert (false);
    }

    spl = apply_shape (spl, _det_shape);
    // logaritmic space envelopes
    spl           = vec_log (spl + 0.00000000000000001);
    vec_type fast = vec_exp (_followers.tick<flw_fast> (spl));
    vec_type slow = vec_exp (_followers.tick<flw_slow> (spl));

    vec_type gate = (fast - slow) / fast; // A ratio. 0 to 1

    gate = _followers.tick<flw_gate> (gate);
    gate = _gate_filter.tick (gate);
    gate = apply_shape (gate, _curve_shape);

    return in * gate;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  vec_type apply_shape (vec_type in, uint mode)
  {
    // TODO: ADAA for the aliasing?
    switch (mode) {
    case 0:
      break;
    case 1:
      in *= in;
      break;
    case 2:
      in *= in * in;
      break;
    case 3:
      in *= in * in * in;
      break;
    case 4:
      in *= in * in * in * in;
      break;
    default:
      assert (false);
    }
    return in;
  }
  //----------------------------------------------------------------------------
  static constexpr double attack_sec = 1. / 200.; // 200Hz seconds

  enum followers { flw_gate, flw_fast, flw_slow, n_flw };
  part_class_array<slew_limiter, vec_type, n_flw> _followers;
  part_class_array<onepole_highpass, vec_type>    _in_filter;
  part_class_array<onepole_lowpass, vec_type>     _gate_filter;

  double _sample_rate;
  uint   _det_shape;
  uint   _curve_shape;

  detector_mode _det_mode;
};
//------------------------------------------------------------------------------
}; // namespace artv
