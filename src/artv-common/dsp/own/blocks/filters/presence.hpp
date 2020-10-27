#pragma once

// Port of Presence e.q by liteon

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

//------------------------------------------------------------------------------
namespace artv { namespace liteon {
// Presence (moorer) JSFX
struct presence_high_shelf {
  //----------------------------------------------------------------------------
  enum coeffs { a0, a1, a2, b1, b2, n_coeffs };
  enum state { y1, y2, x1, x2, n_states };
  //----------------------------------------------------------------------------
  // BW on the original JSFX is unitless BW from 0.007 to 0.4 and the lower end
  // makes it narrower (?). I scale it from 0 to 1.
  static void high_shelf (
    crange<double> co,
    double         freq,
    double         bogus_q,
    double         gain_db,
    double         sr)
  {
    double a        = 0.;
    double a2plus1  = 0.;
    double alphad   = 0.;
    double alphan   = 0.;
    double as2      = 0.;
    double as4      = 0.;
    double asnd     = 0.;
    double b0       = 0.;
    double boost    = 0.;
    double bw       = 0.;
    double c        = 0.;
    double ca       = 0.;
    double cf       = 0.;
    double cs       = 0.;
    double d        = 0.;
    double delta    = 0.;
    double f        = 0.;
    double f2       = 0.;
    double ma2plus1 = 0.;
    double mag      = 0.;
    double recipb0  = 0.;
    double sn       = 0.;
    double t        = 0.;
    double theta    = 0.;
    double tmp      = 0.;
    double xfmbw    = 0.;

    freq = std::max (freq, 3100.);
    freq = std::min (freq, 18500.);

    cf    = freq / sr;
    boost = gain_db;

    // scaling to the orignal range
    assert (bogus_q <= 1.);
    bw = ((1. - bogus_q) * 0.33) + 0.07;

    ca = std::tan (3.141592653589793 * (cf - 0.25));
    a  = pow (10., (boost / 20.));
    if ((boost < 6.0) && (boost > -6.0)) {
      (f = std::sqrt (a));
    }
    else {
      if (a > 1.) {
        (f = a / std::sqrt (2.));
      }
      else {
        (f = a * std::sqrt (2.));
      }
    }
    t     = std::tan (2. * 3.141592653589793 * bw);
    as2   = ca * ca;
    as4   = as2 * as2;
    d     = 2. * as2 * t;
    sn    = (1. + as4) * t;
    cs    = (1. - as4);
    mag   = std::sqrt (sn * sn + cs * cs);
    d     = mag;
    delta = std::atan2 (sn, cs);
    asnd  = (d <= 1. && d >= -1.) ? asin (d) : 0;
    theta = 0.5 * (3.141592653589793 - asnd - delta);
    tmp   = 0.5 * (asnd - delta);
    if (tmp > 0. && tmp < theta) {
      (theta = tmp);
    }
    xfmbw = theta / (2. * 3.141592653589793);
    c     = 1. / std::tan (2. * 3.141592653589793 * xfmbw);
    f2    = f * f;
    tmp   = a * a - f2;
#if 0
    if (std::fabs (tmp) <= spn) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#else
    if (std::fabs (tmp) <= 0.) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#endif
    alphan   = a * alphad;
    a2plus1  = 1. + as2;
    ma2plus1 = 1. - as2;
    co[a0]   = a2plus1 + alphan * ma2plus1;
    co[a1]   = 4.0 * ca;
    co[a2]   = a2plus1 - alphan * ma2plus1;
    b0       = a2plus1 + alphad * ma2plus1;
    co[b2]   = a2plus1 - alphad * ma2plus1;
    recipb0  = 1. / b0;
    co[a0] *= recipb0;
    co[a1] *= recipb0;
    co[a2] *= recipb0;
    co[b1] = co[a1];
    co[b2] *= recipb0;
    co[b1] = -co[b1];
    co[b2] = -co[b2];
    ;
  }
  //----------------------------------------------------------------------------
  static void repair_unsmoothable_coeffs (crange<double>, crange<const double>)
  {}
  //----------------------------------------------------------------------------
  static void reset_states (crange<double> s)
  {
    assert (s.size() >= n_states);
    memset (s.data(), 0, sizeof s[0] * n_states);
  }
  //----------------------------------------------------------------------------
  static double tick (
    crange<const double> co, // coeffs
    crange<double>       st, // state
    double               in)
  {
    assert (st.size() >= n_states);
    assert (co.size() >= n_coeffs);

    double out = co[a0] * in + co[a1] * st[x1] + co[a2] * st[x2]
      + co[b1] * st[y1] + co[b2] * st[y2];
    st[x2] = st[x1];
    st[x1] = in;
    st[y2] = st[y1];
    st[y1] = out;
    return out;
  }
  //----------------------------------------------------------------------------
  static simd_dbl tick (
    crange<const double>          co, // coeffs
    std::array<crange<double>, 2> st, // state
    std::array<double, 2>         in)
  {
    assert (st.size() >= 2);
    assert (co.size() >= n_coeffs);
    assert (st[0].size() >= n_states);
    assert (st[1].size() >= n_states);

    simd_dbl a1v {co[a1]};
    simd_dbl a2v {co[a2]};
    simd_dbl b1v {co[b1]};
    simd_dbl b2v {co[b2]};

    simd_dbl x1v {st[0][x1], st[1][x1]};
    simd_dbl x2v {st[0][x2], st[1][x2]};
    simd_dbl y1v {st[0][y1], st[1][y1]};
    simd_dbl y2v {st[0][y2], st[1][y2]};

    simd_dbl out {in[0], in[1]};
    out *= co[a0];

    // TODO: worth all that load and storing for this (?)
    out += a1v * x1v + a2v * x2v + b1v * y1v + b2v * y2v;

    for (uint i = 0; i < 2; ++i) {
      st[i][x2] = st[i][x1];
      st[i][x1] = in[i];
      st[i][y2] = st[i][y1];
      st[i][y1] = out[i];
    }
    return out;
  }
};

}} // namespace artv::liteon
#if 0
#include <algorithm>
#include <cmath>
#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace liteon {

struct presence_eq {
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::eq;
  //----------------------------------------------------------------------------
#if 0
  /* "mempool" emulates the JSFX per instance memory pool.
  It can be removed if unused. On most cases the length of this array can be
  reduced. It is set to a very conservative (high) value matching what JSFX
  provides. */
  double mempool[8 * 1024 * 1024] = {};
#endif
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::fabs (lhs - rhs) < 0.00001);
  }
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }

  //----------------------------------------------------------------------------
  // stubs for sliders
  // slider1:0<0,1,1{Stereo,Mono}>Processing
  constexpr double get_slider_slider1()
  {
    // TODO: stub, add code for getting "slider1"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    return 0.;
  }
  //----------------------------------------------------------------------------
  // slider2:7700<3100,18500,0.1>Frequency (Hz)
  double get_slider_slider2()
  {
    // TODO: stub, add code for getting "slider2"
    // Range: min:3100.0, max:18500.0, default: 7700.0, step: 0.1
    return frequency_p;
  }
  float frequency_p = 7700.;

  struct frequency_tag {};
  void set (frequency_tag, float v)
  {
    if (v == frequency_p) {
      return;
    }
    frequency_p = v;
    slider();
  }

  static constexpr auto get_parameter (frequency_tag)
  {
    return float_param ("Hz", 3100.f, 18500.f, 7700.f, 0.01f, 0.5f);
  }
  //----------------------------------------------------------------------------
  // slider3:0<-15,15,0.01>Cut/Boost (dB)
  double get_slider_slider3()
  {
    // TODO: stub, add code for getting "slider3"
    // Range: min:-15.0, max:15.0, default: 0.0, step: 0.01
    return gain_p;
  }
  float gain_p = 0;

  struct gain_tag {};
  void set (gain_tag, float v)
  {
    if (v == gain_p) {
      return;
    }
    gain_p = v;
    slider();
  }

  static constexpr auto get_parameter (gain_tag)
  {
    return float_param ("dB", -15.f, 15.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  // slider4:0.20<0.07,0.40,0.0001>BW (Min/Max)
  double get_slider_slider4()
  {
    // TODO: stub, add code for getting "slider4"
    // Range: min:0.07, max:0.4, default: 0.2, step: 0.0001
    return bandwidth_p;
  }
  float bandwidth_p = 0.07;

  struct bandwidth_tag {};
  void set (bandwidth_tag, float v)
  {
    if (v == bandwidth_p) {
      return;
    }
    bandwidth_p = v;
    slider();
  }

  static constexpr auto get_parameter (bandwidth_tag)
  {
    return float_param ("", 0.07f, 0.4f, 0.2f, 0.0001f);
  }
  //----------------------------------------------------------------------------
  // slider5:0<-25,25,0.05>Output (dB)
  constexpr double get_slider_slider5()
  {
    // TODO: stub, add code for getting "slider5"
    // Range: min:-25.0, max:25.0, default: 0.0, step: 0.05
    return 0.;
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<frequency_tag, bandwidth_tag, gain_tag>;
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double spn;
  double x1l;
  double x1r;
  double x2l;
  double x2r;
  double y1l;
  double y1r;
  double y2l;
  double y2r;
  double yl;
  double yr;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    spn = 0;
    x1l = 0;
    x1r = 0;
    x2l = 0;
    x2r = 0;
    y1l = 0;
    y1r = 0;
    y2l = 0;
    y2r = 0;
    yl  = 0;
    yr  = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double a0;
  double a1;
  double a2;
  double b1;
  double b2;
  double mono;
  double outgain;
  float  samplerate;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    a0      = 0;
    a1      = 0;
    a2      = 0;
    b1      = 0;
    b2      = 0;
    mono    = 0;
    outgain = 0;
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    init_init_variables();
    init_slider_variables();
    spn        = 0.;
    y2r        = spn;
    y1r        = y2r;
    x2r        = y1r;
    x1r        = x2r;
    yr         = x1r;
    y2l        = yr;
    y1l        = y2l;
    x2l        = y1l;
    x1l        = x2l;
    yl         = x1l;
    samplerate = pc.get_sample_rate();
    slider();
  }
  //----------------------------------------------------------------------------
  void slider()
  {
    double a        = 0.;
    double a2plus1  = 0.;
    double alphad   = 0.;
    double alphan   = 0.;
    double as2      = 0.;
    double as4      = 0.;
    double asnd     = 0.;
    double b0       = 0.;
    double boost    = 0.;
    double bw       = 0.;
    double c        = 0.;
    double ca       = 0.;
    double cf       = 0.;
    double cs       = 0.;
    double d        = 0.;
    double delta    = 0.;
    double f        = 0.;
    double f2       = 0.;
    double ma2plus1 = 0.;
    double mag      = 0.;
    double recipb0  = 0.;
    double sn       = 0.;
    double t        = 0.;
    double theta    = 0.;
    double tmp      = 0.;
    double xfmbw    = 0.;
    mono            = get_slider_slider1();
    cf              = get_slider_slider2() / samplerate;
    boost           = get_slider_slider3();
    bw              = get_slider_slider4();
    outgain         = eel2_pow (10., (get_slider_slider5() / 20.));
    ca              = std::tan (3.141592653589793 * (cf - 0.25));
    a               = eel2_pow (10., (boost / 20.));
    if ((boost < 6.0) && (boost > -6.0)) {
      (f = std::sqrt (a));
    }
    else {
      if (a > 1.) {
        (f = a / std::sqrt (2.));
      }
      else {
        (f = a * std::sqrt (2.));
      }
    }
    t     = std::tan (2. * 3.141592653589793 * bw);
    as2   = ca * ca;
    as4   = as2 * as2;
    d     = 2. * as2 * t;
    sn    = (1. + as4) * t;
    cs    = (1. - as4);
    mag   = std::sqrt (sn * sn + cs * cs);
    d     = mag;
    delta = atan2 (sn, cs) /*TODO: unknown call */;
    asnd  = (d <= 1. && d >= -1.) ? asin (d) : 0;
    theta = 0.5 * (3.141592653589793 - asnd - delta);
    tmp   = 0.5 * (asnd - delta);
    if (tmp > 0. && tmp < theta) {
      (theta = tmp);
    }
    xfmbw = theta / (2. * 3.141592653589793);
    c     = 1. / std::tan (2. * 3.141592653589793 * xfmbw);
    f2    = f * f;
    tmp   = a * a - f2;
#if 0
    if (std::fabs (tmp) <= spn) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#else
    if (std::fabs (tmp) <= 0.) {
      (alphad = c);
    }
    else {
      (alphad = std::sqrt (c * c * (f2 - 1.) / tmp));
    }
#endif

    alphan   = a * alphad;
    a2plus1  = 1. + as2;
    ma2plus1 = 1. - as2;
    a0       = a2plus1 + alphan * ma2plus1;
    a1       = 4.0 * ca;
    a2       = a2plus1 - alphan * ma2plus1;
    b0       = a2plus1 + alphad * ma2plus1;
    b2       = a2plus1 - alphad * ma2plus1;
    recipb0  = 1. / b0;
    a0 *= recipb0;
    a1 *= recipb0;
    a2 *= recipb0;
    b1 = a1;
    b2 *= recipb0;
    b1 = -b1;
    b2 = -b2;
    ;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    double xl = 0.;
    double xr = 0.;

    for (uint i = 0; i < samples; ++i) {
      auto& spl0 = chnls[0][i];
      auto& spl1 = chnls[1][i];
      if (eel2_eq (mono, 1.)) {
        xl   = (spl0 + spl1) / 2.;
        yl   = a0 * xl + a1 * x1l + a2 * x2l + b1 * y1l + b2 * y2l;
        x2l  = x1l;
        x1l  = xl;
        y2l  = y1l;
        y1l  = yl;
        spl1 = yl * outgain;
        spl0 = spl1;
      }
      else {
        xl   = spl0;
        xr   = spl1;
        yl   = a0 * xl + a1 * x1l + a2 * x2l + b1 * y1l + b2 * y2l;
        x2l  = x1l;
        x1l  = xl;
        y2l  = y1l;
        y1l  = yl;
        yr   = a0 * xr + a1 * x1r + a2 * x2r + b1 * y1r + b2 * y2r;
        x2r  = x1r;
        x1r  = xr;
        y2r  = y1r;
        y1r  = yr;
        spl0 = yl * outgain;
        spl1 = yr * outgain;
      }
    }
  }
}; /* presence_eq */
}} // namespace artv::liteon
#endif
