// Ported from https://github.com/JoepVanlier/JSFX.git
// commit sha: 5d54a165b806b2773792623376f6a51cc957368b
#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
// slider1:current_gain=0<-6,24,1>Gain (dB)
// slider2:current_ceiling=0<-18,0,1>Ceiling (dB)
// slider3:slew=0<0,1,.001>Slew
// slider4:warmth=0<-12,12,.001>Warmth (dB)

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace saike {

struct smooth {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::distortion;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
private:
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
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  static bool   eel2_ne (double lhs, double rhs) { return !eel2_eq (lhs, rhs); }
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }

  //----------------------------------------------------------------------------
  double jsfx_specialvar_get_srate() { return sample_rate; }
  //----------------------------------------------------------------------------
  // stubs for sliders
public:
  struct ceiling_tag {};
  void set (ceiling_tag, float v)
  {
    ceilingp = v;
    slider();
  }
  static constexpr auto get_parameter (ceiling_tag)
  {
    return float_param ("dB", -18.f, 0.f, 0.f, 0.1f);
  }
  float ceilingp;

  double get_slider_current_ceiling()
  {
    // TODO: stub, add code for getting "current_ceiling"
    // Range: min:-18.0, max:0.0, default: 0.0, step: 1.0
    return ceilingp;
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    drivep = v;
    slider();
  }
  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -6.f, 24.f, 0.f, 0.1f);
  }
  float  drivep;
  double get_slider_current_gain()
  {
    // TODO: stub, add code for getting "current_gain"
    // Range: min:-6.0, max:24.0, default: 0.0, step: 1.0
    return drivep;
  }
  //----------------------------------------------------------------------------
  struct slew_tag {};
  void set (slew_tag, float v)
  {
    slewp = v;
    slider();
  }
  static constexpr auto get_parameter (slew_tag)
  {
    return float_param ("", -0.f, 1.f, 0.f, 0.001f);
  }
  float  slewp;
  double get_slider_slew()
  {
    // TODO: stub, add code for getting "slew"
    // Range: min:0.0, max:1.0, default: 0.0, step: 0.001
    return slewp;
  }
  //----------------------------------------------------------------------------
  struct warmth_tag {};
  void set (warmth_tag, float v)
  {
    warmthp = v;
    slider();
  }
  static constexpr auto get_parameter (warmth_tag)
  {
    return float_param ("dB", -12.f, 12.f, 0.f, 0.01f);
  }
  float  warmthp;
  double get_slider_warmth()
  {
    // TODO: stub, add code for getting "warmth"
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.001
    return warmthp;
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<drive_tag, ceiling_tag, slew_tag, warmth_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double dc0$a1_1 {};
  double dc0$a1_2 {};
  double dc0$a1_3 {};
  double dc0$a2_1 {};
  double dc0$a2_2 {};
  double dc0$a2_3 {};
  double dc0$g {};
  double dc0$k1 {};
  double dc0$k2 {};
  double dc0$k3 {};
  double dc1$a1_1 {};
  double dc1$a1_2 {};
  double dc1$a1_3 {};
  double dc1$a2_1 {};
  double dc1$a2_2 {};
  double dc1$a2_3 {};
  double dc1$g {};
  double dc1$k1 {};
  double dc1$k2 {};
  double dc1$k3 {};
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double ceilingv {};
  double inv_ceiling {};
  double inverse_tilt_gain {};
  double l_itilt$a1 {};
  double l_itilt$a2 {};
  double l_itilt$a3 {};
  double l_itilt$m0 {};
  double l_itilt$m1 {};
  double l_itilt$m2 {};
  double l_tilt$a1 {};
  double l_tilt$a2 {};
  double l_tilt$a3 {};
  double l_tilt$m0 {};
  double l_tilt$m1 {};
  double l_tilt$m2 {};
  double preamp {};
  double r_itilt$a1 {};
  double r_itilt$a2 {};
  double r_itilt$a3 {};
  double r_itilt$m0 {};
  double r_itilt$m1 {};
  double r_itilt$m2 {};
  double r_tilt$a1 {};
  double r_tilt$a2 {};
  double r_tilt$a3 {};
  double r_tilt$m0 {};
  double r_tilt$m1 {};
  double r_tilt$m2 {};
  double slew_target {};
  double tilt_gain {};
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double dc0$ic1eq {};
  double dc0$ic2eq {};
  double dc0$ic3eq {};
  double dc0$ic4eq {};
  double dc0$ic5eq {};
  double dc0$ic6eq {};
  double l_itilt$ic1eq {};
  double l_itilt$ic2eq {};
  double l_tilt$ic1eq {};
  double l_tilt$ic2eq {};
  double r_itilt$ic1eq {};
  double r_itilt$ic2eq {};
  double r_tilt$ic1eq {};
  double r_tilt$ic2eq {};
  double slew_left$lx {};
  double slew_left$slew_factor {};
  double slew_right$lx {};
  double slew_right$slew_factor {};

  float sample_rate = 44100;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    jsfx_process_reset();
  }
  //----------------------------------------------------------------------------
private:
  void jsfx_process_reset()
  {
    init$init_HP6 (
      7. / jsfx_specialvar_get_srate(),
      dc0$g,
      dc0$k1,
      dc0$a1_1,
      dc0$k2,
      dc0$a1_2,
      dc0$k3,
      dc0$a1_3,
      dc0$a2_1,
      dc0$a2_2,
      dc0$a2_3);
    init$init_HP6 (
      7. / jsfx_specialvar_get_srate(),
      dc1$g,
      dc1$k1,
      dc1$a1_1,
      dc1$k2,
      dc1$a1_2,
      dc1$k3,
      dc1$a1_3,
      dc1$a2_1,
      dc1$a2_2,
      dc1$a2_3);
    ;
  }
  //----------------------------------------------------------------------------
  void slider()
  {
    double log10d20_conversion = 0.;
    double omega_tilt          = 0.;
    log10d20_conversion        = .11512925464970228420089957273422;
    preamp   = std::exp (log10d20_conversion * get_slider_current_gain());
    ceilingv = std::exp (-log10d20_conversion * get_slider_current_ceiling());
    slew_target
      = 1.0 - std::exp (-log10d20_conversion * 15. * get_slider_slew());
    inv_ceiling       = 1.0 / ceilingv;
    tilt_gain         = eel2_pow (10., (get_slider_warmth() / 20.));
    inverse_tilt_gain = 1.0 / tilt_gain;
    omega_tilt        = 3200. / jsfx_specialvar_get_srate();
    init$init_tilt (
      omega_tilt,
      0.,
      tilt_gain,
      l_tilt$a1,
      l_tilt$a2,
      l_tilt$a3,
      l_tilt$m0,
      l_tilt$m1,
      l_tilt$m2);
    init$init_tilt (
      omega_tilt,
      0.,
      tilt_gain,
      r_tilt$a1,
      r_tilt$a2,
      r_tilt$a3,
      r_tilt$m0,
      r_tilt$m1,
      r_tilt$m2);
    init$init_tilt (
      omega_tilt,
      0.,
      inverse_tilt_gain,
      l_itilt$a1,
      l_itilt$a2,
      l_itilt$a3,
      l_itilt$m0,
      l_itilt$m1,
      l_itilt$m2);
    init$init_tilt (
      omega_tilt,
      0.,
      inverse_tilt_gain,
      r_itilt$a1,
      r_itilt$a2,
      r_itilt$a3,
      r_itilt$m0,
      r_itilt$m1,
      r_itilt$m2);
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    T l = 0.;
    T r = 0.;

    for (int i = 0; i < samples; ++i) {
      auto& spl0 = outs[0][i];
      auto& spl1 = outs[1][i];
      spl0       = ins[0][i];
      spl1       = ins[1][i];
      l          = spl0;
      r          = spl1;
      if (get_slider_slew() > 0.) {
        l = block$slew_buffer (l, slew_left$slew_factor, slew_left$lx);
        r = block$slew_buffer (r, slew_right$slew_factor, slew_right$lx);
      }
      l = ceilingv * preamp * l;
      r = ceilingv * preamp * r;
      if (eel2_ne (get_slider_warmth(), 0.)) {
        l = init$eval_tilt (
              l,
              l_itilt$ic2eq,
              l_itilt$a1,
              l_itilt$ic1eq,
              l_itilt$a2,
              l_itilt$a3,
              l_itilt$m1,
              l_itilt$m2)
          * tilt_gain;
        r = init$eval_tilt (
              r,
              r_itilt$ic2eq,
              r_itilt$a1,
              r_itilt$ic1eq,
              r_itilt$a2,
              r_itilt$a3,
              r_itilt$m1,
              r_itilt$m2)
          * tilt_gain;
      }
      l = block$processChannel (l) * inv_ceiling;
      r = block$processChannel (r) * inv_ceiling;
      if (eel2_ne (get_slider_warmth(), 0.)) {
        l = init$eval_tilt (
              l,
              l_tilt$ic2eq,
              l_tilt$a1,
              l_tilt$ic1eq,
              l_tilt$a2,
              l_tilt$a3,
              l_tilt$m1,
              l_tilt$m2)
          * inverse_tilt_gain;
        r = init$eval_tilt (
              r,
              r_tilt$ic2eq,
              r_tilt$a1,
              r_tilt$ic1eq,
              r_tilt$a2,
              r_tilt$a3,
              r_tilt$m1,
              r_tilt$m2)
          * inverse_tilt_gain;
      }
      spl0 = init$eval_HP6 (
        l,
        dc0$a1_1,
        dc0$ic1eq,
        dc0$a2_1,
        dc0$ic2eq,
        dc0$g,
        dc0$k1,
        dc0$a1_2,
        dc0$ic3eq,
        dc0$a2_2,
        dc0$ic4eq,
        dc0$k2,
        dc0$a1_3,
        dc0$ic5eq,
        dc0$a2_3,
        dc0$ic6eq,
        dc0$k3);
      spl1 = init$eval_HP6 (
        r,
        dc0$a1_1,
        dc0$ic1eq,
        dc0$a2_1,
        dc0$ic2eq,
        dc0$g,
        dc0$k1,
        dc0$a1_2,
        dc0$ic3eq,
        dc0$a2_2,
        dc0$ic4eq,
        dc0$k2,
        dc0$a1_3,
        dc0$ic5eq,
        dc0$a2_3,
        dc0$ic6eq,
        dc0$k3);
      ;
      // addon, keep gain constant
      spl0 /= preamp;
      spl1 /= preamp;
    }
  }
  // functions for section "init"
private:
  //----------------------------------------------------------------------------
  double init$eval_HP6 (
    double  v0,
    double& a1_1,
    double& ic1eq,
    double& a2_1,
    double& ic2eq,
    double& g,
    double& k1,
    double& a1_2,
    double& ic3eq,
    double& a2_2,
    double& ic4eq,
    double& k2,
    double& a1_3,
    double& ic5eq,
    double& a2_3,
    double& ic6eq,
    double& k3)
  {
    double v1 = 0.;
    double v2 = 0.;
    double hp = 0.;
    v1        = a1_1 * ic1eq + a2_1 * (v0 - ic2eq);
    v2        = ic2eq + g * v1;
    ic1eq     = 2. * v1 - ic1eq;
    ic2eq     = 2. * v2 - ic2eq;
    hp        = v0 - k1 * v1 - v2;
    v1        = a1_2 * ic3eq + a2_2 * (hp - ic4eq);
    v2        = ic4eq + g * v1;
    ic3eq     = 2. * v1 - ic3eq;
    ic4eq     = 2. * v2 - ic4eq;
    hp        = hp - k2 * v1 - v2;
    v1        = a1_3 * ic5eq + a2_3 * (hp - ic6eq);
    v2        = ic6eq + g * v1;
    ic5eq     = 2. * v1 - ic5eq;
    ic6eq     = 2. * v2 - ic6eq;
    hp        = hp - k3 * v1 - v2;
    return hp;
  }
  //----------------------------------------------------------------------------
  double init$eval_tilt (
    double  v0,
    double& ic2eq,
    double& a1,
    double& ic1eq,
    double& a2,
    double& a3,
    double& m1,
    double& m2)
  {
    double v3 = 0.;
    double v1 = 0.;
    double v2 = 0.;
    v3        = v0 - ic2eq;
    v1        = a1 * ic1eq + a2 * v3;
    v2        = ic2eq + a2 * ic1eq + a3 * v3;
    ic1eq     = 2. * v1 - ic1eq;
    ic2eq     = 2. * v2 - ic2eq;
    return (v0 + m1 * v1 + m2 * v2);
  }
  //----------------------------------------------------------------------------
  double init$init_HP6 (
    double  freq,
    double& g,
    double& k1,
    double& a1_1,
    double& k2,
    double& a1_2,
    double& k3,
    double& a1_3,
    double& a2_1,
    double& a2_2,
    double& a2_3)
  {
    double res = 0.;
    g          = std::tan (3.141592653589793 * freq);
    k1         = 1.93185165257814;
    a1_1       = 1. / (1. + g * (g + k1));
    k2         = 1.41421356474619;
    a1_2       = 1. / (1. + g * (g + k2));
    k3         = 0.517638090205042;
    a1_3       = 1. / (1. + g * (g + k3));
    a2_1       = g * a1_1;
    a2_2       = g * a1_2;
    a2_3       = g * a1_3;
    return a2_3;
  }
  //----------------------------------------------------------------------------
  double init$init_tilt (
    double  freq,
    double  res,
    double  A,
    double& a1,
    double& a2,
    double& a3,
    double& m0,
    double& m1,
    double& m2)
  {
    double k = 0.;
    double g = 0.;
    g        = std::tan (.5 * 3.141592653589793 * freq) / std::sqrt (A);
    k        = 2. - 2. * res;
    a1       = 1. / (1. + g * (g + k));
    a2       = g * a1;
    a3       = g * a2;
    m0       = 1.;
    m1       = k * (A - 1.);
    m2       = (A * A - 1.);
    return m2;
  }
  // functions for section "block"

  //----------------------------------------------------------------------------
  double block$processChannel (double x)
  {
    double p1    = 0.;
    double p2    = 0.;
    double p6    = 0.;
    double p4    = 0.;
    double p5    = 0.;
    double p0    = 0.;
    double indep = 0.;
    double cterm = 0.;
    double p7    = 0.;
    double p3    = 0.;
    p0           = -0.43497619473852767;
    p1           = -0.5460978062661744;
    p2           = -0.9288977638330088;
    p3           = -1.8931725230583496;
    p4           = -1.901829522587536;
    p5           = -0.7951494375999134;
    p6           = -0.35701520925330404;
    indep        = p2 + p3 * x + p6 * x * x;
    cterm        = [&] {
      if (indep < 0.) {
        return 1.;
      }
      else {
        return [&] {
          if (cterm > 3.141592653589793) {
            return -1.;
          }
          else {
            return std::cos (indep);
          }
        }();
      }
    }();
    return p4 * block$tanh (p0 + p1 * x) + p5 * cterm;
  }
  //----------------------------------------------------------------------------
  double block$slew_buffer (double x, double& slew_factor, double& lx)
  {
    double diff = 0.;
    slew_factor = .999 * slew_factor
      + .001 * (0.00075 * jsfx_specialvar_get_srate() * slew_target);
    diff = block$tanh (slew_factor * (x - lx)) / slew_factor;
    lx += diff;
    return lx;
  }
  //----------------------------------------------------------------------------
  double block$tanh (double x) { return 2. / (1. + std::exp (-2. * x)) - 1.; }
}; /* jsfx_process */

}} // namespace artv::saike
