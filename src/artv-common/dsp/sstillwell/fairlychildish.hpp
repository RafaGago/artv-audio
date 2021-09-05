#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace sstillwell {

struct fairly_childish {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::dynamics;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double eel2_and (double lhs, double rhs)
  {
    return (double) ((uint64_t) lhs & (uint64_t) rhs);
  }
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }

  //----------------------------------------------------------------------------
  // stubs for JSFX special variables

  double jsfx_specialvar_get_srate() { return plugcontext->get_sample_rate(); }

  //----------------------------------------------------------------------------
public:
#if 0
  void set_slider1_slider (float v)
  {
    // Original slider line: slider1:0<-60,0,0.1>Threshold (dB)
    // Range: min:-60.0, max:0.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct threshold_tag {};

  void set (threshold_tag, float v)
  {
    // Original slider line: slider1:0<-60,0,0.1>Threshold (dB)
    // Range: min:-60.0, max:0.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }

  static constexpr auto get_parameter (threshold_tag)
  {
    // Original slider line: slider1:0<-60,0,0.1>Threshold (dB)
    return float_param ("dB", -60.0, 0.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:70<0.1,100,0.1>Bias
    // Range: min:0.1, max:100.0, default: 70.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct bias_tag {};

  void set (bias_tag, float v)
  {
    // Original slider line: slider2:70<0.1,100,0.1>Bias
    // Range: min:0.1, max:100.0, default: 70.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (bias_tag)
  {
    // Original slider line: slider2:70<0.1,100,0.1>Bias
    return float_param ("", 0.1, 100.0, 70.0, 0.1);
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<-30,30,0.1>Makeup Gain
    // Range: min:-30.0, max:30.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct makeup_tag {};

  void set (makeup_tag, float v)
  {
    // Original slider line: slider3:0<-30,30,0.1>Makeup Gain
    // Range: min:-30.0, max:30.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (makeup_tag)
  {
    // Original slider line: slider3:0<-30,30,0.1>Makeup Gain
    return float_param ("dB", -30.0, 30.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider4_slider (float v)
  {
    // Original slider line: slider4:2<0,3,1{Left/Right (Blown
    // Capacitor),Lat/Vert (Blown Capacitor),Left/Right,Lat/Vert}>AGC Range:
    // min:0.0, max:3.0, default: 2.0, step: 1.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct agc_range_tag {};

  void set (agc_range_tag, int v)
  {
    v += 2;
    // Original slider line: slider4:2<0,3,1{Left/Right (Blown
    // Capacitor),Lat/Vert (Blown Capacitor),Left/Right,Lat/Vert}>AGC Range:
    // min:0.0, max:3.0, default: 2.0, step: 1.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }

  static constexpr auto get_parameter (agc_range_tag)
  {
    // Original slider line: slider4:2<0,3,1{Left/Right (Blown
    // Capacitor),Lat/Vert (Blown Capacitor),Left/Right,Lat/Vert}>AGC
    return choice_param (0, make_cstr_array ("Left/Right", "Lat/Vert"));
  }

#endif
#if 0
  void set_slider5_slider (float v)
  {
    // Original slider line: slider5:1<1,6,1>Time Constant
    // Range: min:1.0, max:6.0, default: 1.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct time_constant_tag {};

  void set (time_constant_tag, float v)
  {
    // Original slider line: slider5:1<1,6,1>Time Constant
    // Range: min:1.0, max:6.0, default: 1.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }

  static constexpr auto get_parameter (time_constant_tag)
  {
    // Original slider line: slider5:1<1,6,1>Time Constant
    return float_param ("", 1.0, 6.0, 1.0, 1.0);
  }

#endif
#if 0
  void set_slider6_slider (float v)
  {
    // Original slider line: slider6:100<1,10000,1>Level Detector RMS Window
    // Range: min:1.0, max:10000.0, default: 100.0, step: 1.0
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct rms_window_tag {};

  void set (rms_window_tag, float v)
  {
    // Original slider line: slider6:100<1,10000,1>Level Detector RMS Window
    // Range: min:1.0, max:10000.0, default: 100.0, step: 1.0
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }

  static constexpr auto get_parameter (rms_window_tag)
  {
    // Original slider line: slider6:100<1,10000,1>Level Detector RMS Window
    return float_param ("", 1.0, 10000.0, 100.0, 1.0);
  }

#endif
#if 0 // these are meters
#if 0
  void set_slider7_slider (float v)
  {
    // Original slider line: slider7:1<1,50,0.1>Current Compression Ratio
    // Range: min:1.0, max:50.0, default: 1.0, step: 0.1
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct ratio_tag {};

  void set (ratio_tag, float v)
  {
    // Original slider line: slider7:1<1,50,0.1>Current Compression Ratio
    // Range: min:1.0, max:50.0, default: 1.0, step: 0.1
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }

  static constexpr auto get_parameter (ratio_tag)
  {
    // Original slider line: slider7:1<1,50,0.1>Current Compression Ratio
    return float_param ("1:", 1.0, 50.0, 1.0, 0.01, 0.5);
  }

#endif
#if 0
  void set_slider8_slider (float v)
  {
    // Original slider line: slider8:0<-90,0,0.1>Gain Reduction
    // Range: min:-90.0, max:0.0, default: 0.0, step: 0.1
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct reduction_tag {};

  void set (reduction_tag, float v)
  {
    // Original slider line: slider8:0<-90,0,0.1>Gain Reduction
    // Range: min:-90.0, max:0.0, default: 0.0, step: 0.1
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }

  static constexpr auto get_parameter (reduction_tag)
  {
    // Original slider line: slider8:0<-90,0,0.1>Gain Reduction
    return float_param ("", -90.0, 0.0, 0.0, 0.1);
  }

#endif
#endif // these are meters
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    threshold_tag,
    bias_tag,
    makeup_tag,
    agc_range_tag,
    time_constant_tag,
    rms_window_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double atcoef;
  double attime;
  double cratio;
  double db2log;
  double i;
  double latvert;
  double leftright;
  double log2db;
  double maxover;
  double overdb;
  double ratio;
  double relcoef;
  double reltime;
  double rmscoef;
  double rmstime;
  double rundb;
  double slider1;
  double slider2;
  double slider3;
  double slider4;
  double slider5;
  double slider6;
#if 0
  double slider7;
  double slider8;
#endif
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    atcoef    = 0;
    attime    = 0;
    cratio    = 0;
    db2log    = 0;
    i         = 0;
    latvert   = 0;
    leftright = 0;
    log2db    = 0;
    maxover   = 0;
    overdb    = 0;
    ratio     = 0;
    relcoef   = 0;
    reltime   = 0;
    rmscoef   = 0;
    rmstime   = 0;
    rundb     = 0;
    slider1   = 0;
    slider2   = 0;
    slider3   = 0;
    slider4   = 0;
    slider5   = 0;
    slider6   = 0;
#if 0
    slider7   = 0;
    slider8   = 0;
#endif
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double agc;
  double bias;
  double capsc;
  double makeupv;
  double threshv;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    agc     = 0;
    bias    = 0;
    capsc   = 0;
    makeupv = 0;
    threshv = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double dcoffset;
  double runave;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    dcoffset = 0;
    runave   = 0;
  }
  plugin_context* plugcontext = nullptr;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_slider_variables();
    init_block_variables();
    slider1 = 0.0;
    slider2 = 70.0;
    slider3 = 0.0;
    slider4 = 2.0;
    slider5 = 1.0;
    slider6 = 100.0;
#if 0
    slider7   = 1.0;
    slider8   = 0.0;
#endif
    log2db    = 8.6858896380650365530225783783321;
    db2log    = 0.11512925464970228420089957273422;
    i         = 0.;
    attime    = 0.0002;
    reltime   = 0.300;
    rmstime   = 0.000050;
    maxover   = 0.;
    ratio     = 0.;
    cratio    = 0.;
    rundb     = 0.;
    overdb    = 0.;
    maxover   = 0.;
    atcoef    = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef   = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    rmscoef   = std::exp (-1. / (rmstime * jsfx_specialvar_get_srate()));
    leftright = 0.;
    latvert   = 1.;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double cthresh      = 0.;
    double cthreshv     = 0.;
    double makeup       = 0.;
    double thresh       = 0.;
    double timeconstant = 0.;
    thresh              = slider1;
    threshv             = std::exp (thresh * db2log);
    ratio               = 20.;
    bias                = 80. * slider2 / 100.;
    cthresh             = thresh - bias;
    cthreshv            = std::exp (cthresh * db2log);
    makeup              = slider3;
    makeupv             = std::exp (makeup * db2log);
    agc                 = eel2_and (slider4, 1.);
    capsc               = [&] {
      if (eel2_and (slider4, 2.)) {
        return log2db;
      }
      else {
        return log2db * 2.08136898;
      }
    }();
    timeconstant = slider5;
    if (eel2_eq (timeconstant, 1.)) {
      attime  = 0.0002;
      reltime = 0.300;
    }
    if (eel2_eq (timeconstant, 2.)) {
      attime  = 0.0002;
      reltime = 0.800;
    }
    if (eel2_eq (timeconstant, 3.)) {
      attime  = 0.0004;
      reltime = 2.000;
    }
    if (eel2_eq (timeconstant, 4.)) {
      attime  = 0.0008;
      reltime = 5.000;
    }
    if (eel2_eq (timeconstant, 5.)) {
      attime  = 0.0002;
      reltime = 10.000;
    }
    if (eel2_eq (timeconstant, 6.)) {
      attime  = 0.0004;
      reltime = 25.000;
    }
    atcoef  = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    rmstime = slider6 / 1000000.;
    rmscoef = std::exp (-1. / (rmstime * jsfx_specialvar_get_srate()));
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double aspl0  = 0.;
    double aspl1  = 0.;
    double det    = 0.;
    double gr     = 0.;
    double grv    = 0.;
    double maxspl = 0.;
    double sav0   = 0.;
    double sav1   = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];
      if (eel2_eq (agc, leftright)) {
        aspl0 = std::abs (spl0);
        aspl1 = std::abs (spl1);
      }
      else {
        aspl0 = std::abs (spl0 + spl1) / 2.;
        aspl1 = std::abs (spl0 - spl1) / 2.;
      }
      maxspl = std::max (aspl0, aspl1);
      maxspl *= maxspl;
      runave = maxspl + rmscoef * (runave - maxspl);
      det    = std::sqrt (std::max (0., runave));
      overdb = capsc * std::log (det / threshv);
      overdb = std::max (0., overdb);
      if (overdb > rundb) {
        (rundb = overdb + atcoef * (rundb - overdb));
      }
      else {
        (rundb = overdb + relcoef * (rundb - overdb));
      }
      overdb = std::max (rundb, 0.);
      if (eel2_eq (bias, 0.)) {
        (cratio = ratio);
      }
      else {
        (cratio = 1.
           + (ratio - 1.)
             * std::sqrt ((overdb + dcoffset) / (bias + dcoffset)));
      }
#if 0
      slider7 = cratio;
#endif
      gr = -overdb * (cratio - 1.) / cratio;
#if 0
      slider8 = -gr;
#endif
      grv = std::exp (gr * db2log);
      if (eel2_eq (agc, leftright)) {
        spl0 *= grv * makeupv;
        spl1 *= grv * makeupv;
      }
      else {
        sav0 = (spl0 + spl1) * grv;
        sav1 = (spl0 - spl1) * grv;
        spl0 = makeupv * (sav0 + sav1) * 0.5;
        spl1 = makeupv * (sav0 - sav1) * 0.5;
      };
    }
  }
}; /* jsfx_process */
}} // namespace artv::sstillwell
