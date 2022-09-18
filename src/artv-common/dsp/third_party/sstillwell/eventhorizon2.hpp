#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/third_party/jsfx_engine/jsfx_engine.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace sstillwell {

class event_horizon_2 {
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
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::fabs (lhs - rhs) < 0.00001);
  }
  double jsfx_sign (double value)
  {
    auto v = *((uint64_t*) ((void*) &value));
    return (v == 0) ? 0. : (v & (1ull << 63)) ? -1. : 1.;
  }

  //----------------------------------------------------------------------------
  double jsfx_specialvar_get_srate() { return plugcontext->get_sample_rate(); }
  //----------------------------------------------------------------------------
public:
#if 0
  void set_slider1_slider (float v)
  {
    // Original slider line: slider1:0.0<-30.0,0.0,0.1>Threshold
    // Range: min:-30.0, max:0.0, default: 0.0, step: 0.1
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
    // Original slider line: slider1:0.0<-30.0,0.0,0.1>Threshold
    // Range: min:-30.0, max:0.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }

  static constexpr auto get_parameter (threshold_tag)
  {
    // Original slider line: slider1:0.0<-30.0,0.0,0.1>Threshold
    return float_param ("dB", -30.0, 0.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:-0.1<-20.0,0.0,0.1>Ceiling
    // Range: min:-20.0, max:0.0, default: -0.1, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct ceiling_tag {};

  void set (ceiling_tag, float v)
  {
    // Original slider line: slider2:-0.1<-20.0,0.0,0.1>Ceiling
    // Range: min:-20.0, max:0.0, default: -0.1, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (ceiling_tag)
  {
    // Original slider line: slider2:-0.1<-20.0,0.0,0.1>Ceiling
    return float_param ("dB", -20.0, 0.0, -0.1, 0.1);
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<0,1200,1>Release (ms)
    // Range: min:0.0, max:1200.0, default: 0.0, step: 1.0
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct release_tag {};

  void set (release_tag, float v)
  {
    // Original slider line: slider3:0<0,1200,1>Release (ms)
    // Range: min:0.0, max:1200.0, default: 0.0, step: 1.0
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (release_tag)
  {
    // Original slider line: slider3:0<0,1200,1>Release (ms)
    return float_param ("ms", 0.0, 1200.0, 0.0, 1.0);
  }

#endif
#if 0
#else
  // Snippet for parameter boilerplate in the authors framework....
  using parameters = mp11::mp_list<threshold_tag, ceiling_tag, release_tag>;
#endif
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
private:
  double atcoef;
  double attime;
  double db2log;
  double log2db;
  double pi;
  double relcoef;
  double reltime;
  double slider1;
  double slider2;
  double slider3;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    atcoef  = 0;
    attime  = 0;
    db2log  = 0;
    log2db  = 0;
    pi      = 0;
    relcoef = 0;
    reltime = 0;
    slider1 = 0;
    slider2 = 0;
    slider3 = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
private:
  double ceiling;
  double makeup;
  double overdb;
  double release;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    ceiling = 0;
    makeup  = 0;
    overdb  = 0;
    release = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
private:
  double ceilingdb;
  double maxover;
  double rundb;
  double runmax;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    ceilingdb = 0;
    maxover   = 0;
    rundb     = 0;
    runmax    = 0;
  }
  //----------------------------------------------------------------------------
  plugin_context* plugcontext;

public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_block_variables();
    init_slider_variables();

    slider1 = 0.0;
    slider2 = -0.1;
    slider3 = 0.0;
    pi      = 3.1415926535;
    log2db  = 8.6858896380650365530225783783321;
    db2log  = 0.11512925464970228420089957273422;
    attime  = 0.004;
    reltime = 0.200;
    atcoef  = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double ceildb   = 0.;
    double makeupdb = 0.;
    double peakdb   = 0.;
    double peaklvl  = 0.;
    double thresh   = 0.;
    double threshdb = 0.;
    thresh          = std::exp (slider1 * db2log);
    threshdb        = slider1;
    ceiling         = std::exp (slider2 * db2log);
    ceildb          = slider2;
    makeup          = std::exp ((ceildb - threshdb) * db2log);
    makeupdb        = ceildb - threshdb;
    peakdb          = ceildb + 25.;
    peaklvl         = std::exp (peakdb * db2log);
    release         = slider3 / 1000.;
    if (eel2_eq (release, 0.)) {
      (reltime = overdb / 125.);
    }
    else {
      (reltime = release);
    }
    relcoef = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double abs0    = 0.;
    double abs1    = 0.;
    double gr      = 0.;
    double overdbv = 0.;

    for (uint $$i = 0; $$i < samples; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];
      spl0 *= makeup;
      spl1 *= makeup;
      abs0    = std::abs (spl0);
      abs1    = std::abs (spl1);
      overdbv = std::max (abs0, abs1);
      overdb  = 2.08136898 * std::log (overdbv) * log2db - ceilingdb;
      if (overdb > rundb) {
        (rundb = overdb + atcoef * (rundb - overdb));
      }
      else {
        (rundb = overdb + relcoef * (rundb - overdb));
      }
      overdb = std::max (0., rundb);
      if (eel2_eq (release, 0.)) {
        if (overdb > maxover) {
          maxover = overdb;
          reltime = overdb / 125.;
          relcoef = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
        }
        runmax  = maxover + relcoef * (runmax - maxover);
        maxover = runmax;
      }
      gr = std::exp (overdb * db2log);
      spl0 *= gr;
      spl1 *= gr;
      spl0 = std::min<T> (ceiling, std::abs (spl0)) * jsfx_sign (spl0);
      spl1 = std::min<T> (ceiling, std::abs (spl1)) * jsfx_sign (spl1);
      ;
    }
  }
}; /* jsfx_process */

}} // namespace artv::sstillwell
