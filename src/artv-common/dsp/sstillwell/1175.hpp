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

struct _1175 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::dynamics;
  //----------------------------------------------------------------------------
private:
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
    return float_param ("dB", -60.0, 0.0, 0.0, 0.1, 1.6);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:5<0,9,1{Blown Capacitor 4
    // (Deprecated),Blown Capacitor 8 (Deprecated),Blown Capacitor 12
    // (Deprecated),Blown Capacitor 20 (Deprecated),Blown Capacitor All
    // (Deprecated),4,8,12,20,All}>Ratio Range: min:0.0, max:9.0, default: 5.0,
    // step: 1.0
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct ratio_tag {};

  void set (ratio_tag, int v)
  {
    v += 5;
    // Original slider line: slider2:5<0,9,1{Blown Capacitor 4
    // (Deprecated),Blown Capacitor 8 (Deprecated),Blown Capacitor 12
    // (Deprecated),Blown Capacitor 20 (Deprecated),Blown Capacitor All
    // (Deprecated),4,8,12,20,All}>Ratio Range: min:0.0, max:9.0, default: 5.0,
    // step: 1.0
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (ratio_tag)
  {
    // Original slider line: slider2:5<0,9,1{Blown Capacitor 4
    // (Deprecated),Blown Capacitor 8 (Deprecated),Blown Capacitor 12
    // (Deprecated),Blown Capacitor 20 (Deprecated),Blown Capacitor All
    // (Deprecated),4,8,12,20,All}>Ratio
    return choice_param (
      0, make_cstr_array ("1:4", "1:8", "1:12", "1:20", "All"));
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<-20,20,0.1>Gain (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gain_tag {};

  void set (gain_tag, float v)
  {
    // Original slider line: slider3:0<-20,20,0.1>Gain (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (gain_tag)
  {
    // Original slider line: slider3:0<-20,20,0.1>Gain (dB)
    return float_param ("dB", -20.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider4_slider (float v)
  {
    // Original slider line: slider4:20<20,2000,10>Attack (uS)
    // Range: min:20.0, max:2000.0, default: 20.0, step: 10.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct attack_tag {};

  void set (attack_tag, float v)
  {
    // Original slider line: slider4:20<20,2000,10>Attack (uS)
    // Range: min:20.0, max:2000.0, default: 20.0, step: 10.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }

  static constexpr auto get_parameter (attack_tag)
  {
    // Original slider line: slider4:20<20,2000,10>Attack (uS)
    return float_param ("uS", 20.0, 2000.0, 20.0, 1.0, 0.7);
  }

#endif
#if 0
  void set_slider5_slider (float v)
  {
    // Original slider line: slider5:250<20,1000,1>Release (mS)
    // Range: min:20.0, max:1000.0, default: 250.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct release_tag {};

  void set (release_tag, float v)
  {
    // Original slider line: slider5:250<20,1000,1>Release (mS)
    // Range: min:20.0, max:1000.0, default: 250.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }

  static constexpr auto get_parameter (release_tag)
  {
    // Original slider line: slider5:250<20,1000,1>Release (mS)
    return float_param ("mS", 20.0, 1000.0, 250.0, 1.0, 0.7);
  }

#endif
#if 0
  void set_slider6_slider (float v)
  {
    // Original slider line: slider6:100<0,100,0.1>Mix (%)
    // Range: min:0.0, max:100.0, default: 100.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mix_tag {};

  void set (mix_tag, float v)
  {
    // Original slider line: slider6:100<0,100,0.1>Mix (%)
    // Range: min:0.0, max:100.0, default: 100.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }

  static constexpr auto get_parameter (mix_tag)
  {
    // Original slider line: slider6:100<0,100,0.1>Mix (%)
    return float_param ("%", 0.0, 100.0, 100.0, 0.1);
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    threshold_tag,
    ratio_tag,
    gain_tag,
    attack_tag,
    release_tag,
    mix_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double atcoef;
  double attime;
  double cratio;
  double db2log;
  double gr_meter;
  double gr_meter_decay;
  double log2db;
  double mix;
  double overdb;
  double ratatcoef;
  double ratio;
  double ratrelcoef;
  double relcoef;
  double reltime;
  double rundb;
  double slider1;
  double slider2;
  double slider3;
  double slider4;
  double slider5;
  double slider6;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    atcoef         = 0;
    attime         = 0;
    cratio         = 0;
    db2log         = 0;
    gr_meter       = 0;
    gr_meter_decay = 0;
    log2db         = 0;
    mix            = 0;
    overdb         = 0;
    ratatcoef      = 0;
    ratio          = 0;
    ratrelcoef     = 0;
    relcoef        = 0;
    reltime        = 0;
    rundb          = 0;
    slider1        = 0;
    slider2        = 0;
    slider3        = 0;
    slider4        = 0;
    slider5        = 0;
    slider6        = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double allin;
  double autogain;
  double capsc;
  double cthreshv;
  double makeupv;
  double softknee;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    allin    = 0;
    autogain = 0;
    capsc    = 0;
    cthreshv = 0;
    makeupv  = 0;
    softknee = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double maxover;
  double rmscoef;
  double runave;
  double runmax;
  double runratio;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    maxover  = 0;
    rmscoef  = 0;
    runave   = 0;
    runmax   = 0;
    runratio = 0;
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
    slider1        = 0.0;
    slider2        = 5.0;
    slider3        = 0.0;
    slider4        = 20.0;
    slider5        = 250.0;
    slider6        = 100.0;
    log2db         = 8.6858896380650365530225783783321;
    db2log         = 0.11512925464970228420089957273422;
    attime         = 0.010;
    reltime        = 0.100;
    ratio          = 0.;
    cratio         = 0.;
    rundb          = 0.;
    overdb         = 0.;
    ratatcoef      = std::exp (-1. / (0.00001 * jsfx_specialvar_get_srate()));
    ratrelcoef     = std::exp (-1. / (0.5 * jsfx_specialvar_get_srate()));
    atcoef         = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef        = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    mix            = 1.;
    gr_meter       = 1.;
    gr_meter_decay = std::exp (1. / (1. * jsfx_specialvar_get_srate()));
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double cthresh = 0.;
    double makeup  = 0.;
    double rpos    = 0.;
    double thresh  = 0.;
    double threshv = 0.;
    thresh         = slider1;
    threshv        = std::exp (thresh * db2log);
    capsc          = log2db;
    if ((rpos = slider2) > 4.) {
      (rpos -= 5.);
    }
    else {
      (capsc *= 2.08136898);
    }
    ratio = [&] {
      if (eel2_eq (rpos, 0.)) {
        return 4.;
      }
      else {
        return [&] {
          if (eel2_eq (rpos, 1.)) {
            return 8.;
          }
          else {
            return [&] {
              if (eel2_eq (rpos, 2.)) {
                return 12.;
              }
              else {
                return [&] {
                  if (eel2_eq (rpos, 3.)) {
                    return 20.;
                  }
                  else {
                    return 20.;
                  }
                }();
              }
            }();
          }
        }();
      }
    }();
    if (eel2_eq (rpos, 4.)) {
      allin  = 1.;
      cratio = 20.;
    }
    else {
      allin  = 0.;
      cratio = ratio;
    }
    cthresh = [&] {
      if (softknee) {
        return (thresh - 3.);
      }
      else {
        return thresh;
      }
    }();
    cthreshv = std::exp (cthresh * db2log);
    makeup   = slider3;
    makeupv  = std::exp ((makeup + autogain) * db2log);
    attime   = slider4 / 1000000.;
    reltime  = slider5 / 1000.;
    atcoef   = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef  = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    mix      = slider6 / 100.;
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint block_samples)
  {

    double aspl0    = 0.;
    double aspl1    = 0.;
    double averatio = 0.;
    double det      = 0.;
    double gr       = 0.;
    double grv      = 0.;
    double maxspl   = 0.;
    double ospl0    = 0.;
    double ospl1    = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = chnls[0][$$i];
      auto& spl1 = chnls[1][$$i];
      ospl0      = spl0;
      ospl1      = spl1;
      aspl0      = std::abs (spl0);
      aspl1      = std::abs (spl1);
      maxspl     = std::max (aspl0, aspl1);
      maxspl *= maxspl;
      runave = maxspl + rmscoef * (runave - maxspl);
      det    = std::sqrt (std::max (0., runave));
      overdb = std::max (0., capsc * std::log (det / cthreshv));
      if (overdb - rundb > 5.) {
        (averatio = 4.);
      }
      if (overdb > rundb) {
        rundb    = overdb + atcoef * (rundb - overdb);
        runratio = averatio + ratatcoef * (runratio - averatio);
      }
      else {
        rundb    = overdb + relcoef * (rundb - overdb);
        runratio = averatio + ratrelcoef * (runratio - averatio);
      }
      overdb   = rundb;
      averatio = runratio;
      if (allin) {
        (cratio = 12. + averatio);
      }
      else {
        (cratio = ratio);
      }
      gr      = -overdb * (cratio - 1.) / cratio;
      grv     = std::exp (gr * db2log);
      runmax  = maxover + relcoef * (runmax - maxover);
      maxover = runmax;
      if (grv < gr_meter) {
        gr_meter = grv;
      }
      else {
        gr_meter *= gr_meter_decay;
        if (gr_meter > 1.) {
          gr_meter = 1.;
        }
      }
      spl0 *= grv * makeupv * mix;
      spl1 *= grv * makeupv * mix;
      spl0 += ospl0 * (1. - mix);
      spl1 += ospl1 * (1. - mix);
      ;
    }
  }
}; /* jsfx_process */
}} // namespace artv::sstillwell
