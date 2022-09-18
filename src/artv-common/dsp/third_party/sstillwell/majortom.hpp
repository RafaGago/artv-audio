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
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace sstillwell {

struct major_tom {
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
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }
  std::vector<float> heapmem;
  inline float&      heap (std::size_t value) { return heapmem[value]; }
  void               heap_reset (std::size_t s)
  {
    heapmem.resize (s);
    std::memset (heapmem.data(), 0, heapmem.size() * sizeof heapmem[0]);
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
    return float_param ("dB", -60.0, 0.0, 0.0, 0.1, 1.4);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:1<1,20,0.1>Ratio
    // Range: min:1.0, max:20.0, default: 1.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct ratio_tag {};

  void set (ratio_tag, float v)
  {
    // Original slider line: slider2:1<1,20,0.1>Ratio
    // Range: min:1.0, max:20.0, default: 1.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (ratio_tag)
  {
    // Original slider line: slider2:1<1,20,0.1>Ratio
    return float_param ("", 1.0, 20.0, 1.0, 0.01, 0.4);
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<-20,20,0.1>Gain
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
    // Original slider line: slider3:0<-20,20,0.1>Gain
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (gain_tag)
  {
    // Original slider line: slider3:0<-20,20,0.1>Gain
    return float_param ("dB", -20.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider4_slider (float v)
  {
    // Original slider line: slider4:2<0,3,1{Hard (Blown Capacitor),Soft (Blown
    // Capacitor),Hard,Soft}>Knee Range: min:0.0, max:3.0, default: 2.0,
    // step: 1.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct knee_tag {};

  void set (knee_tag, float v)
  {
    // Original slider line: slider4:2<0,3,1{Hard (Blown Capacitor),Soft (Blown
    // Capacitor),Hard,Soft}>Knee Range: min:0.0, max:3.0, default: 2.0,
    // step: 1.0
    v += 2;
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }

  static constexpr auto get_parameter (knee_tag)
  {
    // Original slider line: slider4:2<0,3,1{Hard (Blown Capacitor),Soft (Blown
    // Capacitor),Hard,Soft}>Knee
    return choice_param (0, make_cstr_array ("Hard", "Soft"));
  }

#endif

#if 0 // sidechain NA

#if 0
  void set_slider5_slider (float v)
  {
    // Original slider line: slider5:0<0,1,1{Normal,Sidechain}>Detector Input
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....

  struct slider5_tag {};

  void set (slider5_tag, float v)
  {
    // Original slider line: slider5:0<0,1,1{Normal,Sidechain}>Detector Input
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }

  static constexpr auto get_parameter (slider5_tag)
  {
    // Original slider line: slider5:0<0,1,1{Normal,Sidechain}>Detector Input
    return float_param ("", 0.0, 1.0, 0.0, 1.0);
  }

#endif

#endif // sidechain NA

#if 0
  void set_slider6_slider (float v)
  {
    // Original slider line: slider6:0<0,1,1{No,Yes}>Automatic Make-Up
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct agc_tag {};

  void set (agc_tag, float v)
  {
    // Original slider line: slider6:0<0,1,1{No,Yes}>Automatic Make-Up
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }

  static constexpr auto get_parameter (agc_tag)
  {
    // Original slider line: slider6:0<0,1,1{No,Yes}>Automatic Make-Up
    return choice_param (0, make_cstr_array ("No", "Yes"));
  }

#endif
#if 0
  void set_slider7_slider (float v)
  {
    // Original slider line: slider7:0<0,1,1{Peak,RMS}>Detection
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct detection_tag {};

  void set (detection_tag, float v)
  {
    // Original slider line: slider7:0<0,1,1{Peak,RMS}>Detection
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }

  static constexpr auto get_parameter (detection_tag)
  {
    // Original slider line: slider7:0<0,1,1{Peak,RMS}>Detection
    return choice_param (0, make_cstr_array ("Peak", "RMS"));
  }

#endif
#if 0
  void set_slider8_slider (float v)
  {
    // Original slider line: slider8:0<0,1,1{Feedforward,Feedback}>Detection
    // Source Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct detection_src_tag {};

  void set (detection_src_tag, float v)
  {
    // Original slider line: slider8:0<0,1,1{Feedforward,Feedback}>Detection
    // Source Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }

  static constexpr auto get_parameter (detection_src_tag)
  {
    // Original slider line: slider8:0<0,1,1{Feedforward,Feedback}>Detection
    // Source
    return choice_param (0, make_cstr_array ("Feedforward", "Feedback"));
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    threshold_tag,
    ratio_tag,
    gain_tag,
    knee_tag,
    agc_tag,
    detection_tag,
    detection_src_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double atcoef;
  double attime;
  double attimes;
  double automakeup;
  double cratio;
  double db2log;
  double fbacoef;
  double fbrcoef;
  double i;
  double log2db;
  double maxover;
  double overdb;
  double ratio;
  double relcoef;
  double reltime;
  double rundb;
  double sidechain;
  double slider1;
  double slider2;
  double slider3;
  double slider4;
  double slider5;
  double slider6;
  double slider7;
  double slider8;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    atcoef     = 0;
    attime     = 0;
    attimes    = 0;
    automakeup = 0;
    cratio     = 0;
    db2log     = 0;
    fbacoef    = 0;
    fbrcoef    = 0;
    i          = 0;
    log2db     = 0;
    maxover    = 0;
    overdb     = 0;
    ratio      = 0;
    relcoef    = 0;
    reltime    = 0;
    rundb      = 0;
    sidechain  = 0;
    slider1    = 0;
    slider2    = 0;
    slider3    = 0;
    slider4    = 0;
    slider5    = 0;
    slider6    = 0;
    slider7    = 0;
    slider8    = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double capsc;
  double cthreshv;
  double makeupv;
  double opto;
  double rmscoef;
  double rmsdet;
  double softknee;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    capsc    = 0;
    cthreshv = 0;
    makeupv  = 0;
    opto     = 0;
    rmscoef  = 0;
    rmsdet   = 0;
    softknee = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double ospl0;
  double ospl1;
  double runave;
  double runmax;
  double runospl;
  double spl2;
  double spl3;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    ospl0   = 0;
    ospl1   = 0;
    runave  = 0;
    runmax  = 0;
    runospl = 0;
    spl2    = 0;
    spl3    = 0;
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
    slider2 = 1.0;
    slider3 = 0.0;
    slider4 = 2.0;
    slider5 = 0.0;
    slider6 = 0.0;
    slider7 = 0.0;
    slider8 = 0.0;
    log2db  = 8.6858896380650365530225783783321;
    db2log  = 0.11512925464970228420089957273422;
    i       = 0.;
    heap_reset (attimes + 120);
    for (int $$i = 0, $$end = std::max (0, (int) (120.)); $$i < $$end; ++$$i) {
      heap (attimes + i)
        = ((0.08924 / i) + (0.60755 / eel2_pow (i, 2.)) - 0.00006);
      i += 1.;
    }
    attime     = 0.010;
    reltime    = 0.100;
    maxover    = 0.;
    ratio      = 0.;
    cratio     = 0.;
    rundb      = 0.;
    overdb     = 0.;
    maxover    = 0.;
    atcoef     = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
    relcoef    = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
    fbacoef    = std::exp (-1000. / (2. * jsfx_specialvar_get_srate()));
    fbrcoef    = std::exp (-1000. / (200. * jsfx_specialvar_get_srate()));
    sidechain  = 0.;
    automakeup = 0.;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double autogain = 0.;
    double cthresh  = 0.;
    double makeup   = 0.;
    double thresh   = 0.;
    double threshv  = 0.;
    thresh          = slider1;
    threshv         = std::exp (thresh * db2log);
    ratio           = slider2;
    softknee        = eel2_and (slider4, 1.);
    capsc           = [&] {
      if (eel2_and (slider4, 2.)) {
        return log2db;
      }
      else {
        return log2db * 2.08136898;
      }
    }();
    cthresh = [&] {
      if (softknee) {
        return (thresh - 3.);
      }
      else {
        return thresh;
      }
    }();
    cthreshv   = std::exp (cthresh * db2log);
    sidechain  = slider5;
    automakeup = slider6;
    if (automakeup) {
      (autogain
       = (std::abs (thresh) - (std::abs (thresh) / std::max (1., ratio - 1.)))
         / 2.);
    }
    else {
      (autogain = 0.);
    }
    makeup  = slider3;
    makeupv = std::exp ((makeup + autogain) * db2log);
    rmsdet  = slider7;
    if (rmsdet) {
      (rmscoef = std::exp (-1000. / (10. * jsfx_specialvar_get_srate())));
    }
    else {
      (rmscoef = std::exp (-1000. / (0.0025 * jsfx_specialvar_get_srate())));
    }
    opto = slider8;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double aspl0  = 0.;
    double aspl1  = 0.;
    double ave    = 0.;
    double det    = 0.;
    double gr     = 0.;
    double grv    = 0.;
    double maxspl = 0.;
    double ospl   = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];
      if (sidechain) {
        aspl0 = std::abs (spl2);
        aspl1 = std::abs (spl3);
      }
      else {
        if (opto) {
          ospl = ospl0 * ospl0 + ospl1 * ospl1;
          if (ospl > runospl) {
            (runospl = ospl + atcoef * (runospl - ospl));
          }
          else {
            (runospl = ospl + relcoef * (runospl - ospl));
          }
          ospl = std::sqrt (std::max (0., runospl));
          ospl *= 0.5;
          aspl0 = std::abs (ospl);
          aspl1 = std::abs (ospl);
        }
        else {
          aspl0 = std::abs (spl0);
          aspl1 = std::abs (spl1);
        }
      }
      if (rmsdet) {
        ave    = (aspl0 * aspl0) + (aspl1 * aspl1);
        runave = ave + rmscoef * (runave - ave);
        det    = std::sqrt (std::max (0., runave));
      }
      else {
        maxspl = std::max (aspl0, aspl1);
        maxspl *= maxspl;
        runave = maxspl + rmscoef * (runave - maxspl);
        det    = std::sqrt (std::max (0., runave));
      }
      overdb = std::log (det / cthreshv) * capsc;
      if (overdb > maxover) {
        maxover = overdb;
        attime
          = heap (attimes + (std::max (0., std::floor (std::abs (overdb)))));
        atcoef  = std::exp (-1. / (attime * jsfx_specialvar_get_srate()));
        reltime = overdb / 125.;
        relcoef = std::exp (-1. / (reltime * jsfx_specialvar_get_srate()));
      }
      overdb = std::max (0., overdb);
      if (overdb > rundb) {
        (rundb = overdb + atcoef * (rundb - overdb));
      }
      else {
        (rundb = overdb + relcoef * (rundb - overdb));
      }
      overdb = rundb;
      cratio = [&] {
        if (softknee) {
          return (1. + (ratio - 1.) * std::min (overdb, 6.) / 6.);
        }
        else {
          return ratio;
        }
      }();
      gr      = -overdb * (cratio - 1.) / cratio;
      grv     = std::exp (gr * db2log);
      runmax  = maxover + relcoef * (runmax - maxover);
      maxover = runmax;
      spl0 *= grv * makeupv;
      spl1 *= grv * makeupv;
      ospl0 = spl0;
      ospl1 = spl1;
      ;
    }
  }
}; /* jsfx_process */
}} // namespace artv::sstillwell
