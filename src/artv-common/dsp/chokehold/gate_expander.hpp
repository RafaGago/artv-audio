#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls

#include "artv-common/dsp/jsfx_engine/jsfx_engine.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace chokehold {
class gate_expander {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::dynamics;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }

  //----------------------------------------------------------------------------
  // stubs for JSFX special variables

  double jsfx_specialvar_get_srate() { return plugcontext->get_sample_rate(); }

  //----------------------------------------------------------------------------
public:
#if 0
  void set_operation_slider (float v)
  {
    // Original slider line: slider1:operation=0<0,1,{Gate  [fade to
    // silence],Expander  [fade to range]}> Operation Range: min:0.0, max:1.0,
    // default: 0.0, step: None
    if (v == operation) {
      return;
    }
    operation = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct operation_tag {};

  void set (operation_tag, int v)
  {
    // Original slider line: slider1:operation=0<0,1,{Gate  [fade to
    // silence],Expander  [fade to range]}> Operation Range: min:0.0, max:1.0,
    // default: 0.0, step: None
    if (v == (int) operation) {
      return;
    }
    operation = (float) v;
    slider();
  }

  static constexpr auto get_parameter (operation_tag)
  {
    // Original slider line: slider1:operation=0<0,1,{Gate  [fade to
    // silence],Expander  [fade to range]}> Operation
    return choice_param (1, make_cstr_array ("Gate", "Expander"));
  }

#endif
#if 0
  void set_gatethresh_slider (float v)
  {
    // Original slider line: slider2:gateThresh=-60<-60,0,0.1> Threshold [dB]
    // Range: min:-60.0, max:0.0, default: -60.0, step: 0.1
    if (v == gatethresh) {
      return;
    }
    gatethresh = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gatethresh_tag {};

  void set (gatethresh_tag, float v)
  {
    // Original slider line: slider2:gateThresh=-60<-60,0,0.1> Threshold [dB]
    // Range: min:-60.0, max:0.0, default: -60.0, step: 0.1
    if (v == gatethresh) {
      return;
    }
    gatethresh = v;
    slider();
  }

  static constexpr auto get_parameter (gatethresh_tag)
  {
    // Original slider line: slider2:gateThresh=-60<-60,0,0.1> Threshold [dB]
    return float_param ("dB", -60.0, 0.0, -60.0, 0.1);
  }

#endif
#if 0
  void set_gaterange_slider (float v)
  {
    // Original slider line: slider3:gateRange=-40<-40,0,0.1> Exp. Range [dB]
    // Range: min:-40.0, max:0.0, default: -40.0, step: 0.1
    if (v == gaterange) {
      return;
    }
    gaterange = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gaterange_tag {};

  void set (gaterange_tag, float v)
  {
    // Original slider line: slider3:gateRange=-40<-40,0,0.1> Exp. Range [dB]
    // Range: min:-40.0, max:0.0, default: -40.0, step: 0.1
    if (v == gaterange) {
      return;
    }
    gaterange = v;
    slider();
  }

  static constexpr auto get_parameter (gaterange_tag)
  {
    // Original slider line: slider3:gateRange=-40<-40,0,0.1> Exp. Range [dB]
    return float_param ("dB", -60.0, 0.0, -40.0, 0.01);
  }

#endif
#if 0
  void set_gatehyst_slider (float v)
  {
    // Original slider line: slider4:gateHyst=-3<-12,0,0.01> Hysteresis [dB]
    // Range: min:-12.0, max:0.0, default: -3.0, step: 0.01
    if (v == gatehyst) {
      return;
    }
    gatehyst = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gatehyst_tag {};

  void set (gatehyst_tag, float v)
  {
    // Original slider line: slider4:gateHyst=-3<-12,0,0.01> Hysteresis [dB]
    // Range: min:-12.0, max:0.0, default: -3.0, step: 0.01
    if (v == gatehyst) {
      return;
    }
    gatehyst = v;
    slider();
  }

  static constexpr auto get_parameter (gatehyst_tag)
  {
    // Original slider line: slider4:gateHyst=-3<-12,0,0.01> Hysteresis [dB]
    return float_param ("dB", -12.0, 0.0, -3.0, 0.01);
  }

#endif
#if 0
  void set_gateattack_slider (float v)
  {
    // Original slider line: slider5:gateAttack=5<1,50,0.1> Attack [ms]
    // Range: min:1.0, max:50.0, default: 5.0, step: 0.1
    if (v == gateattack) {
      return;
    }
    gateattack = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gateattack_tag {};

  void set (gateattack_tag, float v)
  {
    // Original slider line: slider5:gateAttack=5<1,50,0.1> Attack [ms]
    // Range: min:1.0, max:50.0, default: 5.0, step: 0.1
    if (v == gateattack) {
      return;
    }
    gateattack = v;
    slider();
  }

  static constexpr auto get_parameter (gateattack_tag)
  {
    // Original slider line: slider5:gateAttack=5<1,50,0.1> Attack [ms]
    return float_param ("ms", 1, 50.0, 5.0, 0.1);
  }

#endif
#if 0
  void set_gaterelease_slider (float v)
  {
    // Original slider line: slider6:gateRelease=300<50, 2500, 0.1> Release [ms]
    // Range: min:50.0, max:2500.0, default: 300.0, step: 0.1
    if (v == gaterelease) {
      return;
    }
    gaterelease = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gaterelease_tag {};

  void set (gaterelease_tag, float v)
  {
    // Original slider line: slider6:gateRelease=300<50, 2500, 0.1> Release [ms]
    // Range: min:50.0, max:2500.0, default: 300.0, step: 0.1
    if (v == gaterelease) {
      return;
    }
    gaterelease = v;
    slider();
  }

  static constexpr auto get_parameter (gaterelease_tag)
  {
    // Original slider line: slider6:gateRelease=300<50, 2500, 0.1> Release [ms]
    return float_param ("ms", 10.0, 2500.0, 300.0, 0.1);
  }

#endif
#if 0
  void set_linkamount_slider (float v)
  {
    // Original slider line: slider7:linkAmount=0<0,100,1> Stereo Link [%]
    // Range: min:0.0, max:100.0, default: 0.0, step: 1.0
    if (v == linkamount) {
      return;
    }
    linkamount = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct linkamount_tag {};

  void set (linkamount_tag, float v)
  {
    // Original slider line: slider7:linkAmount=0<0,100,1> Stereo Link [%]
    // Range: min:0.0, max:100.0, default: 0.0, step: 1.0
    if (v == linkamount) {
      return;
    }
    linkamount = v;
    slider();
  }

  static constexpr auto get_parameter (linkamount_tag)
  {
    // Original slider line: slider7:linkAmount=0<0,100,1> Stereo Link [%]
    return float_param ("%", 0.0, 100.0, 0.0, 1.0);
  }

#endif
#if 0
  void set_scfreq_slider (float v)
  {
    // Original slider line: slider8:scFreq=70<20,350,1> SC High Pass [Hz]
    // Range: min:20.0, max:350.0, default: 70.0, step: 1.0
    if (v == scfreq) {
      return;
    }
    scfreq = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct scfreq_tag {};

  void set (scfreq_tag, float v)
  {
    // Original slider line: slider8:scFreq=70<20,350,1> SC High Pass [Hz]
    // Range: min:20.0, max:350.0, default: 70.0, step: 1.0
    v = midi_note_to_hz (v);
    if (v == scfreq) {
      return;
    }
    scfreq = v;
    slider();
  }

  static constexpr auto get_parameter (scfreq_tag)
  {
    // Original slider line: slider8:scFreq=70<20,350,1> SC High Pass [Hz]
    return frequency_parameter (20.0, 350.0, 70.0);
  }

#endif
#if 0
  void set_routing_slider (float v)
  {
    // Original slider line: slider9:routing=0<0,3,{Stereo (internal SC),Stereo
    // (external SC 3-4),Mono (L),Mono (R),Mono (signal L / key R),Mono (key L /
    // signal R)}> Routing Range: min:0.0, max:3.0, default: 0.0, step: None
    if (v == routing) {
      return;
    }
    routing = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct routing_tag {};

  void set (routing_tag, int v)
  {
// Original slider line: slider9:routing=0<0,3,{Stereo (internal SC),Stereo
// (external SC 3-4),Mono (L),Mono (R),Mono (signal L / key R),Mono (key L /
// signal R)}> Routing Range: min:0.0, max:3.0, default: 0.0, step: None
#if 1
    v += v != 0; // skip sidechain
#endif
    if (v == (int) routing) {
      return;
    }
    routing = (float) v;
    slider();
  }

  static constexpr auto get_parameter (routing_tag)
  {
    // Original slider line: slider9:routing=0<0,3,{Stereo (internal SC),Stereo
    // (external SC 3-4),Mono (L),Mono (R),Mono (signal L / key R),Mono (key L /
    // signal R)}> Routing
    return choice_param (
      0,
      make_cstr_array (
        "LR",
#if 0
        "Sidechain",
#endif
        "L",
        "R",
        "Mono (Sig L, Key R)",
        "Mono (Sig R, Key L)"));
  }

#endif
#if 0
#else
  // Snippet for parameter boilerplate in the authors framework....
  using parameters = mp_list<
    operation_tag,
    gatethresh_tag,
    gaterange_tag,
    gatehyst_tag,
    gateattack_tag,
    gaterelease_tag,
    linkamount_tag,
    scfreq_tag,
    routing_tag>;
#endif

private:
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double gateattack;
  double gatehyst;
  double gaterange;
  double gaterelease;
  double gatethresh;
  double halfpi;
  double keyl;
  double keyr;
  double linkamount;
  double msrate;
  double operation;
  double routing;
  double scfreq;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    gateattack  = 0;
    gatehyst    = 0;
    gaterange   = 0;
    gaterelease = 0;
    gatethresh  = 0;
    halfpi      = 0;
    keyl        = 0;
    keyr        = 0;
    linkamount  = 0;
    msrate      = 0;
    operation   = 1;
    routing     = 0;
    scfreq      = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double filterl$a1;
  double filterl$a2;
  double filterl$a3;
  double filterl$m0;
  double filterl$m1;
  double filterl$m2;
  double filterr$a1;
  double filterr$a2;
  double filterr$a3;
  double filterr$m0;
  double filterr$m1;
  double filterr$m2;
  double gatel$attack;
  double gatel$hysteresis;
  double gatel$range;
  double gatel$release;
  double gatel$threshold;
  double gater$attack;
  double gater$hysteresis;
  double gater$range;
  double gater$release;
  double gater$threshold;
  double lnkmix;
  double splmix;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    filterl$a1       = 0;
    filterl$a2       = 0;
    filterl$a3       = 0;
    filterl$m0       = 0;
    filterl$m1       = 0;
    filterl$m2       = 0;
    filterr$a1       = 0;
    filterr$a2       = 0;
    filterr$a3       = 0;
    filterr$m0       = 0;
    filterr$m1       = 0;
    filterr$m2       = 0;
    gatel$attack     = 0;
    gatel$hysteresis = 0;
    gatel$range      = 0;
    gatel$release    = 0;
    gatel$threshold  = 0;
    gater$attack     = 0;
    gater$hysteresis = 0;
    gater$range      = 0;
    gater$release    = 0;
    gater$threshold  = 0;
    lnkmix           = 0;
    splmix           = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double filterl$ic1eq;
  double filterl$ic2eq;
  double filterl$v1;
  double filterl$v2;
  double filterl$v3;
  double filterr$ic1eq;
  double filterr$ic2eq;
  double filterr$v1;
  double filterr$v2;
  double filterr$v3;
  double gatel$envelope;
  double gatel$gain;
  double gater$envelope;
  double gater$gain;
  double spl2;
  double spl3;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    filterl$ic1eq  = 0;
    filterl$ic2eq  = 0;
    filterl$v1     = 0;
    filterl$v2     = 0;
    filterl$v3     = 0;
    filterr$ic1eq  = 0;
    filterr$ic2eq  = 0;
    filterr$v1     = 0;
    filterr$v2     = 0;
    filterr$v3     = 0;
    gatel$envelope = 0;
    gatel$gain     = 0;
    gater$envelope = 0;
    gater$gain     = 0;
    spl2           = 0;
    spl3           = 0;
  }
  //----------------------------------------------------------------------------
  plugin_context* plugcontext;

public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    operation   = 0.0;
    gatethresh  = -60.0;
    gaterange   = -40.0;
    gatehyst    = -3.0;
    gateattack  = 5.0;
    gaterelease = 300.0;
    linkamount  = 0.0;
    scfreq      = 70.0;
    routing     = 0.0;
    halfpi      = 3.141592653589793 * 0.5;
    msrate      = jsfx_specialvar_get_srate() * 0.001;
    keyl        = 0.0;
    keyr        = 0.0;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    init$gatesetup (
      gateattack,
      gaterelease,
      gatethresh,
      gaterange,
      gatehyst,
      gatel$range,
      gatel$attack,
      gatel$release,
      gatel$hysteresis,
      gatel$threshold);
    init$gatesetup (
      gateattack,
      gaterelease,
      gatethresh,
      gaterange,
      gatehyst,
      gater$range,
      gater$attack,
      gater$release,
      gater$hysteresis,
      gater$threshold);
    lnkmix = linkamount * 0.01;
    splmix = 1.0 - lnkmix;
    init$eqhp (
      scfreq,
      1.5,
      filterl$a1,
      filterl$a2,
      filterl$a3,
      filterl$m0,
      filterl$m1,
      filterl$m2);
    init$eqhp (
      scfreq,
      1.5,
      filterr$a1,
      filterr$a2,
      filterr$a3,
      filterr$m0,
      filterr$m1,
      filterr$m2);
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    double linked = 0.;

    for (uint $$i = 0; $$i < samples; ++$$i) {
      auto& spl0 = chnls[0][$$i];
      auto& spl1 = chnls[1][$$i];
      switch (uint (routing)) {
      case 0:
        keyl = init$eqtick (
          spl0,
          filterl$v3,
          filterl$ic2eq,
          filterl$v1,
          filterl$a1,
          filterl$ic1eq,
          filterl$a2,
          filterl$v2,
          filterl$a3,
          filterl$m0,
          filterl$m1,
          filterl$m2);
        keyr = init$eqtick (
          spl1,
          filterr$v3,
          filterr$ic2eq,
          filterr$v1,
          filterr$a1,
          filterr$ic1eq,
          filterr$a2,
          filterr$v2,
          filterr$a3,
          filterr$m0,
          filterr$m1,
          filterr$m2);
        if ((lnkmix > 0.)) {
          linked = std::max (keyl, keyr) * lnkmix;
          keyl *= splmix;
          keyr *= splmix;
          keyl += linked;
          keyr += linked;
        }
        spl0 *= init$tickgate (
          keyl,
          gatel$threshold,
          gatel$hysteresis,
          gatel$range,
          gatel$gain,
          gatel$envelope,
          gatel$attack,
          gatel$release);
        spl1 *= init$tickgate (
          keyr,
          gater$threshold,
          gater$hysteresis,
          gater$range,
          gater$gain,
          gater$envelope,
          gater$attack,
          gater$release);
        break;
      case 1:
        keyl = init$eqtick (
          spl2,
          filterl$v3,
          filterl$ic2eq,
          filterl$v1,
          filterl$a1,
          filterl$ic1eq,
          filterl$a2,
          filterl$v2,
          filterl$a3,
          filterl$m0,
          filterl$m1,
          filterl$m2);
        keyr = init$eqtick (
          spl3,
          filterr$v3,
          filterr$ic2eq,
          filterr$v1,
          filterr$a1,
          filterr$ic1eq,
          filterr$a2,
          filterr$v2,
          filterr$a3,
          filterr$m0,
          filterr$m1,
          filterr$m2);
        if (lnkmix > 0.) {
          linked = std::max (keyl, keyr) * lnkmix;
          keyl *= splmix;
          keyr *= splmix;
          keyl += linked;
          keyr += linked;
        }
        spl0 *= init$tickgate (
          keyl,
          gatel$threshold,
          gatel$hysteresis,
          gatel$range,
          gatel$gain,
          gatel$envelope,
          gatel$attack,
          gatel$release);
        spl1 *= init$tickgate (
          keyr,
          gater$threshold,
          gater$hysteresis,
          gater$range,
          gater$gain,
          gater$envelope,
          gater$attack,
          gater$release);
        break;
      case 2:
        (spl1 = spl0
         = spl0
           * init$tickgate (
             init$eqtick (
               spl0,
               filterl$v3,
               filterl$ic2eq,
               filterl$v1,
               filterl$a1,
               filterl$ic1eq,
               filterl$a2,
               filterl$v2,
               filterl$a3,
               filterl$m0,
               filterl$m1,
               filterl$m2),
             gatel$threshold,
             gatel$hysteresis,
             gatel$range,
             gatel$gain,
             gatel$envelope,
             gatel$attack,
             gatel$release));
        break;
      case 3:

        (spl0 = spl1
         = spl1
           * init$tickgate (
             init$eqtick (
               spl1,
               filterr$v3,
               filterr$ic2eq,
               filterr$v1,
               filterr$a1,
               filterr$ic1eq,
               filterr$a2,
               filterr$v2,
               filterr$a3,
               filterr$m0,
               filterr$m1,
               filterr$m2),
             gater$threshold,
             gater$hysteresis,
             gater$range,
             gater$gain,
             gater$envelope,
             gater$attack,
             gater$release));
        break;
      case 4:
        (spl1 = spl0
         = spl0
           * init$tickgate (
             init$eqtick (
               spl1,
               filterl$v3,
               filterl$ic2eq,
               filterl$v1,
               filterl$a1,
               filterl$ic1eq,
               filterl$a2,
               filterl$v2,
               filterl$a3,
               filterl$m0,
               filterl$m1,
               filterl$m2),
             gatel$threshold,
             gatel$hysteresis,
             gatel$range,
             gatel$gain,
             gatel$envelope,
             gatel$attack,
             gatel$release));
        break;
      case 5:
        (spl0 = spl1
         = spl1
           * init$tickgate (
             init$eqtick (
               spl0,
               filterr$v3,
               filterr$ic2eq,
               filterr$v1,
               filterr$a1,
               filterr$ic1eq,
               filterr$a2,
               filterr$v2,
               filterr$a3,
               filterr$m0,
               filterr$m1,
               filterr$m2),
             gater$threshold,
             gater$hysteresis,
             gater$range,
             gater$gain,
             gater$envelope,
             gater$attack,
             gater$release));
        break;
      default:
        break;
      }
    }
  }

private:
  // functions for section "init"
  //----------------------------------------------------------------------------
  double init$dbtogain (double decibels)
  {
    return eel2_pow (10.0, (decibels / 20.0));
  }
  //----------------------------------------------------------------------------
  double init$envfollow (
    double  signal,
    double& $envelope,
    double& $attack,
    double& $release)
  {
    $envelope
      = (((signal > $envelope) * $attack) + ((signal <= $envelope) * $release))
        * ($envelope - signal)
      + signal;
    return $envelope;
  }
  //----------------------------------------------------------------------------
  double init$envsetup (
    double  msattack,
    double  msrelease,
    double  dbthreshold,
    double  dbhysteresis,
    double& $attack,
    double& $release,
    double& $hysteresis,
    double& $threshold)
  {
    $attack = std::exp (
      -1000. * init$fastreciprocal (msattack * jsfx_specialvar_get_srate()));
    $release = std::exp (
      -1000. * init$fastreciprocal (msrelease * jsfx_specialvar_get_srate()));
    $hysteresis = init$dbtogain (dbhysteresis);
    $threshold  = init$dbtogain (dbthreshold);
    return $threshold;
  }
  //----------------------------------------------------------------------------
  double init$eqhp (
    double  hz,
    double  q,
    double& $a1,
    double& $a2,
    double& $a3,
    double& $m0,
    double& $m1,
    double& $m2)
  {
    double k = 0.;
    double g = 0.;
    g        = std::tan (halfpi * (hz / jsfx_specialvar_get_srate()));
    k        = 1.0 / q;
    $a1      = 1.0 / (1.0 + g * (g + k));
    $a2      = $a1 * g;
    $a3      = $a2 * g;
    $m0      = 1.0;
    $m1      = -k;
    $m2      = -1.0;
    return $m2;
  }
  //----------------------------------------------------------------------------
  double init$eqtick (
    double  sample,
    double& $v3,
    double& $ic2eq,
    double& $v1,
    double& $a1,
    double& $ic1eq,
    double& $a2,
    double& $v2,
    double& $a3,
    double& $m0,
    double& $m1,
    double& $m2)
  {
    $v3    = sample - $ic2eq;
    $v1    = $a1 * $ic1eq + $a2 * $v3;
    $v2    = $ic2eq + $a2 * $ic1eq + $a3 * $v3;
    $ic1eq = 2.0 * $v1 - $ic1eq;
    $ic2eq = 2.0 * $v2 - $ic2eq;
    return ($m0 * sample + $m1 * $v1 + $m2 * $v2);
  }
#if 1
  // Was missing, typical quake impl
  double invsqrt (double number)
  {
    double       y  = number;
    double       x2 = y * 0.5;
    std::int64_t i  = *(std::int64_t*) &y;
    // The magic number is for doubles is from
    // https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double*) &i;
    y = y * (1.5 - (x2 * y * y)); // 1st iteration
    y = y * (1.5 - (x2 * y * y)); // 2nd iteration, this can be removed
    return y;
  }
#endif
  //----------------------------------------------------------------------------
  double init$fastreciprocal (double value)
  {
    double v = invsqrt (value);
    return v * v;
  }
  //----------------------------------------------------------------------------
  double init$gatesetup (
    double  msattack,
    double  msrelease,
    double  dbthreshold,
    double  dbrange,
    double  dbhysteresis,
    double& $range,
    double& $attack,
    double& $release,
    double& $hysteresis,
    double& $threshold)
  {
    init$envsetup (
      msattack,
      msrelease,
      dbthreshold,
      dbhysteresis,
      $attack,
      $release,
      $hysteresis,
      $threshold);
    $range = [&] {
      if (eel2_eq (operation, 1.)) {
        return init$dbtogain (dbrange);
      }
      else {
        return 0.;
      }
    }();
    return $range;
  }
  //----------------------------------------------------------------------------
  double init$tickgate (
    double  sample,
    double& $threshold,
    double& $hysteresis,
    double& $range,
    double& $gain,
    double& $envelope,
    double& $attack,
    double& $release)
  {
    init$envfollow (
      [&] {
        if ((std::abs (sample) < ($threshold * $hysteresis))) {
          return $range;
        }
        else {
          return 1.;
        }
      }(),
      $envelope,
      $attack,
      $release);
    $gain = std::max ($envelope, $range);
    return $gain;
  }
}; /* jsfx_process */
}} // namespace artv::chokehold
