// Ported from https://github.com/JoepVanlier/JSFX.git
// commit sha: 5d54a165b806b2773792623376f6a51cc957368b
#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
// slider1:Delay=8<1,25,1>-Delay [ms]
// slider2:Strength=0.3<0,1,.0001>-Strength [-]
// slider3:LowCut=0.4<0,1,.0001>-Crossover [log(w)]
// slider4:SideLevel=1<0,2,.0001>-Blend old side level [-]
// slider5:CheckMono=0<0,1,1{No, Yes}>-Check Mono?
// slider6:SideHP=0<0,1,1{No,Yes}>-Pass side through HPF?
// slider7:VibratoSpeed=15<0,30,.1>-Vibrato Speed [-]
// slider8:VibratoAmount=0<0,40,.1>-Vibrato Amount [-]
// slider9:Saturation=0<0,1,0.0000001>-Non-linearity [-]
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace saike {

struct stereo_bub3 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::stereo;
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
#else
  std::vector<float> mempool;
#endif
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  void jsfx_memset (double idx, double val, double size)
  {
    std::memset (&mempool[(size_t) idx], (int) val, (size_t) size);
  }
  double jsfx_sign (double value)
  {
    auto v = *((uint64_t*) ((void*) &value));
    return (v == 0) ? 0. : (v & (1ull << 63)) ? -1. : 1.;
  }

  //----------------------------------------------------------------------------
  // stubs for JSFX special variables
  float  sample_rate = 44100;
  double jsfx_specialvar_get_srate() { return sample_rate; }

  //----------------------------------------------------------------------------
  // stubs for sliders
public:
#if 1
  double get_slider_checkmono()
  {
    // TODO: stub, add code for getting "checkmono"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    // Original line: slider5:CheckMono=0<0,1,1{No, Yes}>-Check Mono?
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_checkmono()
  {
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    return checkmonop;
  }

  float checkmonop = 0.0;

  struct checkmono_tag {};

  void set (checkmono_tag, float v)
  {
    if (v == checkmonop) {
      return;
    }
    checkmonop = v;
    slider();
  }

  static constexpr auto get_parameter (checkmono_tag)
  {
    // Original slider line: slider5:CheckMono=0<0,1,1{No, Yes}>-Check Mono?
    return float_param ("", 0.0, 1.0, 0.0, 1.0);
  }
#endif
#if 0
  double get_slider_delay()
  {
    // TODO: stub, add code for getting "delay"
    // Range: min:1.0, max:25.0, default: 8.0, step: 1.0
    // Original line: slider1:Delay=8<1,25,1>-Delay [ms]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_delay()
  {
    // Range: min:1.0, max:25.0, default: 8.0, step: 1.0
    return delayp;
  }

  float delayp = 8.0;

  struct delay_tag {};

  void set (delay_tag, float v)
  {
    if (v == delayp) {
      return;
    }
    delayp = v;
    slider();
  }

  static constexpr auto get_parameter (delay_tag)
  {
    // Original slider line: slider1:Delay=8<1,25,1>-Delay [ms]
    return float_param ("ms", 1.0, 25.0, 8.0, 1.0);
  }
#endif
#if 0
  double get_slider_lowcut()
  {
    // TODO: stub, add code for getting "lowcut"
    // Range: min:0.0, max:1.0, default: 0.4, step: 0.0001
    // Original line: slider3:LowCut=0.4<0,1,.0001>-Crossover [log(w)]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_lowcut()
  {
    // Range: min:0.0, max:1.0, default: 0.4, step: 0.0001
    return lowcutp;
  }

  float lowcutp = 0.4;

  struct lowcut_tag {};

  void set (lowcut_tag, float v)
  {
    if (v == lowcutp) {
      return;
    }
    lowcutp = v;
    slider();
  }

  static constexpr auto get_parameter (lowcut_tag)
  {
    // Original slider line: slider3:LowCut=0.4<0,1,.0001>-Crossover [log(w)]
    return float_param ("log(w)", 0.0, 1.0, 0.4, 0.0001);
  }
#endif
#if 0
  double get_slider_saturation()
  {
    // TODO: stub, add code for getting "saturation"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1e-07
    // Original line: slider9:Saturation=0<0,1,0.0000001>-Non-linearity [-]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_saturation()
  {
    // Range: min:0.0, max:1.0, default: 0.0, step: 1e-07
    return saturationp;
  }

  float saturationp = 0.0;

  struct saturation_tag {};

  void set (saturation_tag, float v)
  {
    if (v == saturationp) {
      return;
    }
    saturationp = v;
    slider();
  }

  static constexpr auto get_parameter (saturation_tag)
  {
    // Original slider line: slider9:Saturation=0<0,1,0.0000001>-Non-linearity
    // [-]
    return float_param ("", 0.0, 1.0, 0.0, 1e-07);
  }
#endif
#if 0
  double get_slider_sidehp()
  {
    // TODO: stub, add code for getting "sidehp"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    // Original line: slider6:SideHP=0<0,1,1{No,Yes}>-Pass side through HPF?
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_sidehp()
  {
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    return sidehpp;
  }

  float sidehpp = 0.0;

  struct sidehp_tag {};

  void set (sidehp_tag, float v)
  {
    if (v == sidehpp) {
      return;
    }
    sidehpp = v;
    slider();
  }

  static constexpr auto get_parameter (sidehp_tag)
  {
    // Original slider line: slider6:SideHP=0<0,1,1{No,Yes}>-Pass side through
    // HPF?
    return choice_param (0, make_cstr_array ("No", "Yes"));
  }
#endif
#if 0
  double get_slider_sidelevel()
  {
    // TODO: stub, add code for getting "sidelevel"
    // Range: min:0.0, max:2.0, default: 1.0, step: 0.0001
    // Original line: slider4:SideLevel=1<0,2,.0001>-Blend old side level [-]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_sidelevel()
  {
    // Range: min:0.0, max:2.0, default: 1.0, step: 0.0001
    return sidelevelp;
  }

  float sidelevelp = 1.0;

  struct sidelevel_tag {};

  void set (sidelevel_tag, float v)
  {
    if (v == sidelevelp) {
      return;
    }
    sidelevelp = v;
    slider();
  }

  static constexpr auto get_parameter (sidelevel_tag)
  {
    // Original slider line: slider4:SideLevel=1<0,2,.0001>-Blend old side level
    // [-]
    return float_param ("", 0.0, 2.0, 1.0, 0.0001);
  }
#endif
#if 0
  double get_slider_strength()
  {
    // TODO: stub, add code for getting "strength"
    // Range: min:0.0, max:1.0, default: 0.3, step: 0.0001
    // Original line: slider2:Strength=0.3<0,1,.0001>-Strength [-]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_strength()
  {
    // Range: min:0.0, max:1.0, default: 0.3, step: 0.0001
    return strengthp;
  }

  float strengthp = 0.3;

  struct strength_tag {};

  void set (strength_tag, float v)
  {
    if (v == strengthp) {
      return;
    }
    strengthp = v;
    slider();
  }

  static constexpr auto get_parameter (strength_tag)
  {
    // Original slider line: slider2:Strength=0.3<0,1,.0001>-Strength [-]
    return float_param ("", 0.0, 1.0, 0.3, 0.0001);
  }
#endif
#if 0
  double get_slider_vibratoamount()
  {
    // TODO: stub, add code for getting "vibratoamount"
    // Range: min:0.0, max:40.0, default: 0.0, step: 0.1
    // Original line: slider8:VibratoAmount=0<0,40,.1>-Vibrato Amount [-]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_vibratoamount()
  {
    // Range: min:0.0, max:40.0, default: 0.0, step: 0.1
    return vibratoamountp;
  }

  float vibratoamountp = 0.0;

  struct vibratoamount_tag {};

  void set (vibratoamount_tag, float v)
  {
    if (v == vibratoamountp) {
      return;
    }
    vibratoamountp = v;
    slider();
  }

  static constexpr auto get_parameter (vibratoamount_tag)
  {
    // Original slider line: slider8:VibratoAmount=0<0,40,.1>-Vibrato Amount [-]
    return float_param ("", 0.0, 40.0, 0.0, 0.1);
  }
#endif
#if 0
  double get_slider_vibratospeed()
  {
    // TODO: stub, add code for getting "vibratospeed"
    // Range: min:0.0, max:30.0, default: 15.0, step: 0.1
    // Original line: slider7:VibratoSpeed=15<0,30,.1>-Vibrato Speed [-]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_vibratospeed()
  {
    // Range: min:0.0, max:30.0, default: 15.0, step: 0.1
    return vibratospeedp;
  }

  float vibratospeedp = 15.0;

  struct vibratospeed_tag {};

  void set (vibratospeed_tag, float v)
  {
    if (v == vibratospeedp) {
      return;
    }
    vibratospeedp = v;
    slider();
  }

  static constexpr auto get_parameter (vibratospeed_tag)
  {
    // Original slider line: slider7:VibratoSpeed=15<0,30,.1>-Vibrato Speed [-]
    return float_param ("", 0.0, 30.0, 15.0, 0.1);
  }
#endif

  using parameters = mp_list<
    delay_tag,
    lowcut_tag,
    sidehp_tag,
    strength_tag,
    sidelevel_tag,
    saturation_tag,
    vibratoamount_tag,
    vibratospeed_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double buffer$scopebuffer;
  double buffer$scopebuffermax;
  double buffer$scopeptr;
  double bufferdist;
  double chdelay;
  double delaybuf1;
  double delta_t;
#if 0
  double displayfreq;
  double fontface;
  double gfx_a;
  double gfx_x;
  double gfx_y;
  double gonioout$recptr;
  double gonioout$rend;
  double gonioout$rstart;
  double gonioout$samples;
  double gonioout$thisui;
#endif
  double hp$a1;
  double hp$a2;
  double hp$a3;
  double hp$k;
  double hp2$a1;
  double hp2$a2;
  double hp2$a3;
  double hp2$k;
  double lowcutfreq;
  double maxbuffersize;
  double maxdelay;
  double maxlowcut;
  double maxsidelevel;
  double maxstrength;
  double maxvibratoamount;
  double maxvibratospeed;
  double mindelay;
  double minlowcut;
  double minsidelevel;
  double minstrength;
#if 0
  double newui;
  double scaling;
#endif
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    buffer$scopebuffer    = 0;
    buffer$scopebuffermax = 0;
    buffer$scopeptr       = 0;
    bufferdist            = 0;
    chdelay               = 0;
    delaybuf1             = 0;
    delta_t               = 0;
#if 0
    displayfreq           = 0;
    fontface              = 0;
    gfx_a                 = 0;
    gfx_x                 = 0;
    gfx_y                 = 0;
    gonioout$recptr       = 0;
    gonioout$rend         = 0;
    gonioout$rstart       = 0;
    gonioout$samples      = 0;
    gonioout$thisui       = 0;
#endif
    hp$a1            = 0;
    hp$a2            = 0;
    hp$a3            = 0;
    hp$k             = 0;
    hp2$a1           = 0;
    hp2$a2           = 0;
    hp2$a3           = 0;
    hp2$k            = 0;
    lowcutfreq       = 0;
    maxbuffersize    = 0;
    maxdelay         = 0;
    maxlowcut        = 0;
    maxsidelevel     = 0;
    maxstrength      = 0;
    maxvibratoamount = 0;
    maxvibratospeed  = 0;
    mindelay         = 0;
    minlowcut        = 0;
    minsidelevel     = 0;
    minstrength      = 0;
#if 0
    newui            = 0;
    scaling          = 0;
#endif
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double buffer$frac;
  double buffer$readptr;
  double ct;
  double forceupdate;
  double hp$ic10eq;
  double hp$ic11eq;
  double hp$ic12eq;
  double hp$ic1eq;
  double hp$ic2eq;
  double hp$ic3eq;
  double hp$ic4eq;
  double hp$ic5eq;
  double hp$ic6eq;
  double hp$ic7eq;
  double hp$ic8eq;
  double hp$ic9eq;
  double hp2$ic10eq;
  double hp2$ic11eq;
  double hp2$ic12eq;
  double hp2$ic1eq;
  double hp2$ic2eq;
  double hp2$ic3eq;
  double hp2$ic4eq;
  double hp2$ic5eq;
  double hp2$ic6eq;
  double hp2$ic7eq;
  double hp2$ic8eq;
  double hp2$ic9eq;
  double tdelay;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    buffer$frac    = 0;
    buffer$readptr = 0;
    ct             = 0;
    forceupdate    = 0;
    hp$ic10eq      = 0;
    hp$ic11eq      = 0;
    hp$ic12eq      = 0;
    hp$ic1eq       = 0;
    hp$ic2eq       = 0;
    hp$ic3eq       = 0;
    hp$ic4eq       = 0;
    hp$ic5eq       = 0;
    hp$ic6eq       = 0;
    hp$ic7eq       = 0;
    hp$ic8eq       = 0;
    hp$ic9eq       = 0;
    hp2$ic10eq     = 0;
    hp2$ic11eq     = 0;
    hp2$ic12eq     = 0;
    hp2$ic1eq      = 0;
    hp2$ic2eq      = 0;
    hp2$ic3eq      = 0;
    hp2$ic4eq      = 0;
    hp2$ic5eq      = 0;
    hp2$ic6eq      = 0;
    hp2$ic7eq      = 0;
    hp2$ic8eq      = 0;
    hp2$ic9eq      = 0;
    tdelay         = 0;
  }
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    init_init_variables();
    init_block_variables();
    sample_rate      = pc.get_sample_rate();
    delta_t          = 44. / sample_rate;
    mindelay         = 1.;
    maxdelay         = 25.;
    minstrength      = 0.;
    maxstrength      = 1.;
    maxvibratoamount = 40.;
    maxvibratospeed  = 30.;
    minlowcut        = 0.;
    maxlowcut        = 1.;
    minsidelevel     = 0.;
    maxsidelevel     = 2.;
#if 0
    fontface         = 0.; /* jsfx2cpp strings unsupported: was: "Arial" */
#endif

#if 0
    bufferDist = 65536.;
#else
    bufferdist = (sample_rate * 1000.f) / get_parameter (delay_tag {}).max;
    mempool.resize (bufferdist);
    memset (mempool.data(), 0, sizeof mempool[0] * mempool.size());
#endif
    delaybuf1 = 0.;
    init$initbuffer (
      delaybuf1,
      delaybuf1 + bufferdist - 1.,
      buffer$scopebuffer,
      buffer$scopebuffermax,
      buffer$scopeptr);
#if 0
    init$initgoniometer (
      bufferdist,
      1000.,
      gonioout$samples,
      gonioout$rstart,
      gonioout$rend,
      gonioout$recptr,
      gonioout$thisui);
#endif
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider() { init$sliderupdate(); }

public:
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double avg  = 0.;
    double pref = 0.;
    double rb   = 0.;
    double side = 0.;

    delta_t = 44. / jsfx_specialvar_get_srate();

    for (int $$i = 0, $$end = samples; $$i < $$end; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];
      ct += get_slider_vibratospeed() / jsfx_specialvar_get_srate();
      tdelay += chdelay * delta_t - delta_t * tdelay;
      init$setoffset (
        tdelay + get_slider_vibratoamount() * (.5 + .5 * std::sin (ct)),
        buffer$readptr,
        buffer$scopeptr,
        buffer$frac,
        buffer$scopebuffer,
        buffer$scopebuffermax);
      if (eel2_eq (forceupdate, 1.)) {
        init$sliderupdate();
        forceupdate = 0.;
      };
      avg  = .5 * (spl0 + spl1);
      side = .5 * (spl0 - spl1);
      init$updatebuffer (
        avg, buffer$scopeptr, buffer$scopebuffermax, buffer$scopebuffer);
      rb = init$readbuffer (
        buffer$readptr, buffer$scopebuffermax, buffer$scopebuffer, buffer$frac);
      rb = block$eval_linearsvf_hp (
        rb,
        hp$ic2eq,
        hp$a1,
        hp$ic1eq,
        hp$a2,
        hp$a3,
        hp$k,
        hp$ic4eq,
        hp$ic3eq,
        hp$ic6eq,
        hp$ic5eq,
        hp$ic8eq,
        hp$ic7eq,
        hp$ic10eq,
        hp$ic9eq,
        hp$ic12eq,
        hp$ic11eq);
      if (get_slider_sidehp()) {
        (side = block$eval_linearsvf_hp (
           side,
           hp2$ic2eq,
           hp2$a1,
           hp2$ic1eq,
           hp2$a2,
           hp2$a3,
           hp2$k,
           hp2$ic4eq,
           hp2$ic3eq,
           hp2$ic6eq,
           hp2$ic5eq,
           hp2$ic8eq,
           hp2$ic7eq,
           hp2$ic10eq,
           hp2$ic9eq,
           hp2$ic12eq,
           hp2$ic11eq));
      }
      if (get_slider_saturation()) {
        pref = std::exp (get_slider_saturation() * 9. - 6.);
        rb   = init$tanh (pref * rb) / pref;
      }
      spl0 = avg + get_slider_strength() * rb + get_slider_sidelevel() * side;
      spl1 = avg - get_slider_strength() * rb - get_slider_sidelevel() * side;
      spl1 = .5 * (spl0 + spl1);
      if (eel2_eq (get_slider_checkmono(), 1.)) {
        (spl0 = spl1);
      }
#if 0
      init$feedgonio (
        spl0, spl1, gonioout$recptr, gonioout$rend, gonioout$rstart);
      ;
#endif
    }
  }
  // functions for section "init"
private:
  //----------------------------------------------------------------------------
  double init$clearbuffer (double& scopebuffer, double& scopeptr)
  {
    jsfx_memset (scopebuffer, 0., maxbuffersize);
    scopeptr = scopebuffer;
    return scopeptr;
  }
//----------------------------------------------------------------------------
#if 0
  double init$drawgoniometer (
    double  r,
    double  g,
    double  b,
    double  a,
    double& recptr,
    double& x,
    double& w,
    double& y,
    double& h,
    double& samples,
    double& dgonio_zoom,
    double& rend,
    double& rstart)
  {
    double yy     = 0.;
    double angle  = 0.;
    double horiz  = 0.;
    double radius = 0.;
    double step   = 0.;
    double sign0  = 0.;
    double xx     = 0.;
    double cy     = 0.;
    double cx     = 0.;
    double s0     = 0.;
    double s1     = 0.;
    double sign1  = 0.;
    double vert   = 0.;
    double ptr    = 0.;
    ptr           = recptr + 2.;
    cx            = x + 0.5 * w;
    cy            = y + 0.5 * h;
    gfx_set (r, g, b, a) /*TODO: unknown call */;
    return [&] {
      double $$loop_ret_0 = {};
      for (int $$i = 0, $$end = std::max (0, (int) (samples)); $$i < $$end;
           ++$$i) {
        s0 = mempool[(size_t) (ptr)];
        (ptr += 1.);
        s1    = mempool[(size_t) (ptr)];
        sign0 = jsfx_sign (s0);
        sign1 = jsfx_sign (s1);
        angle = std::atan (s0 / s1);
        if (
          (eel2_eq (sign0, 1.) && eel2_eq (sign1, -1.))
          || (eel2_eq (sign0, -1.) && eel2_eq (sign1, -1.))) {
          angle += 3.141592654;
        }
        if ((eel2_eq (sign0, -1.) && eel2_eq (sign1, 1.))) {
          angle += 6.283185307;
        }
        angle = 4.71238898;
        if (eel2_eq (s1, 0.)) {
          if (s0 > 0.) {
            angle = 1.570796327;
          }
          else {
            angle;
          }
        }
        angle = 3.141592654;
        if (eel2_eq (s0, 0.)) {
          if (s1 > 0.) {
            angle = 0.;
          }
          else {
            angle;
          }
        }
        radius = (1. + dgonio_zoom)
          * sqrt (std::sqrt (s0) + std::sqrt (s1)) /*TODO: unknown call */;
        if (radius > 1.) {
          radius = 1.;
        }
        angle -= .25 * 3.141592653589793;
        xx    = cx - .5 * h * std::sin (angle) * radius;
        yy    = cy - .5 * h * std::cos (angle) * radius;
        gfx_a = .1 * a;
        gfx_circle (xx, yy, 1., 1.) /*TODO: unknown call */;
        ptr += 1.;
        $$loop_ret_0 = [&] {
          if (ptr > rend) {
            ptr = rstart;
            return ptr;
          }
          else {
            return 0. / 0.;
          }
        }();
      }
      return $$loop_ret_0;
    }();
  }

  //----------------------------------------------------------------------------
  double init$drawgoniotop (
    double  r,
    double  g,
    double  b,
    double  a,
    double& x,
    double& w,
    double& y,
    double& h,
    double& fw,
    double& fh)
  {
    double r = 0.;
    r        = .25;
    gfx_set (r, g, b, a) /*TODO: unknown call */;
    gfx_line (
      x - r * w + .5 * w,
      y - r * w + .5 * h + 2.,
      x + r * w + .5 * w,
      y + r * w + .5 * h + 2.) /*TODO: unknown call */;
    gfx_line (
      x - r * w + .5 * w,
      y + r * w + .5 * h + 2.,
      x + r * w + .5 * w,
      y - r * w + .5 * h + 2.) /*TODO: unknown call */;
    gfx_setfont (9., fontface, 12. * (1. + scaling)) /*TODO: unknown call */;
    r = .29;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "L" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x - r * w + .5 * w - .5 * fw;
    gfx_y = y - r * w + .5 * h + 2. - .5 * fh;
    gfx_printf (
      0.; /* jsfx2cpp strings unsupported: was: "L" */) /*TODO: unknown call */;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "R" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x + r * w + .5 * w - .5 * fw;
    gfx_y = y - r * w + .5 * h + 2. - .5 * fh;
    gfx_printf (
      0.; /* jsfx2cpp strings unsupported: was: "R" */) /*TODO: unknown call */;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "R" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x + 10.;
    gfx_y = y + .5 * h - .5 * fh;
    gfx_printf (
      0.;
      /* jsfx2cpp strings unsupported: was: "+S" */) /*TODO: unknown call */;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "R" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x + w - 10. - fw;
    gfx_y = y + .5 * h - .5 * fh;
    gfx_printf (
      0.;
      /* jsfx2cpp strings unsupported: was: "-S" */) /*TODO: unknown call */;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "M" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x + .5 * w - .5 * fw;
    gfx_y = y + 5.;
    return gfx_printf (
      0.; /* jsfx2cpp strings unsupported: was: "M" */) /*TODO: unknown call */;
  }
#endif
//----------------------------------------------------------------------------
#if 0
  double init$feedgonio (
    double  l1,
    double  r1,
    double& recptr,
    double& rend,
    double& rstart)
  {
    mempool[(size_t) (recptr)]                 = l1;
    mempool[(size_t) ((recptr = recptr + 1.))] = r1;
    mempool[(size_t) ((recptr = recptr + 1.))] = r1;
    recptr                                     = [&] {
      if ((recptr + 1.) >= rend) {
        return rstart;
      }
      else {
        return recptr + 1.;
      }
    }();
    return recptr;
  }
#endif
  //----------------------------------------------------------------------------
  double init$init_linearsvf (
    double  freq,
    double  res,
    double& k,
    double& a1,
    double& a2,
    double& a3)
  {
    double g = 0.;
    g        = std::tan (.5 * freq);
    k        = 2. - 2. * res;
    a1       = 1. / (1. + g * (g + k));
    a2       = g * a1;
    a3       = g * a2;
    return a3;
  }
  //----------------------------------------------------------------------------
  double init$initbuffer (
    double  scopebuffer_in,
    double  scopebuffermax_in,
    double& scopebuffer,
    double& scopebuffermax,
    double& scopeptr)
  {
    scopebuffer    = scopebuffer_in;
    scopebuffermax = scopebuffermax_in;
    scopeptr       = scopebuffer;
    return [&] {
      if (scopeptr < scopebuffer) {
        (scopeptr = scopebuffer);
        return scopeptr;
      }
      else {
        return [&] {
          if ((scopeptr > scopebuffermax)) {
            return scopeptr;
          }
          else {
            return 0. / 0.;
          }
        }();
      }
    }();
  }
//----------------------------------------------------------------------------
#if 0
  double init$initgoniometer (
    double  in,
    double  samplecount,
    double& samples,
    double& rstart,
    double& rend,
    double& recptr,
    double& thisui)
  {
    samples = samplecount;
    rstart  = in;
    rend    = in + samplecount - 1.;
    recptr  = in;
    newui += 1.;
    thisui = newui;
    return thisui;
  }
#endif
  //----------------------------------------------------------------------------
  double init$readbuffer (
    double& readptr,
    double& scopebuffermax,
    double& scopebuffer,
    double& frac)
  {
    double c1 = 0.;
    double c2 = 0.;
    c1        = mempool[(size_t) (readptr)];
    readptr += 1.;
    if (readptr > scopebuffermax) {
      readptr = scopebuffer;
    }
    c2 = mempool[(size_t) (readptr)];
    return c2 * (1.0 - frac) + c1 * frac;
  }
  //----------------------------------------------------------------------------
  double init$setoffset (
    double  offset,
    double& readptr,
    double& scopeptr,
    double& frac,
    double& scopebuffer,
    double& scopebuffermax)
  {
    readptr = scopeptr;
    frac    = offset - std::floor (offset);
    readptr -= std::floor (offset);
    return [&] {
      if (readptr < scopebuffer) {
        readptr += (scopebuffermax - scopebuffer + 1.);
        return readptr;
      }
      else {
        return 0. / 0.;
      }
    }();
  }
//----------------------------------------------------------------------------
#if 0
  double init$setwindow (
    double  _x,
    double  _y,
    double  _w,
    double  _h,
    double& x,
    double& y,
    double& w,
    double& h)
  {
    x = _x;
    y = _y;
    w = _w;
    h = _h;
    return h;
  }
#endif
  //----------------------------------------------------------------------------
  double init$sliderupdate()
  {
    chdelay = get_slider_delay() * jsfx_specialvar_get_srate() / 2000.;
    lowcutfreq
      = std::exp ((1. - get_slider_lowcut()) * std::log (20. / 22050.));
    init$init_linearsvf (lowcutfreq, 0., hp$k, hp$a1, hp$a2, hp$a3);
    hp2$k  = hp$k;
    hp2$a1 = hp$a1;
    hp2$a2 = hp$a2;
    hp2$a3 = hp$a3;
#if 0
    displayfreq = .5 * jsfx_specialvar_get_srate() * lowcutfreq;
    return displayfreq;
#else
    return 0;
#endif
  }
  //----------------------------------------------------------------------------
  double init$tanh (double x) { return 2. / (std::exp (x * -2.) + 1.) - 1.; }
  //----------------------------------------------------------------------------
  double init$updatebuffer (
    double  m,
    double& scopeptr,
    double& scopebuffermax,
    double& scopebuffer)
  {
    mempool[(size_t) (scopeptr)] = m;
    scopeptr += 1.;
    if (scopeptr > scopebuffermax) {
      scopeptr = scopebuffer;
    }
    return m;
  }
  // functions for section "block"

  //----------------------------------------------------------------------------
  double block$eval_linearsvf_hp (
    double  v0,
    double& ic2eq,
    double& a1,
    double& ic1eq,
    double& a2,
    double& a3,
    double& k,
    double& ic4eq,
    double& ic3eq,
    double& ic6eq,
    double& ic5eq,
    double& ic8eq,
    double& ic7eq,
    double& ic10eq,
    double& ic9eq,
    double& ic12eq,
    double& ic11eq)
  {
    double v2 = 0.;
    double v3 = 0.;
    double v1 = 0.;
    v3        = v0 - ic2eq;
    v1        = a1 * ic1eq + a2 * v3;
    v2        = ic2eq + a2 * ic1eq + a3 * v3;
    ic1eq     = 2. * v1 - ic1eq;
    ic2eq     = 2. * v2 - ic2eq;
    v0        = v0 - k * v1 - v2;
    v3        = v0 - ic4eq;
    v1        = a1 * ic3eq + a2 * v3;
    v2        = ic4eq + a2 * ic3eq + a3 * v3;
    ic3eq     = 2. * v1 - ic3eq;
    ic4eq     = 2. * v2 - ic4eq;
    v0        = v0 - k * v1 - v2;
    v3        = v0 - ic6eq;
    v1        = a1 * ic5eq + a2 * v3;
    v2        = ic6eq + a2 * ic5eq + a3 * v3;
    ic5eq     = 2. * v1 - ic5eq;
    ic6eq     = 2. * v2 - ic6eq;
    v0        = v0 - k * v1 - v2;
    v3        = v0 - ic8eq;
    v1        = a1 * ic7eq + a2 * v3;
    v2        = ic8eq + a2 * ic7eq + a3 * v3;
    ic7eq     = 2. * v1 - ic7eq;
    ic8eq     = 2. * v2 - ic8eq;
    v0        = v0 - k * v1 - v2;
    v3        = v0 - ic10eq;
    v1        = a1 * ic9eq + a2 * v3;
    v2        = ic10eq + a2 * ic9eq + a3 * v3;
    ic9eq     = 2. * v1 - ic9eq;
    ic10eq    = 2. * v2 - ic10eq;
    v0        = v0 - k * v1 - v2;
    v3        = v0 - ic12eq;
    v1        = a1 * ic11eq + a2 * v3;
    v2        = ic12eq + a2 * ic11eq + a3 * v3;
    ic11eq    = 2. * v1 - ic11eq;
    ic12eq    = 2. * v2 - ic12eq;
    return v0 - k * v1 - v2;
  }
}; /* stereo_bub3 */
}} // namespace artv::saike
