// Ported from https://github.com/Sonic-Anomaly/Sonic-Anomaly-JSFX.git
// commit sha: 663735021e59d788e8ebacc4b9678569378deec9

// Mods/notes:
// - Required less memory on the delay buffers

#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>
// slider2:0<0,1,1{Off,On}>-SC Filter
// slider3:0<0,8,1{Off,20Hz,40Hz,80Hz,100Hz,150Hz,200Hz,350Hz,500Hz}>-LF Filter
// slider4:0<0,11,1{Zero,250us,500us,1ms,5ms,10ms,20ms,50ms,100ms,200ms,500ms,1sec}>-Attack
// Mode slider5:0<-30,0,1>-Push Down [dB] slider6:0<0,30,1>-Pull Up [dB]
// slider7:-5<-10,0,1>-Recovery
// slider8:0<0,10,1>-Emphasis
// slider9:0<-20,20,1>-Output [dB]
// slider10:10<1,10,0.1>-[Display Speed]

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/jsfx_engine/jsfx_engine.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace sonic_anomaly {

struct vola2 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::dynamics;
  //----------------------------------------------------------------------------
private:
  // definitions for environment function calls
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }

  std::vector<float> heapmem;

  inline float& heap (std::size_t value) { return heapmem[value]; }

  void heap_reset (std::size_t s)
  {
    heapmem.resize (s);
    std::memset (heapmem.data(), 0, heapmem.size() * sizeof heapmem[0]);
  }

  void jsfx_memset (size_t idx, int val, size_t size)
  {
    std::memset (&heapmem[idx], val, size);
  }

  static double jsfx_rand (double maxv = 1.)
  {
    return jsfx_engine::rand (maxv);
  }
//----------------------------------------------------------------------------
// stubs for JSFX special variables
#if 0
  double jsfx_specialvar_get_play_state() { return 1.; }
#endif

  void jsfx_specialvar_set_pdc_bot_ch (double val)
  {
    /* TODO: stub for setting JSFX var "pdc_bot_ch" */
    // latency compensation ignored for now...
  }

  void jsfx_specialvar_set_pdc_delay (double val)
  {
    plugcontext->set_delay_compensation ((uint) val);
  }

  void jsfx_specialvar_set_pdc_top_ch (double val)
  {
    /* TODO: stub for setting JSFX var "pdc_top_ch" */
    // latency compensation ignored for now...
  }
  //----------------------------------------------------------------------------
  // stubs for sliders
public:
#if 1
  double get_slider_slider10()
  {
    // TODO: stub, add code for getting "slider10"
    // Range: min:1.0, max:10.0, default: 10.0, step: 0.1
    // Original line: slider10:10<1,10,0.1>-[Display Speed]
    return 1.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider10()
  {
    // Range: min:1.0, max:10.0, default: 10.0, step: 0.1
    return slider10p;
  }

  float slider10p = 10.0;

  struct slider10_tag {};

  void set (slider10_tag, float v)
  {
    if (v == slider10p) {
      return;
    }
    slider10p = v;
    slider();
  }

  static constexpr auto get_parameter (slider10_tag)
  {
    // Original slider line: slider10:10<1,10,0.1>-[Display Speed]
    return float_param ("", 1.0, 10.0, 10.0, 0.1);
  }
#endif
#if 0
  double get_slider_slider2()
  {
    // TODO: stub, add code for getting "slider2"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    // Original line: slider2:0<0,1,1{Off,On}>-SC Filter
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider2()
  {
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    return slider2p;
  }

  float slider2p = 0.0;

  struct sc_filter_tag {};

  void set (sc_filter_tag, float v)
  {
    if (v == slider2p) {
      return;
    }
    slider2p = v;
    slider();
  }

  static constexpr auto get_parameter (sc_filter_tag)
  {
    // Original slider line: slider2:0<0,1,1{Off,On}>-SC Filter
    return float_param ("", 0.0, 1.0, 0.0, 1.0);
  }
#endif
#if 0
  double get_slider_slider3()
  {
    // TODO: stub, add code for getting "slider3"
    // Range: min:0.0, max:8.0, default: 0.0, step: 1.0
    // Original line:
    // slider3:0<0,8,1{Off,20Hz,40Hz,80Hz,100Hz,150Hz,200Hz,350Hz,500Hz}>-LF
    // Filter
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider3()
  {
    // Range: min:0.0, max:8.0, default: 0.0, step: 1.0
    return slider3p;
  }

  float slider3p = 0.0;

  struct lf_tag {};

  static constexpr auto get_parameter (lf_tag)
  {
    // Original slider line:
    // slider3:0<0,8,1{Off,20Hz,40Hz,80Hz,100Hz,150Hz,200Hz,350Hz,500Hz}>-LF
    // Filter
    return choice_param (
      0,
      make_cstr_array (
        "Off",
        "20Hz",
        "40Hz",
        "80Hz",
        "100Hz",
        "150Hz",
        "200Hz",
        "350Hz",
        "500Hz"),
      20);
  }

  void set (lf_tag, float v)
  {
    if (v == slider3p || v >= get_parameter (lf_tag {}).choices.size()) {
      return;
    }
    slider3p = v;
    slider();
  }

#endif
#if 0
  double get_slider_slider4()
  {
    // TODO: stub, add code for getting "slider4"
    // Range: min:0.0, max:11.0, default: 0.0, step: 1.0
    // Original line:
    // slider4:0<0,11,1{Zero,250us,500us,1ms,5ms,10ms,20ms,50ms,100ms,200ms,500ms,1sec}>-Attack
    // Mode
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider4()
  {
    // Range: min:0.0, max:11.0, default: 0.0, step: 1.0
    return slider4p;
  }

  float slider4p = 0.0;

  struct attack_tag {};

  static constexpr auto get_parameter (attack_tag)
  {
    // Original slider line:
    // slider4:0<0,11,1{Zero,250us,500us,1ms,5ms,10ms,20ms,50ms,100ms,200ms,500ms,1sec}>-Attack
    // Mode
    return choice_param (
      0,
      make_cstr_array (
        "Zero",
        "250us",
        "500us",
        "1ms",
        "5ms",
        "10ms",
        "20ms",
        "50ms",
        "100ms",
        "200ms",
        "500ms",
        "1sec"),
      20);
  }

  void set (attack_tag, float v)
  {
    if (v == slider4p || v >= get_parameter (attack_tag {}).choices.size()) {
      return;
    }
    slider4p = v;
    slider();
  }

#endif
#if 0
  double get_slider_slider5()
  {
    // TODO: stub, add code for getting "slider5"
    // Range: min:-30.0, max:0.0, default: 0.0, step: 1.0
    // Original line: slider5:0<-30,0,1>-Push Down [dB]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider5()
  {
    // Range: min:-30.0, max:0.0, default: 0.0, step: 1.0
    return slider5p;
  }

  float slider5p = 0.0;

  struct push_down_tag {};

  void set (push_down_tag, float v)
  {
    if (v == slider5p) {
      return;
    }
    slider5p = v;
    slider();
  }

  static constexpr auto get_parameter (push_down_tag)
  {
    // Original slider line: slider5:0<-30,0,1>-Push Down [dB]
    return float_param ("dB", -30.0, 0.0, 0.0, 1.0);
  }
#endif
#if 0
  double get_slider_slider6()
  {
    // TODO: stub, add code for getting "slider6"
    // Range: min:0.0, max:30.0, default: 0.0, step: 1.0
    // Original line: slider6:0<0,30,1>-Pull Up [dB]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider6()
  {
    // Range: min:0.0, max:30.0, default: 0.0, step: 1.0
    return slider6p;
  }

  float slider6p = 0.0;

  struct pull_up_tag {};

  void set (pull_up_tag, float v)
  {
    if (v == slider6p) {
      return;
    }
    slider6p = v;
    slider();
  }

  static constexpr auto get_parameter (pull_up_tag)
  {
    // Original slider line: slider6:0<0,30,1>-Pull Up [dB]
    return float_param ("dB", 0.0, 30.0, 0.0, 1.0);
  }
#endif
#if 0
  double get_slider_slider7()
  {
    // TODO: stub, add code for getting "slider7"
    // Range: min:-10.0, max:0.0, default: -5.0, step: 1.0
    // Original line: slider7:-5<-10,0,1>-Recovery
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider7()
  {
    // Range: min:-10.0, max:0.0, default: -5.0, step: 1.0
    return slider7p;
  }

  float slider7p = -5.0;

  struct recovery_tag {};

  void set (recovery_tag, float v)
  {
    if (v == slider7p) {
      return;
    }
    slider7p = v;
    slider();
  }

  static constexpr auto get_parameter (recovery_tag)
  {
    // Original slider line: slider7:-5<-10,0,1>-Recovery
    return float_param ("", -10.0, 0.0, -5.0, 1.0);
  }
#endif
#if 0
  double get_slider_slider8()
  {
    // TODO: stub, add code for getting "slider8"
    // Range: min:0.0, max:10.0, default: 0.0, step: 1.0
    // Original line: slider8:0<0,10,1>-Emphasis
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider8()
  {
    // Range: min:0.0, max:10.0, default: 0.0, step: 1.0
    return slider8p;
  }

  float slider8p = 0.0;

  struct emphasis_tag {};

  void set (emphasis_tag, float v)
  {
    if (v == slider8p) {
      return;
    }
    slider8p = v;
    slider();
  }

  static constexpr auto get_parameter (emphasis_tag)
  {
    // Original slider line: slider8:0<0,10,1>-Emphasis
    return float_param ("", 0.0, 10.0, 0.0, 1.0);
  }
#endif
#if 1
  double get_slider_slider9()
  {
    // TODO: stub, add code for getting "slider9"
    // Range: min:-20.0, max:20.0, default: 0.0, step: 1.0
    // Original line: slider9:0<-20,20,1>-Output [dB]
    return 0.;
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  double get_slider_slider9()
  {
    // Range: min:-20.0, max:20.0, default: 0.0, step: 1.0
    return slider9p;
  }

  float slider9p = 0.0;

  struct slider9_tag {};

  void set (slider9_tag, float v)
  {
    if (v == slider9p) {
      return;
    }
    slider9p = v;
    slider();
  }

  static constexpr auto get_parameter (slider9_tag)
  {
    // Original slider line: slider9:0<-20,20,1>-Output [dB]
    return float_param ("", -20.0, 20.0, 0.0, 1.0);
  }
#endif
  using parameters = mp_list<
    sc_filter_tag,
    lf_tag,
    attack_tag,
    push_down_tag,
    pull_up_tag,
    recovery_tag,
    emphasis_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double _sliderdirty;
  double att;
#if 0
  double buffsize;
  double buildstr;
#endif
  double d0$len;
  double d0$sloop;
  double d0$splay;
  double d1$len;
  double d1$sloop;
  double d1$splay;
  double emphasis;
  double emul;
  double flt0$a1;
  double flt0$a2;
  double flt0$a3;
  double flt0$b1;
  double flt0$b2;
  double flt0$c;
  double flt1$a1;
  double flt1$a2;
  double flt1$a3;
  double flt1$b1;
  double flt1$b2;
  double flt1$c;
  double gfxmem_start;
  double gr_pos_x;
  double gr_pos_y;
  double i;
  double inertia;
#if 0
  double loopsize;
#endif
  double maxbcount;
  double mgain;
  double mlim;
#if 0
  double mtr_x;
  double mtr_y;
  double mtrspeed;
#endif
  double output;
  double p_tmp;
  double p_tmp2;
  double prefilt;
  double recovery;
  double rel;
  double rmm;
  double rmss;
  double rmt;
  double s;
  double scf;
  double scf$a1;
  double scf$a2;
  double scf$a3;
  double scf$b1;
  double scf$b2;
  double scf$c;
  double scfilt;
  double sniffer$ismono;
  double ss;
  double ssf;
  double ssr;
#if 0
  double step;
  double strmem;
#endif
  double t20;
  double t200;
#if 0
  double x_pos;
  double y_mem;
  double y_mem2;
  double y_pos;
#endif
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    _sliderdirty = 0;
    att          = 0;
#if 0
    buffsize       = 0;
    buildstr       = 0;
#endif
    d0$len   = 0;
    d0$sloop = 0;
    d0$splay = 0;
    d1$len   = 0;
    d1$sloop = 0;
    d1$splay = 0;
    emphasis = 0;
    emul     = 0;
    flt0$a1  = 0;
    flt0$a2  = 0;
    flt0$a3  = 0;
    flt0$b1  = 0;
    flt0$b2  = 0;
    flt0$c   = 0;
    flt1$a1  = 0;
    flt1$a2  = 0;
    flt1$a3  = 0;
    flt1$b1  = 0;
    flt1$b2  = 0;
    flt1$c   = 0;
#if 0
    gfxmem_start   = 0;
    gr_pos_x       = 0;
    gr_pos_y       = 0;
#endif
    i       = 0;
    inertia = 0;
#if 0
    loopsize  = 0;
    maxbcount = 0;
#endif
    mgain = 0;
    mlim  = 0;
#if 0
    mtr_x          = 0;
    mtr_y          = 0;
    mtrspeed       = 0;
#endif
    output         = 0;
    p_tmp          = 0;
    p_tmp2         = 0;
    prefilt        = 0;
    recovery       = 0;
    rel            = 0;
    rmm            = 0;
    rmss           = 0;
    rmt            = 0;
    s              = 0;
    scf            = 0;
    scf$a1         = 0;
    scf$a2         = 0;
    scf$a3         = 0;
    scf$b1         = 0;
    scf$b2         = 0;
    scf$c          = 0;
    scfilt         = 0;
    sniffer$ismono = 0;
    ss             = 0;
    ssf            = 0;
    ssr            = 0;
#if 0
    step           = 0;
    strmem         = 0;
#endif
    t20  = 0;
    t200 = 0;
#if 0
    x_pos          = 0;
    y_mem          = 0;
    y_mem2         = 0;
    y_pos          = 0;
#endif
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double blocksniffer$mspl;
  double d0$sindex;
  double d1$sindex;
  double e0$env;
  double e0$tmp;
  double e1$env;
  double e1$tmp;
  double e2$env;
  double e2$tmp;
  double e3$env;
  double e3$tmp;
  double flt0$mx1;
  double flt0$mx2;
  double flt0$my1;
  double flt0$my2;
  double flt0$output;
  double flt1$mx1;
  double flt1$mx2;
  double flt1$my1;
  double flt1$my2;
  double flt1$output;
  double genv;
  double init;
  double scf$mx1;
  double scf$mx2;
  double scf$my1;
  double scf$my2;
  double scf$output;
  double snl;
  double snr;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    blocksniffer$mspl = 0;
    d0$sindex         = 0;
    d1$sindex         = 0;
    e0$env            = 0;
    e0$tmp            = 0;
    e1$env            = 0;
    e1$tmp            = 0;
    e2$env            = 0;
    e2$tmp            = 0;
    e3$env            = 0;
    e3$tmp            = 0;
    flt0$mx1          = 0;
    flt0$mx2          = 0;
    flt0$my1          = 0;
    flt0$my2          = 0;
    flt0$output       = 0;
    flt1$mx1          = 0;
    flt1$mx2          = 0;
    flt1$my1          = 0;
    flt1$my2          = 0;
    flt1$output       = 0;
    genv              = 0;
    init              = 0;
    scf$mx1           = 0;
    scf$mx2           = 0;
    scf$my1           = 0;
    scf$my2           = 0;
    scf$output        = 0;
    snl               = 0;
    snr               = 0;
  }
  //----------------------------------------------------------------------------
  plugin_context* plugcontext = nullptr;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_block_variables();
#if 0
    buildstr    = 0.; /* jsfx2cpp strings unsupported: was: "Build 170206" */
#endif

    p_tmp2    = 0.;
    p_tmp     = p_tmp2;
    t20       = std::exp (-1. / (plugcontext->get_sample_rate() * 0.02));
    t200      = std::exp (-1. / (plugcontext->get_sample_rate() * 0.2));
    rmm       = (44100. / plugcontext->get_sample_rate()) * 0.0003;
    maxbcount = 32.;
    init$filterhp_init (
      2000., 5., scf$c, scf$a1, scf$a2, scf$a3, scf$b1, scf$b2);
#if 1
    heap_reset (128);
#endif
    init$delay_init (64., 0., d0$len, d0$splay, d0$sloop);
    init$delay_init (64., 1., d1$len, d1$splay, d1$sloop);
    jsfx_specialvar_set_pdc_delay (64.);
    jsfx_specialvar_set_pdc_bot_ch (0.);
    jsfx_specialvar_set_pdc_top_ch (2.);
#if 0
    mtr_x        = 17.;
    mtr_y        = 290.;
    y_mem2       = 360.;
    y_mem        = y_mem2;
#endif
    rmss = 10000.;
#if 0
    x_pos        = mtr_x + 90.;
    y_pos        = mtr_y + 190.;
    gr_pos_x     = mtr_x + 10.;
    gr_pos_y     = mtr_y + 10.;
    gfxmem_start = plugcontext->get_sample_rate() * 2. + 1.;
    loopsize = 360.;
#endif
#if 0
    strmem       = 600000.;
    heap (strmem) = 0.; /* jsfx2cpp strings unsupported: was: "Off" */
    ;
    heap (strmem + 1.) = 0.; /* jsfx2cpp strings unsupported: was: "20 Hz" */
    ;
    heap (strmem + 2.) = 0.; /* jsfx2cpp strings unsupported: was: "40 Hz" */
    ;
    heap (strmem + 3.) = 0.; /* jsfx2cpp strings unsupported: was: "80 Hz" */
    ;
    heap (strmem + 4.) = 0.; /* jsfx2cpp strings unsupported: was: "100 Hz" */
    ;
    heap (strmem + 5.) = 0.; /* jsfx2cpp strings unsupported: was: "150 Hz" */
    ;
    heap (strmem + 6.) = 0.; /* jsfx2cpp strings unsupported: was: "200 Hz" */
    ;
    heap (strmem + 7.) = 0.; /* jsfx2cpp strings unsupported: was: "350 Hz" */
    ;
    heap (strmem + 8.) = 0.; /* jsfx2cpp strings unsupported: was: "500 Hz" */
    ;
    strmem        = 600010.;
    heap (strmem) = 0.; /* jsfx2cpp strings unsupported: was: "Zero" */
    ;
    heap (strmem + 1.) = 0.; /* jsfx2cpp strings unsupported: was: "250us" */
    ;
    heap (strmem + 2.) = 0.; /* jsfx2cpp strings unsupported: was: "500 us" */
    ;
    heap (strmem + 3.) = 0.; /* jsfx2cpp strings unsupported: was: "1 ms" */
    ;
    heap (strmem + 4.) = 0.; /* jsfx2cpp strings unsupported: was: "5 ms" */
    ;
    heap (strmem + 5.) = 0.; /* jsfx2cpp strings unsupported: was: "10 ms" */
    ;
    heap (strmem + 6.) = 0.; /* jsfx2cpp strings unsupported: was: "20 ms" */
    ;
    heap (strmem + 7.) = 0.; /* jsfx2cpp strings unsupported: was: "50 ms" */
    ;
    heap (strmem + 8.) = 0.; /* jsfx2cpp strings unsupported: was: "100 ms" */
    ;
    heap (strmem + 9.) = 0.; /* jsfx2cpp strings unsupported: was: "200 ms" */
    ;
    heap (strmem + 10.) = 0.; /* jsfx2cpp strings unsupported: was: "500 ms" */
    ;
    heap (strmem + 11.) = 0.; /* jsfx2cpp strings unsupported: was: "1 s" */
    ;
#endif
    _sliderdirty = 1.;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double s10$value  = 0.;
    double s2$checked = 0.;
    double s3$value   = 0.;
    double s4$value   = 0.;
    double s5$value   = 0.;
    double s6$value   = 0.;
    double s7$value   = 0.;
    double s8$value   = 0.;
    double s9$value   = 0.;
    s2$checked        = get_slider_slider2();
    s3$value          = get_slider_slider3();
    s4$value          = get_slider_slider4();
    s5$value          = get_slider_slider5();
    s6$value          = get_slider_slider6();
    s7$value          = get_slider_slider7();
    s8$value          = get_slider_slider8();
    s9$value          = get_slider_slider9();
    s10$value         = get_slider_slider10();
    init$process_slider();
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    double env    = 0.;
    double gfxmem = 0.;
    double inl    = 0.;
    double inr    = 0.;
    double key    = 0.;
    double norm   = 0.;
    double outl   = 0.;
    double outr   = 0.;
    double rmod   = 0.;

#if 0
    if (!init) {
      init$resetdisplay();
      init = 1.;
    }

    init$blocksniffer (snl, snr, blocksniffer$mspl);
    norm = (jsfx_rand (1.) * 0.000000000000001);
    if (jsfx_specialvar_get_play_state() && genv) {
      gfxmem            = gfxmem_start;
      heap (gfxmem + i) = 200.
        - std::ceil (100.
                     - std::log (
                         1. / std::max (std::min (genv * output, 10.), 0.1))
                       * 43.3);
      i += 1.;
      if (i >= buffsize / samples) {
        i = 0.;
      }
      else {
        i;
      }
    };
#endif
    for (int $$i = 0, $$end = std::max (0, (int) (samples)); $$i < $$end;
         ++$$i) {
      T& spl0 = chnls[0][$$i];
      T& spl1 = chnls[1][$$i];
      snl     = spl0 + norm;
      inl     = snl;
      snr     = spl1 + norm;
      inr     = snr;
      if (prefilt) {
        inl = init$filterhp (
          inl,
          flt0$output,
          flt0$a1,
          flt0$a2,
          flt0$mx1,
          flt0$a3,
          flt0$mx2,
          flt0$b1,
          flt0$my1,
          flt0$b2,
          flt0$my2);
        (inr = inl);
        if (!sniffer$ismono) {
          (inr = init$filterhp (
             inr,
             flt1$output,
             flt1$a1,
             flt1$a2,
             flt1$mx1,
             flt1$a3,
             flt1$mx2,
             flt1$b1,
             flt1$my1,
             flt1$b2,
             flt1$my2));
        }
        else {
          inr;
        }
      }
      key = [&] {
        if (std::abs (inr) > std::abs (inl)) {
          return inr;
        }
        else {
          return inl;
        }
      }();
      key =
        [&] {
          if (emphasis) {
            return init$filterhp (
                     key,
                     scf$output,
                     scf$a1,
                     scf$a2,
                     scf$mx1,
                     scf$a3,
                     scf$mx2,
                     scf$b1,
                     scf$my1,
                     scf$b2,
                     scf$my2)
              * emphasis;
          }
          else {
#if 0
            return 0. / 0.;
#else
            // this is a jsfx2cpp bug. To be honest what happens on this line
            // is completely obscure and insane. Luckily jsfx2cpp generated a
            // Nan. Original expression:
            //
            // key = (emphasis ? scf.filterHP(key) * emphasis) + key;
            //
            // Any sane language would require you to explicitly write the
            // else branch explicitly...
            return 0.;
#endif
          }
        }()
        + key;
      if (scfilt) {
        (key = init$filterhp (
           key,
           flt0$output,
           flt0$a1,
           flt0$a2,
           flt0$mx1,
           flt0$a3,
           flt0$mx2,
           flt0$b1,
           flt0$my1,
           flt0$b2,
           flt0$my2));
      }
      rmod = init$follower (std::abs (key), 0., t200 - rmt, e3$tmp, e3$env);
      rmod = rel / (rmm * rmod + 1.);
      env  = init$follower (std::abs (key), att, rmod, e0$tmp, e0$env);
      env  = init$follower (env, 0., t20, e1$tmp, e1$env);
      env  = init$follower (env, 0., t20, e2$tmp, e2$env);
      env  = mlim / (env + mgain);
      outl = env * init$delay (inl, d0$sloop, d0$sindex, d0$len, d0$splay);
      outr = env * init$delay (inr, d1$sloop, d1$sindex, d1$len, d1$splay);
      genv = env;
      spl0 = outl * output;
      spl1 = outr * output;
      ;
    }
  }
  // functions for section "init"
private:
//----------------------------------------------------------------------------
#if 0
  double init$blocksniffer (double inl, double inr, double& mspl)
  {
    return [&] {
      if ((mspl < maxbcount)) {
        return [&] {
          if (eel2_eq (inl, inr)) {
            (mspl += 1.);
            return mspl;
          }
          else {
            mspl           = 0.;
            sniffer$ismono = 0.;
            return sniffer$ismono;
          }
        }();
      }
      else {
        sniffer$ismono = 1.;
        mspl           = 0.;
        return mspl;
      }
    }();
  }
#endif
  //----------------------------------------------------------------------------
  double init$delay (
    double  input,
    double& sloop,
    double& sindex,
    double& len,
    double& splay)
  {
    heap (sloop + sindex) = input;
    sindex += 1.;
#if 0
    if (sindex > len) {
#else
    if (sindex >= len) {
#endif
    sindex = 0.;
  }
  return heap (splay + sindex);
}
//----------------------------------------------------------------------------
double
init$delay_init (
  double  samples,
  double  index,
  double& len,
  double& splay,
  double& sloop)
{
#if 0
    len = [&] {
      if (samples > plugcontext->get_sample_rate()) {
        return (double) plugcontext->get_sample_rate();
      }
      else {
        return samples;
      }
    }();
    splay = plugcontext->get_sample_rate() * index;
    sloop = splay;
    return sloop;
#else
    // this function is called with 64 samples delay, we are going to play at
    // least at 44100
    len   = samples;
    splay = len * index;
    sloop = splay;
    return sloop;
#endif
}
//----------------------------------------------------------------------------
double init$filterhp (
  double  input,
  double& output,
  double& a1,
  double& a2,
  double& mx1,
  double& a3,
  double& mx2,
  double& b1,
  double& my1,
  double& b2,
  double& my2)
{
  output = a1 * input + a2 * mx1 + a3 * mx2 - b1 * my1 - b2 * my2;
  mx2    = mx1;
  mx1    = input;
  my2    = my1;
  my1    = output;
  return my1;
}
//----------------------------------------------------------------------------
double init$filterhp_init (
  double  f,
  double  r,
  double& c,
  double& a1,
  double& a2,
  double& a3,
  double& b1,
  double& b2)
{
  c  = std::tan (3.141592653589793 * f / plugcontext->get_sample_rate());
  a1 = 1.0 / (1.0 + r * c + c * c);
  a2 = -2. * a1;
  a3 = a1;
  b1 = 2.0 * (c * c - 1.0) * a1;
  b2 = (1.0 - r * c + c * c) * a1;
  return b2;
}
//----------------------------------------------------------------------------
double init$follower (
  double  input,
  double  att,
  double  rel,
  double& tmp,
  double& env)
{
  tmp = input;
  (env = rel * (env - tmp) + tmp);
  return [&] {
    if ((tmp > env)) {
      (env = att * (env - tmp) + tmp);
      return env;
    }
    else {
      return env;
    }
  }();
}
//----------------------------------------------------------------------------
double init$process_slider()
{
  p_tmp2 = 0.;
  p_tmp  = p_tmp2;
  mlim   = eel2_pow (10., (get_slider_slider5() / 20.));
  mgain
    = eel2_pow (10., ((-get_slider_slider6() + get_slider_slider5()) / 20.));
  if (eel2_eq (get_slider_slider7(), -10.)) {
    recovery = 100.;
  }
  if (eel2_eq (get_slider_slider7(), -9.)) {
    recovery = 10.;
  }
  if (eel2_eq (get_slider_slider7(), -8.)) {
    recovery = 5.;
  }
  if (eel2_eq (get_slider_slider7(), -7.)) {
    recovery = 1.;
  }
  if (eel2_eq (get_slider_slider7(), -6.)) {
    recovery = 0.5;
  }
  if (eel2_eq (get_slider_slider7(), -5.)) {
    recovery = 0.25;
  }
  if (eel2_eq (get_slider_slider7(), -4.)) {
    recovery = 0.15;
  }
  if (eel2_eq (get_slider_slider7(), -3.)) {
    recovery = 0.1;
  }
  if (eel2_eq (get_slider_slider7(), -2.)) {
    recovery = 0.08;
  }
  if (eel2_eq (get_slider_slider7(), -1.)) {
    recovery = 0.06;
  }
  if (eel2_eq (get_slider_slider7(), 0.)) {
    recovery = 0.03;
  }
  rel = std::exp (-1. / (plugcontext->get_sample_rate() * recovery));
  rmt = [&] {
    if (eel2_eq (recovery, 100.)) {
      return 0.00003;
    }
    else {
      return 0.;
    }
  }();
  output = eel2_pow (10., (get_slider_slider9() / 20.));
  ss     = get_slider_slider3();
  if (eel2_eq (ss, 1.)) {
    ssf = 20.;
    ssr = 1.42;
  }
  if (eel2_eq (ss, 2.)) {
    ssf = 40.;
    ssr = 1.42;
  }
  if (eel2_eq (ss, 3.)) {
    ssf = 60.;
    ssr = 1.75;
  }
  if (eel2_eq (ss, 4.)) {
    ssf = 60.;
    ssr = 2.;
  }
  if (eel2_eq (ss, 5.)) {
    ssf = 60.;
    ssr = 3.;
  }
  if (eel2_eq (ss, 6.)) {
    ssf = 60.;
    ssr = 4.;
  }
  if (eel2_eq (ss, 7.)) {
    ssf = 60.;
    ssr = 6.;
  }
  if (eel2_eq (ss, 8.)) {
    ssf = 60.;
    ssr = 8.;
  }
  init$filterhp_init (
    ssf, ssr, flt0$c, flt0$a1, flt0$a2, flt0$a3, flt0$b1, flt0$b2);
  init$filterhp_init (
    ssf, ssr, flt1$c, flt1$a1, flt1$a2, flt1$a3, flt1$b1, flt1$b2);
  if (eel2_eq (get_slider_slider4(), 0.)) {
    inertia = 0.;
  }
  if (eel2_eq (get_slider_slider4(), 1.)) {
    inertia = 0.00025;
  }
  if (eel2_eq (get_slider_slider4(), 2.)) {
    inertia = 0.0005;
  }
  if (eel2_eq (get_slider_slider4(), 3.)) {
    inertia = 0.001;
  }
  if (eel2_eq (get_slider_slider4(), 4.)) {
    inertia = 0.005;
  }
  if (eel2_eq (get_slider_slider4(), 5.)) {
    inertia = 0.01;
  }
  if (eel2_eq (get_slider_slider4(), 6.)) {
    inertia = 0.02;
  }
  if (eel2_eq (get_slider_slider4(), 7.)) {
    inertia = 0.05;
  }
  if (eel2_eq (get_slider_slider4(), 8.)) {
    inertia = 0.1;
  }
  if (eel2_eq (get_slider_slider4(), 9.)) {
    inertia = 0.2;
  }
  if (eel2_eq (get_slider_slider4(), 10.)) {
    inertia = 0.5;
  }
  if (eel2_eq (get_slider_slider4(), 11.)) {
    inertia = 1.;
  }
  att      = std::exp (-1. / (plugcontext->get_sample_rate() * inertia));
  emul     = eel2_pow ((get_slider_slider8() * 0.1), 0.65) * 3.16;
  emphasis = std::min ((inertia * 2000.) + emul, 10.);
#if 0
    mtrspeed = get_slider_slider10();
    buffsize = plugcontext->get_sample_rate() * mtrspeed;
    step     = buffsize / jsfx_specialvar_get_samplesblock() / loopsize;
#endif
  scf    = get_slider_slider2();
  scfilt = [&] {
    if (ss && scf) {
      return 1.;
    }
    else {
      return 0.;
    }
  }();
  prefilt = [&] {
    if (ss && !scf) {
      return 1.;
    }
    else {
      return 0.;
    }
  }();
  return prefilt;
}
//----------------------------------------------------------------------------
#if 0
  void init$resetdisplay()
  {
    i = 0.;
    s = 10000.;
    jsfx_memset (gfxmem_start, 100., buffsize);
  }
#endif
}; }
} // namespace artv::sonic_anomaly
