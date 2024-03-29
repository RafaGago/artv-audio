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

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace saike {

// slider1:Delay=8<1,25,1>-Delay [ms]
// slider2:Strength=0.3<0,1,.0001>-Strength [-]
// slider3:LowCut=0.4<0,1,.0001>-Crossover [log(w)]
// slider4:SideLevel=1<0,2,.0001>-Blend old side level [-]
// slider5:CheckMono=0<0,1,1{No, Yes}>-Check Mono?
// slider6:SideHP=0<0,1,1{No,Yes}>-Pass side through HPF?
struct stereo_bub2 {
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
  // 25ms buffer (0.025s). (1/40s); 192000 / 40 = 4800 we round to 5000 samples
  // static constexpr int mempool_size          = 5000;
  // float                mempool[mempool_size] = {};
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

  double jsfx_specialvar_get_srate()
  {
    return sample_rate;
  }

public:
  //----------------------------------------------------------------------------
  // stubs for sliders

  double get_slider_CheckMono()
  {
    // TODO: stub, add code for getting "CheckMono"
    return 0.;
  }

  //----------------------------------------------------------------------------
  struct delay_tag {};
  void set (delay_tag, float v)
  {
    // from init$sliderUpdate()
    chDelay = v * jsfx_specialvar_get_srate() / 2000.;
  }
  static constexpr auto get_parameter (delay_tag)
  {
    return float_param ("ms", 1.f, 25.f, 8.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct low_cut_tag {};
  void set (low_cut_tag, float v)
  {
    v = midi_note_to_hz (v);
    v /= 22000.f;
    if (v == low_cutp) {
      // don't reset the filter...
      return;
    }
    low_cutp = v;
    // from init$sliderUpdate()
    lowCutFreq
      = std::exp ((1. - get_slider_LowCut()) * std::log (20. / 22050.));
    init$init_linearSVF (lowCutFreq, 0., hp$k, hp$a1, hp$a2, hp$a3);
    hp2$k  = hp$k;
    hp2$a1 = hp$a1;
    hp2$a2 = hp$a2;
    hp2$a3 = hp$a3;
  }
  static constexpr auto get_parameter (low_cut_tag)
  {
    return frequency_parameter (19.0, 22000.0, 350.0);
  }

  float  low_cutp = -1.f;
  double get_slider_LowCut()
  {
    return low_cutp;
  }
  //----------------------------------------------------------------------------
  struct side_level_tag {};
  void set (side_level_tag, float v)
  {
    sideLevel = v;
  }
  static constexpr auto get_parameter (side_level_tag)
  {
    return float_param ("", 0.f, 2.f, 1.f, 0.01f);
  }

  double get_slider_SideHP()
  {
    return sideLevel;
  }
  //----------------------------------------------------------------------------
  struct strength_tag {};
  void set (strength_tag, float v)
  {
    strengthp = v;
  }
  static constexpr auto get_parameter (strength_tag)
  {
    return float_param ("", 0.f, 1.f, 0.4f, 0.0001f);
  }

  float  strengthp = 0.4f;
  double get_slider_Strength()
  {
    // TODO: stub, add code for getting "Strength"
    return strengthp;
  }
  //----------------------------------------------------------------------------
  using parameters
    = mp_list<delay_tag, low_cut_tag, side_level_tag, strength_tag>;
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
private:
  double MAXBUFFERSIZE = {};
#if 1
  double MaxDelay {};
  double MaxLowCut {};
  double MaxSideLevel {};
  double MaxStrength {};
  double MinDelay {};
  double MinLowCut {};
  double MinSideLevel {};
  double MinStrength {};
  double buffer$scopebuffer {};
  double buffer$scopebuffermax {};
  double buffer$scopeptr {};
  double bufferDist {};
#endif
  double chDelay {};
  double delayBuf1 {};
  double delta_t {};
  double displayFreq {};
#if 0
  double fontface         = {};
  double gfx_a            = {};
  double gfx_x            = {};
  double gfx_y            = {};
  double gonioOut$rEnd    = {};
  double gonioOut$rStart  = {};
  double gonioOut$recPtr  = {};
  double gonioOut$samples = {};
  double gonioOut$thisUI  = {};
#endif
  double hp$a1 {};
  double hp$a2 {};
  double hp$a3 {};
  double hp$k {};
  double hp2$a1 {};
  double hp2$a2 {};
  double hp2$a3 {};
  double hp2$k {};
  double lowCutFreq {};
#if 0
  double newUI            = {};
  double scaling          = {};
#endif
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double buffer$frac {};
  double buffer$readptr {};
  double forceUpdate {};
  double hp$ic10eq {};
  double hp$ic11eq {};
  double hp$ic12eq {};
  double hp$ic1eq {};
  double hp$ic2eq {};
  double hp$ic3eq {};
  double hp$ic4eq {};
  double hp$ic5eq {};
  double hp$ic6eq {};
  double hp$ic7eq {};
  double hp$ic8eq {};
  double hp$ic9eq {};
  double hp2$ic10eq {};
  double hp2$ic11eq {};
  double hp2$ic12eq {};
  double hp2$ic1eq {};
  double hp2$ic2eq {};
  double hp2$ic3eq {};
  double hp2$ic4eq {};
  double hp2$ic5eq {};
  double hp2$ic6eq {};
  double hp2$ic7eq {};
  double hp2$ic8eq {};
  double hp2$ic9eq {};
  double sideLevel {};
  double tDelay {};
  //----------------------------------------------------------------------------
  // added while porting
  float sample_rate = 44100.f;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    sample_rate = pc.get_sample_rate();
    delta_t     = 44. / jsfx_specialvar_get_srate();
#if 0
    MinDelay     = 1.;
    MaxDelay     = 25.;
    MinStrength  = 0.;
    MaxStrength  = 1.;
    MinLowCut    = 0.;
    MaxLowCut    = 1.;
    MinSideLevel = 0.;
    MaxSideLevel = 2.;
    fontface     = 0.; /* jsfx2cpp strings unsupported: was: "Arial" */
#endif
#if 0
    bufferDist = 65536.;
#else
    bufferDist = (samplerate * 1000.f) / get_parameter (delay_tag {}).max;
    mempool.resize (bufferDist);
    memset (mempool.data(), 0, sizeof mempool[0] * mempool.size());
#endif
    delayBuf1 = 0.;
    init$initBuffer (
      delayBuf1,
      delayBuf1 + bufferDist - 1.,
      buffer$scopebuffer,
      buffer$scopebuffermax,
      buffer$scopeptr);
#if 0
    init$initGonioMeter (
      bufferDist,
      1000.,
      gonioOut$samples,
      gonioOut$rStart,
      gonioOut$rEnd,
      gonioOut$recPtr,
      gonioOut$thisUI);
    ;
#endif
  }
//----------------------------------------------------------------------------
#if 0
  void slider() { init$sliderUpdate(); }
#endif
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double avg  = 0.;
    double rb   = 0.;
    double side = 0.;
    delta_t     = 44. / jsfx_specialvar_get_srate();

    for (int i = 0; i < samples; ++i) {
      auto& spl0 = outs[0][i];
      auto& spl1 = outs[1][i];
      spl0       = ins[0][i];
      spl1       = ins[1][i];

      tDelay += chDelay * delta_t - delta_t * tDelay;
      init$setOffset (
        tDelay,
        buffer$readptr,
        buffer$scopeptr,
        buffer$frac,
        buffer$scopebuffer,
        buffer$scopebuffermax);
#if 0
      if (eel2_eq (forceUpdate, 1.)) {
        init$sliderUpdate();
        forceUpdate = 0.;
      }
#endif
      avg  = .5 * (spl0 + spl1);
      side = .5 * (spl0 - spl1);

      init$updateBuffer (
        avg, buffer$scopeptr, buffer$scopebuffermax, buffer$scopebuffer);
      rb = init$readBuffer (
        buffer$readptr, buffer$scopebuffermax, buffer$scopebuffer, buffer$frac);
      rb = block$eval_linearSVF_HP (
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
#if 0
      side = block$eval_linearSVF_HP (
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
        hp2$ic11eq);
#endif
      spl0 = avg + get_slider_Strength() * rb + sideLevel * side;
      spl1 = avg - get_slider_Strength() * rb - sideLevel * side;
#if 0
      if (eel2_eq (get_slider_CheckMono(), 1.)) {
        spl1 = .5 * (spl0 + spl1);
        (spl0 = spl1);
      }
      init$feedGonio (
        spl0, spl1, gonioOut$recPtr, gonioOut$rEnd, gonioOut$rStart);
      ;
#endif
    }
  }

private:
  // functions for section "init"
#if 0 // GFX pollution
  //----------------------------------------------------------------------------
  double init$clearBuffer (double& scopebuffer, double& scopeptr)
  {
    jsfx_memset (scopebuffer, 0., MAXBUFFERSIZE);
    scopeptr = scopebuffer;
    return scopeptr;
  }
  //----------------------------------------------------------------------------
  double init$drawGonioMeter (
    double  r,
    double  g,
    double  b,
    double  a,
    double& recPtr,
    double& x,
    double& w,
    double& y,
    double& h,
    double& samples,
    double& dgonio_zoom,
    double& rEnd,
    double& rStart)
  {
    double cy     = 0.;
    double radius = 0.;
    double s1     = 0.;
    double angle  = 0.;
    double xx     = 0.;
    double horiz  = 0.;
    double sign1  = 0.;
    double step   = 0.;
    double s0     = 0.;
    double yy     = 0.;
    double ptr    = 0.;
    double sign0  = 0.;
    double cx     = 0.;
    double vert   = 0.;
    ptr           = recPtr + 2.;
    cx            = x + 0.5 * w;
    cy            = y + 0.5 * h;
    // gfx_set (r, g, b, a) /*TODO: unknown call */;
    return [&] {
      double $$loop_ret_0 = {};
      for (int $$i = 0, $$end = std::max (0, (int) (samples)); $$i < $$end;
           ++$$i) {
        s0 = mempool[(size_t) (ptr)];
        (ptr += 1.);
        s1    = mempool[(size_t) (ptr)];
        sign0 = jsfx_sign (s0);
        sign1 = jsfx_sign (s1);
        angle = atan (s0 / s1) /*TODO: unknown call */;
        if (
          (eel2_eq (sign0, 1.) && eel2_eq (sign1, -1.))
          || (eel2_eq (sign0, -1.) && eel2_eq (sign1, -1.))) {
          angle += 3.141592654;
        }
        if ((eel2_eq (sign0, -1.) && eel2_eq (sign1, 1.))) {
          angle += 6.283185307;
        }
        if (eel2_eq (s1, 0.)) {
          if (s0 > 0.) {
            angle = 1.570796327;
          }
          else {
            angle = 4.71238898;
          }
        }
        if (eel2_eq (s0, 0.)) {
          if (s1 > 0.) {
            angle = 0.;
          }
          else {
            angle = 3.141592654;
          }
        }
        radius = (1. + dgonio_zoom)
          * std::sqrt (
                   sqr (s0) /*TODO: unknown call */
                   + sqr (s1) /*TODO: unknown call */);
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
          if (ptr > rEnd) {
            ptr = rStart;
            return ptr;
          }
        }();
      }
      return $$loop_ret_0;
    }();
  }
  //----------------------------------------------------------------------------
  double init$drawGonioTop (
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
      0.;
      /* jsfx2cpp strings unsupported: was: "L" */) /*TODO: unknown call */;
    gfx_measurestr (0.; /* jsfx2cpp strings unsupported: was: "R" */,
                        fw,
                        fh) /*TODO: unknown call */;
    gfx_x = x + r * w + .5 * w - .5 * fw;
    gfx_y = y - r * w + .5 * h + 2. - .5 * fh;
    gfx_printf (
      0.;
      /* jsfx2cpp strings unsupported: was: "R" */) /*TODO: unknown call */;
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
      0.;
      /* jsfx2cpp strings unsupported: was: "M" */) /*TODO: unknown call */;
  }
  //----------------------------------------------------------------------------
  double init$feedGonio (
    double  l1,
    double  r1,
    double& recPtr,
    double& rEnd,
    double& rStart)
  {
    mempool[(size_t) (recPtr)]                 = l1;
    mempool[(size_t) ((recPtr = recPtr + 1.))] = r1;
    mempool[(size_t) ((recPtr = recPtr + 1.))] = r1;
    recPtr                                     = [&] {
      if ((recPtr + 1.) >= rEnd) {
        return rStart;
      }
      else {
        return recPtr + 1.;
      }
    }();
    return recPtr;
  }
  //----------------------------------------------------------------------------
#endif
  double init$initBuffer (
    double  scopebuffer_in,
    double  scopebuffermax_in,
    double& scopebuffer,
    double& scopebuffermax,
    double& scopeptr)
  {
    scopebuffer    = scopebuffer_in;
    scopebuffermax = scopebuffermax_in;
    if (scopeptr < scopebuffer) {
      scopeptr = scopebuffer;
    }
    else {
      if (scopeptr > scopebuffermax) {
        scopeptr = scopebuffer;
      }
    }
    return scopeptr;
  }
#if 0
  //----------------------------------------------------------------------------
  double init$initGonioMeter (
    double  in,
    double  sampleCount,
    double& samples,
    double& rStart,
    double& rEnd,
    double& recPtr,
    double& thisUI)
  {
    samples = sampleCount;
    rStart  = in;
    rEnd    = in + sampleCount - 1.;
    recPtr  = in;
    newUI += 1.;
    thisUI = newUI;
    return thisUI;
  }
#endif // GFX_STUFF
  //----------------------------------------------------------------------------
  double init$init_linearSVF (
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
  double init$readBuffer (
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
  double init$setOffset (
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
    if (readptr < scopebuffer) {
      readptr += (scopebuffermax - scopebuffer + 1);
    };
    return readptr;
  }
#if 0
  //----------------------------------------------------------------------------
  double init$setWindow (
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
#if 0
  double init$sliderUpdate()
  {
    chDelay = get_slider_Delay() * jsfx_specialvar_get_srate() / 2000.;
    lowCutFreq
      = std::exp ((1. - get_slider_LowCut()) * std::log (20. / 22050.));
    init$init_linearSVF (lowCutFreq, 0., hp$k, hp$a1, hp$a2, hp$a3);
    hp2$k       = hp$k;
    hp2$a1      = hp$a1;
    hp2$a2      = hp$a2;
    hp2$a3      = hp$a3;
    displayFreq = .5 * jsfx_specialvar_get_srate() * lowCutFreq;
    return displayFreq;
  }
#endif
  //----------------------------------------------------------------------------
  double init$updateBuffer (
    double  M,
    double& scopeptr,
    double& scopebuffermax,
    double& scopebuffer)
  {
    mempool[(size_t) (scopeptr)] = M;
    scopeptr += 1.;
    if (scopeptr > scopebuffermax) {
      scopeptr = scopebuffer;
    }
    return M;
  }
  // functions for section "block"
  //----------------------------------------------------------------------------
  double block$eval_linearSVF_HP (
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
    double v1 = 0.;
    double v3 = 0.;
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
};
}} // namespace artv::saike
