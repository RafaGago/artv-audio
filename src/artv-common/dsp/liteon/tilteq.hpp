// Ported from https://github.com/Rcomian/jsfx-vcv
// commit sha: 65bfa7896480f3de7d424e42a6ef2b3737b43a20
#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls
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

struct tilt_eq {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::eq;
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
  static double eel2_pow (double lhs, double rhs)
  {
    return std::pow (lhs, rhs);
  }

  //----------------------------------------------------------------------------
  // stubs for JSFX special variables

  double jsfx_specialvar_get_samplesblock()
  {
    return 0.; /* TODO: stub for getting JSFX var "samplesblock" */
  }

  double jsfx_specialvar_get_srate()
  {
    return 0.; /* TODO: stub for getting JSFX var "srate" */
  }

  //----------------------------------------------------------------------------
public:
  // stubs for sliders
  // slider1:0<0,1,1{Stereo,Mono}>Processing
  double get_slider_slider1()
  {
    // TODO: stub, add code for getting "slider1"
    // Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    return 0.;
  }
  //----------------------------------------------------------------------------
  // slider2:50<0,100,0.05>Center Frequency (Scale)
  double get_slider_slider2()
  {
    // TODO: stub, add code for getting "slider2"
    // Range: min:0.0, max:100.0, default: 50.0, step: 0.05
    return frequencyp;
  }
  float frequencyp = 50.;
  struct frequency_tag {};
  void set (frequency_tag, float v)
  {
    if (v == frequencyp) {
      return;
    }
    frequencyp = v;
    slider();
  }

  static constexpr auto get_parameter (frequency_tag)
  {
    return float_param ("", 0.f, 100.f, 50.f, 0.05f);
  }
  //----------------------------------------------------------------------------
  // slider3:0<-6,6,0.05>Tilt (Low / High) (dB)
  double get_slider_slider3()
  {
    // TODO: stub, add code for getting "slider3"
    // Range: min:-6.0, max:6.0, default: 0.0, step: 0.05
    return tiltp;
  }
  float tiltp = 50.;
  struct tilt_tag {};
  void set (tilt_tag, float v)
  {
    if (v == tiltp) {
      return;
    }
    tiltp = v;
    slider();
  }

  static constexpr auto get_parameter (tilt_tag)
  {
    return float_param ("dB", -6.f, 6.f, 0.f, 0.05f);
  }
  //----------------------------------------------------------------------------
  // slider4:0<-25,25,0.05>Output Gain (dB)
  double get_slider_slider4()
  {
    // TODO: stub, add code for getting "slider4"
    // Range: min:-25.0, max:25.0, default: 0.0, step: 0.05
    return 0.;
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<frequency_tag, tilt_tag>;

private:
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double amp;
  double denorm;
  double pi;
  double sr3;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    amp    = 0;
    denorm = 0;
    pi     = 0;
    sr3    = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double a0;
  double b1;
  double hgain;
  double lgain;
  double mono;
  double outgain;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    a0      = 0;
    b1      = 0;
    hgain   = 0;
    lgain   = 0;
    mono    = 0;
    outgain = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double lp_out;
  double lp_out_r;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    lp_out   = 0;
    lp_out_r = 0;
  }
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    init_init_variables();
    init_slider_variables();
    init_block_variables();
    amp    = 6. / std::log (2.);
    denorm = eel2_pow (10., -30.);
    pi     = 3.141592653589793;
    sr3    = 3. * pc.get_sample_rate();
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double f0          = 0.;
    double g1          = 0.;
    double g2          = 0.;
    double gain        = 0.;
    double gfactor     = 0.;
    double n           = 0.;
    double omega       = 0.;
    double slider2_val = 0.;
    double sx          = 0.;
    mono               = get_slider_slider1();
    gain               = get_slider_slider3();
    outgain            = std::exp (get_slider_slider4() / amp);
    gfactor            = 4.;
    if (gain > 0.) {
      g1 = -gfactor * gain;
      g2 = gain;
    }
    else {
      g1 = -gain;
      g2 = gfactor * gain;
    }
    lgain       = std::exp (g1 / amp) - 1.;
    hgain       = std::exp (g2 / amp) - 1.;
    slider2_val = std::max (std::min (get_slider_slider2(), 100.), 0.);
    sx          = 16. + slider2_val * 1.20103;
    f0          = std::floor (std::exp (sx * std::log (1.059)) * 8.17742);
    omega       = 2. * pi * f0;
    n           = 1. / (sr3 + omega);
    a0          = 2. * omega * n;
    b1          = (sr3 - omega) * n;
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {

    double input    = 0.;
    double input_r  = 0.;
    double output   = 0.;
    double output_r = 0.;
    double spl0     = 0.;
    double spl1     = 0.;

    for (uint i = 0; i < samples; ++i) {
      auto& spl0 = chnls[0][i];
      auto& spl1 = chnls[1][i];
      if (eel2_eq (mono, 1.)) {
        input  = (spl0 + spl1) / 2.;
        lp_out = a0 * input + b1 * lp_out;
        output = input + lgain * lp_out + hgain * (input - lp_out);
        spl1   = output * outgain + denorm;
        spl0   = spl1;
      }
      else {
        input    = spl0;
        lp_out   = a0 * input + b1 * lp_out;
        output   = input + lgain * lp_out + hgain * (input - lp_out);
        spl0     = output * outgain + denorm;
        input_r  = spl1;
        lp_out_r = a0 * input_r + b1 * lp_out_r;
        output_r = input_r + lgain * lp_out_r + hgain * (input_r - lp_out_r);
        spl1     = output_r * outgain + denorm;
      }
    }
  }
}; /* tilt_eq */
}} // namespace artv::liteon
