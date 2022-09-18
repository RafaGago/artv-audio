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

struct huge_booty {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::exciter;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
private:
  // definitions for environment function calls
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
  void set_slider1_slider (float v)
  {
    // Original slider line: slider1:0<0,100,0.1>Mix (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mix_tag {};

  void set (mix_tag, float v)
  {
    // Original slider line: slider1:0<0,100,0.1>Mix (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }

  static constexpr auto get_parameter (mix_tag)
  {
    // Original slider line: slider1:0<0,100,0.1>Mix (%)
    return float_param ("", 0.0, 100.0, 100.0, 0.1);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:0<0,100,0.1>Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct drive_tag {};

  void set (drive_tag, float v)
  {
    // Original slider line: slider2:0<0,100,0.1>Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (drive_tag)
  {
    // Original slider line: slider2:0<0,100,0.1>Drive (%)
    return float_param ("%", 0.0, 100.0, 30.0, 0.1);
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:100<20,200,1>Frequency (Hz)
    // Range: min:20.0, max:200.0, default: 100.0, step: 1.0
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct frequency_tag {};

  void set (frequency_tag, float v)
  {
    // Original slider line: slider3:100<20,200,1>Frequency (Hz)
    // Range: min:20.0, max:200.0, default: 100.0, step: 1.0
    v = midi_note_to_hz (v);
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (frequency_tag)
  {
    // Original slider line: slider3:100<20,200,1>Frequency (Hz)
    return frequency_parameter (20.0, 200.0, 100.0);
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<mix_tag, drive_tag, frequency_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double a01;
  double a03;
  double a1;
  double a11;
  double a13;
  double a21;
  double a23;
  double a3;
  double alpha1;
  double alpha3;
  double b01;
  double b03;
  double b11;
  double b13;
  double b21;
  double b23;
  double cosw01;
  double cosw03;
  double freq;
  double freq1;
  double freq3;
  double gain1;
  double gain3;
  double q1;
  double q3;
  double s1;
  double s3;
  double sinw01;
  double sinw03;
  double slider1;
  double slider2;
  double slider3;
  double w01;
  double w03;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    a01     = 0;
    a03     = 0;
    a1      = 0;
    a11     = 0;
    a13     = 0;
    a21     = 0;
    a23     = 0;
    a3      = 0;
    alpha1  = 0;
    alpha3  = 0;
    b01     = 0;
    b03     = 0;
    b11     = 0;
    b13     = 0;
    b21     = 0;
    b23     = 0;
    cosw01  = 0;
    cosw03  = 0;
    freq    = 0;
    freq1   = 0;
    freq3   = 0;
    gain1   = 0;
    slider1 = 0.0;
    gain3   = 0;
    q1      = 0;
    q3      = 0;
    s1      = 0;
    s3      = 0;
    sinw01  = 0;
    sinw03  = 0;
    slider1 = 0;
    slider2 = 0;
    slider3 = 0;
    w01     = 0;
    w03     = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double drive1;
  double drive2;
  double mix;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    drive1 = 0;
    drive2 = 0;
    mix    = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double xl11;
  double xl13;
  double xl21;
  double xl23;
  double xr11;
  double xr13;
  double xr21;
  double xr23;
  double yl11;
  double yl13;
  double yl21;
  double yl23;
  double yr11;
  double yr13;
  double yr21;
  double yr23;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    xl11 = 0;
    xl13 = 0;
    xl21 = 0;
    xl23 = 0;
    xr11 = 0;
    xr13 = 0;
    xr21 = 0;
    xr23 = 0;
    yl11 = 0;
    yl13 = 0;
    yl21 = 0;
    yl23 = 0;
    yr11 = 0;
    yr13 = 0;
    yr21 = 0;
    yr23 = 0;
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
#if 0
    slider1 = 0.0;
#else
    slider1 = 100.0;
#endif
    slider2 = 0.0;
    slider3 = 100.0;
    gain1   = 0.;
    freq1   = freq;
    a1      = 1.;
    s1      = 1.;
    q1      = 1. / std::sqrt ((a1 + 1. / a1) * (1. / s1 - 1.) + 2.);
    w01     = 2. * 3.141592653589793 * freq1 / jsfx_specialvar_get_srate();
    cosw01  = std::cos (w01);
    sinw01  = std::sin (w01);
    alpha1  = sinw01 / (2. * q1);
    b01     = (1. - cosw01) / 2.;
    b11     = (1. - cosw01);
    b21     = (1. - cosw01) / 2.;
    a01     = 1. + alpha1;
    a11     = -2. * cosw01;
    a21     = 1. - alpha1;
    b01 /= a01;
    b11 /= a01;
    b21 /= a01;
    a11 /= a01;
    a21 /= a01;
    gain3  = 0.;
    freq3  = freq;
    a3     = eel2_pow (10., (gain3 / 40.));
    s3     = 1.;
    q3     = 1. / std::sqrt ((a3 + 1. / a3) * (1. / s3 - 1.) + 2.);
    w03    = 2. * 3.141592653589793 * freq3 / jsfx_specialvar_get_srate();
    cosw03 = std::cos (w03);
    sinw03 = std::sin (w03);
    alpha3 = sinw03 / (2. * q3);
    b03    = (1. + cosw03) / 2.;
    b13    = -(1. + cosw03);
    b23    = (1. + cosw03) / 2.;
    a03    = 1. + alpha3;
    a13    = -2. * cosw03;
    a23    = 1. - alpha3;
    b03 /= a03;
    b13 /= a03;
    b23 /= a03;
    a13 /= a03;
    a23 /= a03;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double drive = 0.;
    double mix1  = 0.;
    mix          = (slider1 * 1.5) / 100.;
    drive        = slider2 / 100.;
    mix1         = 1. - mix;
    drive1       = 1. / (1. - (drive / 2.));
    drive2       = drive / 2.;
    freq         = slider3;
    gain1        = 0.;
    freq1        = freq;
    a1           = 1.;
    s1           = 1.;
    q1           = 1. / std::sqrt ((a1 + 1. / a1) * (1. / s1 - 1.) + 2.);
    w01          = 2. * 3.141592653589793 * freq1 / jsfx_specialvar_get_srate();
    cosw01       = std::cos (w01);
    sinw01       = std::sin (w01);
    alpha1       = sinw01 / (2. * q1);
    b01          = (1. - cosw01) / 2.;
    b11          = (1. - cosw01);
    b21          = (1. - cosw01) / 2.;
    a01          = 1. + alpha1;
    a11          = -2. * cosw01;
    a21          = 1. - alpha1;
    b01 /= a01;
    b11 /= a01;
    b21 /= a01;
    a11 /= a01;
    a21 /= a01;
    gain3  = 0.;
    freq3  = freq;
    a3     = eel2_pow (10., (gain3 / 40.));
    s3     = 1.;
    q3     = 1. / std::sqrt ((a3 + 1. / a3) * (1. / s3 - 1.) + 2.);
    w03    = 2. * 3.141592653589793 * freq3 / jsfx_specialvar_get_srate();
    cosw03 = std::cos (w03);
    sinw03 = std::sin (w03);
    alpha3 = sinw03 / (2. * q3);
    b03    = (1. + cosw03) / 2.;
    b13    = -(1. + cosw03);
    b23    = (1. + cosw03) / 2.;
    a03    = 1. + alpha3;
    a13    = -2. * cosw03;
    a23    = 1. - alpha3;
    b03 /= a03;
    b13 /= a03;
    b23 /= a03;
    a13 /= a03;
    a23 /= a03;
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));
    double dry0  = 0.;
    double dry1  = 0.;
    double ospl0 = 0.;
    double ospl1 = 0.;

    double wet0 = 0.;
    double wet1 = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];
      dry0       = spl0;
      dry1       = spl1;
      ospl0      = dry0;
      dry0  = b01 * dry0 + b11 * xl11 + b21 * xl21 - a11 * yl11 - a21 * yl21;
      xl21  = xl11;
      xl11  = ospl0;
      yl21  = yl11;
      yl11  = dry0;
      ospl1 = dry1;
      dry1  = b01 * dry1 + b11 * xr11 + b21 * xr21 - a11 * yr11 - a21 * yr21;
      xr21  = xr11;
      xr11  = ospl1;
      yr21  = yr11;
      yr11  = dry1;
      wet0  = drive1 * dry0 * (1. - std::abs (dry0 * drive2));
      wet1  = drive1 * dry1 * (1. - std::abs (dry1 * drive2));
      ospl0 = wet0;
      wet0  = b03 * wet0 + b13 * xl13 + b23 * xl23 - a13 * yl13 - a23 * yl23;
      xl23  = xl13;
      xl13  = ospl0;
      yl23  = yl13;
      yl13  = wet0;
      ospl1 = wet1;
      wet1  = b03 * wet1 + b13 * xr13 + b23 * xr23 - a13 * yr13 - a23 * yr23;
      xr23  = xr13;
      xr13  = ospl1;
      yr23  = yr13;
      yr13  = wet1;
      spl0 += mix * wet0;
      spl1 += mix * wet1;
      ;
    }
  }
}; /* jsfx_process */
}} // namespace artv::sstillwell
