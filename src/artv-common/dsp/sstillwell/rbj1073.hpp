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

struct rbj1073 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::eq;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // definitions for environment function calls
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  static bool   eel2_ne (double lhs, double rhs) { return !eel2_eq (lhs, rhs); }
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
    // Original slider line: slider1:0<0,4,1{Off,50,80,160,300}>HPF
    // Range: min:0.0, max:4.0, default: 0.0, step: 1.0
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct hpf_tag {};

  void set (hpf_tag, float v)
  {
    // Original slider line: slider1:0<0,4,1{Off,50,80,160,300}>HPF
    // Range: min:0.0, max:4.0, default: 0.0, step: 1.0
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }

  static constexpr auto get_parameter (hpf_tag)
  {
    // Original slider line: slider1:0<0,4,1{Off,50,80,160,300}>HPF
    return choice_param (
      0, make_cstr_array ("Off", "50Hz", "80Hz", "160Hz", "300Hz"));
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:0<0,4,1{Off,35,60,110,220}>Low Shelf (Hz)
    // Range: min:0.0, max:4.0, default: 0.0, step: 1.0
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct low_shelf_tag {};

  void set (low_shelf_tag, float v)
  {
    // Original slider line: slider2:0<0,4,1{Off,35,60,110,220}>Low Shelf (Hz)
    // Range: min:0.0, max:4.0, default: 0.0, step: 1.0
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (low_shelf_tag)
  {
    // Original slider line: slider2:0<0,4,1{Off,35,60,110,220}>Low Shelf (Hz)
    return choice_param (
      0, make_cstr_array ("Off", "35Hz", "60Hz", "110Hz", "220Hz"));
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<-20,20,0.1>Low Boost/Cut (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct low_gain_tag {};

  void set (low_gain_tag, float v)
  {
    // Original slider line: slider3:0<-20,20,0.1>Low Boost/Cut (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (low_gain_tag)
  {
    // Original slider line: slider3:0<-20,20,0.1>Low Boost/Cut (dB)
    return float_param ("dB", -20.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider4_slider (float v)
  {
    // Original slider line: slider4:0<0,5,1{360,700,1.6k,3.2k,4.8k,7.2k}>Mid
    // Freq (Hz) Range: min:0.0, max:5.0, default: 0.0, step: 1.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mid_freq_tag {};

  void set (mid_freq_tag, float v)
  {
    // Original slider line: slider4:0<0,5,1{360,700,1.6k,3.2k,4.8k,7.2k}>Mid
    // Freq (Hz) Range: min:0.0, max:5.0, default: 0.0, step: 1.0
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }

  static constexpr auto get_parameter (mid_freq_tag)
  {
    // Original slider line: slider4:0<0,5,1{360,700,1.6k,3.2k,4.8k,7.2k}>Mid
    // Freq (Hz)
    return choice_param (
      0,
      make_cstr_array (
        "Off", "360Hz", "700Hz", "1.6kHz", "3.2kHz", "4.8kHz", "7.2kHz"));
  }

#endif
#if 0
  void set_slider5_slider (float v)
  {
    // Original slider line: slider5:0<-20,20,0.1>Mid Boost/Cut (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mid_gain_tag {};

  void set (mid_gain_tag, float v)
  {
    // Original slider line: slider5:0<-20,20,0.1>Mid Boost/Cut (dB)
    // Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }

  static constexpr auto get_parameter (mid_gain_tag)
  {
    // Original slider line: slider5:0<-20,20,0.1>Mid Boost/Cut (dB)
    return float_param ("", -20.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider6_slider (float v)
  {
    // Original slider line: slider6:0<-20,20,0.1>High Shelf (12k) Boost/Cut
    // (dB) Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct high_12k_gain_tag {};

  void set (high_12k_gain_tag, float v)
  {
    // Original slider line: slider6:0<-20,20,0.1>High Shelf (12k) Boost/Cut
    // (dB) Range: min:-20.0, max:20.0, default: 0.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }

  static constexpr auto get_parameter (high_12k_gain_tag)
  {
    // Original slider line: slider6:0<-20,20,0.1>High Shelf (12k) Boost/Cut
    // (dB)
    return float_param ("", -20.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider7_slider (float v)
  {
    // Original slider line: slider7:0<-20,10,0.1>Gain (dB)
    // Range: min:-20.0, max:10.0, default: 0.0, step: 0.1
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct gain_tag {};

  void set (gain_tag, float v)
  {
    // Original slider line: slider7:0<-20,10,0.1>Gain (dB)
    // Range: min:-20.0, max:10.0, default: 0.0, step: 0.1
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }

  static constexpr auto get_parameter (gain_tag)
  {
    // Original slider line: slider7:0<-20,10,0.1>Gain (dB)
    return float_param ("", -20.0, 10.0, 0.0, 0.1);
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    hpf_tag,
    low_shelf_tag,
    low_gain_tag,
    mid_freq_tag,
    mid_gain_tag,
    high_12k_gain_tag,
    gain_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double a01;
  double a03;
  double a05;
  double a07;
  double a1;
  double a11;
  double a13;
  double a15;
  double a17;
  double a21;
  double a23;
  double a25;
  double a27;
  double a3;
  double a5;
  double a7;
  double alpha1;
  double alpha3;
  double alpha5;
  double alpha7;
  double b01;
  double b03;
  double b05;
  double b07;
  double b11;
  double b13;
  double b15;
  double b17;
  double b21;
  double b23;
  double b25;
  double b27;
  double cosw01;
  double cosw03;
  double cosw05;
  double cosw07;
  double freq1;
  double freq3;
  double freq5;
  double freq7;
  double gain;
  double gain1;
  double gain3;
  double gain5;
  double gain7;
  double hpf;
  double lshelf;
  double q1;
  double q3;
  double q5;
  double q7;
  double s1;
  double s3;
  double s7;
  double sinw01;
  double sinw03;
  double sinw05;
  double sinw07;
  double slider1;
  double slider2;
  double slider3;
  double slider4;
  double slider5;
  double slider6;
  double slider7;
  double w01;
  double w03;
  double w05;
  double w07;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    a01     = 0;
    a03     = 0;
    a05     = 0;
    a07     = 0;
    a1      = 0;
    a11     = 0;
    a13     = 0;
    a15     = 0;
    a17     = 0;
    a21     = 0;
    a23     = 0;
    a25     = 0;
    a27     = 0;
    a3      = 0;
    a5      = 0;
    a7      = 0;
    alpha1  = 0;
    alpha3  = 0;
    alpha5  = 0;
    alpha7  = 0;
    b01     = 0;
    b03     = 0;
    b05     = 0;
    b07     = 0;
    b11     = 0;
    b13     = 0;
    b15     = 0;
    b17     = 0;
    b21     = 0;
    b23     = 0;
    b25     = 0;
    b27     = 0;
    cosw01  = 0;
    cosw03  = 0;
    cosw05  = 0;
    cosw07  = 0;
    freq1   = 0;
    freq3   = 0;
    freq5   = 0;
    freq7   = 0;
    gain    = 0;
    gain1   = 0;
    gain3   = 0;
    gain5   = 0;
    gain7   = 0;
    hpf     = 0;
    lshelf  = 0;
    q1      = 0;
    q3      = 0;
    q5      = 0;
    q7      = 0;
    s1      = 0;
    s3      = 0;
    s7      = 0;
    sinw01  = 0;
    sinw03  = 0;
    sinw05  = 0;
    sinw07  = 0;
    slider1 = 0;
    slider2 = 0;
    slider3 = 0;
    slider4 = 0;
    slider5 = 0;
    slider6 = 0;
    slider7 = 0;
    w01     = 0;
    w03     = 0;
    w05     = 0;
    w07     = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double xl11;
  double xl13;
  double xl15;
  double xl17;
  double xl21;
  double xl23;
  double xl25;
  double xl27;
  double xr11;
  double xr13;
  double xr15;
  double xr17;
  double xr21;
  double xr23;
  double xr25;
  double xr27;
  double yl11;
  double yl13;
  double yl15;
  double yl17;
  double yl21;
  double yl23;
  double yl25;
  double yl27;
  double yr11;
  double yr13;
  double yr15;
  double yr17;
  double yr21;
  double yr23;
  double yr25;
  double yr27;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    xl11 = 0;
    xl13 = 0;
    xl15 = 0;
    xl17 = 0;
    xl21 = 0;
    xl23 = 0;
    xl25 = 0;
    xl27 = 0;
    xr11 = 0;
    xr13 = 0;
    xr15 = 0;
    xr17 = 0;
    xr21 = 0;
    xr23 = 0;
    xr25 = 0;
    xr27 = 0;
    yl11 = 0;
    yl13 = 0;
    yl15 = 0;
    yl17 = 0;
    yl21 = 0;
    yl23 = 0;
    yl25 = 0;
    yl27 = 0;
    yr11 = 0;
    yr13 = 0;
    yr15 = 0;
    yr17 = 0;
    yr21 = 0;
    yr23 = 0;
    yr25 = 0;
    yr27 = 0;
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
    slider1 = 0.0;
    slider2 = 0.0;
    slider3 = 0.0;
    slider4 = 0.0;
    slider5 = 0.0;
    slider6 = 0.0;
    slider7 = 0.0;
    hpf     = 0.;
    gain1   = 0.;
    freq1   = 50.;
    a1      = 1.;
    s1      = 1.;
    q1      = 1. / std::sqrt ((a1 + 1. / a1) * (1. / s1 - 1.) + 2.);
    w01     = 2. * 3.141592653589793 * freq1 / jsfx_specialvar_get_srate();
    cosw01  = std::cos (w01);
    sinw01  = std::sin (w01);
    alpha1  = sinw01 / (2. * q1);
    b01     = (1. + cosw01) / 2.;
    b11     = -(1. + cosw01);
    b21     = (1. + cosw01) / 2.;
    a01     = 1. + alpha1;
    a11     = -2. * cosw01;
    a21     = 1. - alpha1;
    b01 /= a01;
    b11 /= a01;
    b21 /= a01;
    a11 /= a01;
    a21 /= a01;
    lshelf = 0.;
    gain3  = 0.;
    freq3  = 35.;
    a3     = eel2_pow (10., (gain3 / 40.));
    s3     = 2.;
    q3     = 1. / std::sqrt ((a3 + 1. / a3) * (1. / s3 - 1.) + 2.);
    w03    = 2. * 3.141592653589793 * freq3 / jsfx_specialvar_get_srate();
    cosw03 = std::cos (w03);
    sinw03 = std::sin (w03);
    alpha3 = sinw03 / (2. * q3);
    b03 = a3 * ((a3 + 1.) - (a3 - 1.) * cosw03 + 2. * std::sqrt (a3) * alpha3);
    b13 = 2. * a3 * ((a3 - 1.) - (a3 + 1.) * cosw03);
    b23 = a3 * ((a3 + 1.) - (a3 - 1.) * cosw03 - 2. * std::sqrt (a3) * alpha3);
    a03 = (a3 + 1.) + (a3 - 1.) * cosw03 + 2. * std::sqrt (a3) * alpha3;
    a13 = -2. * ((a3 - 1.) + (a3 + 1.) * cosw03);
    a23 = (a3 + 1.) + (a3 - 1.) * cosw03 - 2. * std::sqrt (a3) * alpha3;
    b03 /= a03;
    b13 /= a03;
    b23 /= a03;
    a13 /= a03;
    a23 /= a03;
    gain5  = 0.;
    freq5  = 360.;
    a5     = eel2_pow (10., (gain5 / 20.));
    q5     = 1.4;
    w05    = 2. * 3.141592653589793 * freq5 / jsfx_specialvar_get_srate();
    cosw05 = std::cos (w05);
    sinw05 = std::sin (w05);
    alpha5 = sinw05 / (2. * q5);
    b05    = 1. + alpha5 * a5;
    b15    = -2. * cosw05;
    b25    = 1. - alpha5 * a5;
    a05    = 1. + alpha5 / a5;
    a15    = -2. * cosw05;
    a25    = 1. - alpha5 / a5;
    b05 /= a05;
    b15 /= a05;
    b25 /= a05;
    a15 /= a05;
    a25 /= a05;
    gain7  = 0.;
    freq7  = 12000.;
    a7     = eel2_pow (10., (gain7 / 40.));
    s7     = 0.3;
    q7     = 1. / std::sqrt ((a7 + 1. / a7) * (1. / s7 - 1.) + 2.);
    w07    = 2. * 3.141592653589793 * freq7 / jsfx_specialvar_get_srate();
    cosw07 = std::cos (w07);
    sinw07 = std::sin (w07);
    alpha7 = sinw07 / (2. * q7);
    b07 = a7 * ((a7 + 1.) + (a7 - 1.) * cosw07 + 2. * std::sqrt (a7) * alpha7);
    b17 = -2. * a7 * ((a7 - 1.) + (a7 + 1.) * cosw07);
    b27 = a7 * ((a7 + 1.) + (a7 - 1.) * cosw07 - 2. * std::sqrt (a7) * alpha7);
    a07 = (a7 + 1.) - (a7 - 1.) * cosw07 + 2. * std::sqrt (a7) * alpha7;
    a17 = 2. * ((a7 - 1.) - (a7 + 1.) * cosw07);
    a27 = (a7 + 1.) - (a7 - 1.) * cosw07 - 2. * std::sqrt (a7) * alpha7;
    b07 /= a07;
    b17 /= a07;
    b27 /= a07;
    a17 /= a07;
    a27 /= a07;
    gain = 1.;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    freq1 = [&] {
      if (eel2_eq (slider1, 0.)) {
        return 50.;
      }
      else {
        return [&] {
          if (eel2_eq (slider1, 1.)) {
            return 50.;
          }
          else {
            return [&] {
              if (eel2_eq (slider1, 2.)) {
                return 80.;
              }
              else {
                return [&] {
                  if (eel2_eq (slider1, 3.)) {
                    return 160.;
                  }
                  else {
                    return 300.;
                  }
                }();
              }
            }();
          }
        }();
      }
    }();
    freq3 = [&] {
      if (eel2_eq (slider2, 0.)) {
        return 35.;
      }
      else {
        return [&] {
          if (eel2_eq (slider2, 1.)) {
            return 35.;
          }
          else {
            return [&] {
              if (eel2_eq (slider2, 2.)) {
                return 60.;
              }
              else {
                return [&] {
                  if (eel2_eq (slider2, 3.)) {
                    return 110.;
                  }
                  else {
                    return 220.;
                  }
                }();
              }
            }();
          }
        }();
      }
    }();
    gain3 = slider3;
    freq5 = [&] {
      if (eel2_eq (slider4, 0.)) {
        return 360.;
      }
      else {
        return [&] {
          if (eel2_eq (slider4, 1.)) {
            return 700.;
          }
          else {
            return [&] {
              if (eel2_eq (slider4, 2.)) {
                return 1600.;
              }
              else {
                return [&] {
                  if (eel2_eq (slider4, 3.)) {
                    return 3200.;
                  }
                  else {
                    return [&] {
                      if (eel2_eq (slider4, 4.)) {
                        return 4800.;
                      }
                      else {
                        return 7200.;
                      }
                    }();
                  }
                }();
              }
            }();
          }
        }();
      }
    }();
    gain5 = slider5;
    gain7 = slider6;
    gain  = eel2_pow (10., (slider7 / 20.));
    if (eel2_eq (slider1, 0.)) {
      hpf = 0.;
    }
    else {
      hpf = 1.;
    }
    if (eel2_eq (slider2, 0.)) {
      lshelf = 0.;
    }
    else {
      lshelf = 1.;
    }
    a1     = 1.;
    s1     = 1.;
    q1     = 1. / std::sqrt ((a1 + 1. / a1) * (1. / s1 - 1.) + 2.);
    w01    = 2. * 3.141592653589793 * freq1 / jsfx_specialvar_get_srate();
    cosw01 = std::cos (w01);
    sinw01 = std::sin (w01);
    alpha1 = sinw01 / (2. * q1);
    b01    = (1. + cosw01) / 2.;
    b11    = -(1. + cosw01);
    b21    = (1. + cosw01) / 2.;
    a01    = 1. + alpha1;
    a11    = -2. * cosw01;
    a21    = 1. - alpha1;
    b01 /= a01;
    b11 /= a01;
    b21 /= a01;
    a11 /= a01;
    a21 /= a01;
    a3     = eel2_pow (10., (gain3 / 40.));
    s3     = 2.;
    q3     = 1. / std::sqrt ((a3 + 1. / a3) * (1. / s3 - 1.) + 2.);
    w03    = 2. * 3.141592653589793 * freq3 / jsfx_specialvar_get_srate();
    cosw03 = std::cos (w03);
    sinw03 = std::sin (w03);
    alpha3 = sinw03 / (2. * q3);
    b03 = a3 * ((a3 + 1.) - (a3 - 1.) * cosw03 + 2. * std::sqrt (a3) * alpha3);
    b13 = 2. * a3 * ((a3 - 1.) - (a3 + 1.) * cosw03);
    b23 = a3 * ((a3 + 1.) - (a3 - 1.) * cosw03 - 2. * std::sqrt (a3) * alpha3);
    a03 = (a3 + 1.) + (a3 - 1.) * cosw03 + 2. * std::sqrt (a3) * alpha3;
    a13 = -2. * ((a3 - 1.) + (a3 + 1.) * cosw03);
    a23 = (a3 + 1.) + (a3 - 1.) * cosw03 - 2. * std::sqrt (a3) * alpha3;
    b03 /= a03;
    b13 /= a03;
    b23 /= a03;
    a13 /= a03;
    a23 /= a03;
    a5     = eel2_pow (10., (gain5 / 20.));
    q5     = 1.4;
    w05    = 2. * 3.141592653589793 * freq5 / jsfx_specialvar_get_srate();
    cosw05 = std::cos (w05);
    sinw05 = std::sin (w05);
    alpha5 = sinw05 / (2. * q5);
    b05    = 1. + alpha5 * a5;
    b15    = -2. * cosw05;
    b25    = 1. - alpha5 * a5;
    a05    = 1. + alpha5 / a5;
    a15    = -2. * cosw05;
    a25    = 1. - alpha5 / a5;
    b05 /= a05;
    b15 /= a05;
    b25 /= a05;
    a15 /= a05;
    a25 /= a05;
    a7     = eel2_pow (10., (gain7 / 40.));
    freq7  = 12000.;
    s7     = 0.3;
    q7     = 1. / std::sqrt ((a7 + 1. / a7) * (1. / s7 - 1.) + 2.);
    w07    = 2. * 3.141592653589793 * freq7 / jsfx_specialvar_get_srate();
    cosw07 = std::cos (w07);
    sinw07 = std::sin (w07);
    alpha7 = sinw07 / (2. * q7);
    b07 = a7 * ((a7 + 1.) + (a7 - 1.) * cosw07 + 2. * std::sqrt (a7) * alpha7);
    b17 = -2. * a7 * ((a7 - 1.) + (a7 + 1.) * cosw07);
    b27 = a7 * ((a7 + 1.) + (a7 - 1.) * cosw07 - 2. * std::sqrt (a7) * alpha7);
    a07 = (a7 + 1.) - (a7 - 1.) * cosw07 + 2. * std::sqrt (a7) * alpha7;
    a17 = 2. * ((a7 - 1.) - (a7 + 1.) * cosw07);
    a27 = (a7 + 1.) - (a7 - 1.) * cosw07 - 2. * std::sqrt (a7) * alpha7;
    b07 /= a07;
    b17 /= a07;
    b27 /= a07;
    a17 /= a07;
    a27 /= a07;
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint block_samples)
  {
    double ospl0 = 0.;
    double ospl1 = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = chnls[0][$$i];
      auto& spl1 = chnls[1][$$i];
      if (eel2_ne (hpf, 0.)) {
        ospl0 = spl0;
        spl0  = b01 * spl0 + b11 * xl11 + b21 * xl21 - a11 * yl11 - a21 * yl21;
        xl21  = xl11;
        xl11  = ospl0;
        yl21  = yl11;
        yl11  = spl0;
        ospl1 = spl1;
        spl1  = b01 * spl1 + b11 * xr11 + b21 * xr21 - a11 * yr11 - a21 * yr21;
        xr21  = xr11;
        xr11  = ospl1;
        yr21  = yr11;
        yr11  = spl1;
      }
      if (eel2_ne (lshelf, 0.) && eel2_ne (gain3, 0.)) {
        ospl0 = spl0;
        spl0  = b03 * spl0 + b13 * xl13 + b23 * xl23 - a13 * yl13 - a23 * yl23;
        xl23  = xl13;
        xl13  = ospl0;
        yl23  = yl13;
        yl13  = spl0;
        ospl1 = spl1;
        spl1  = b03 * spl1 + b13 * xr13 + b23 * xr23 - a13 * yr13 - a23 * yr23;
        xr23  = xr13;
        xr13  = ospl1;
        yr23  = yr13;
        yr13  = spl1;
      }
      if (eel2_ne (gain5, 0.)) {
        ospl0 = spl0;
        spl0  = b05 * spl0 + b15 * xl15 + b25 * xl25 - a15 * yl15 - a25 * yl25;
        xl25  = xl15;
        xl15  = ospl0;
        yl25  = yl15;
        yl15  = spl0;
        ospl1 = spl1;
        spl1  = b05 * spl1 + b15 * xr15 + b25 * xr25 - a15 * yr15 - a25 * yr25;
        xr25  = xr15;
        xr15  = ospl1;
        yr25  = yr15;
        yr15  = spl1;
      }
      if (eel2_ne (gain7, 0.)) {
        ospl0 = spl0;
        spl0  = b07 * spl0 + b17 * xl17 + b27 * xl27 - a17 * yl17 - a27 * yl27;
        xl27  = xl17;
        xl17  = ospl0;
        yl27  = yl17;
        yl17  = spl0;
        ospl1 = spl1;
        spl1  = b07 * spl1 + b17 * xr17 + b27 * xr27 - a17 * yr17 - a27 * yr27;
        xr27  = xr17;
        xr17  = ospl1;
        yr27  = yr17;
        yr17  = spl1;
      }
      spl0 *= gain;
      spl1 *= gain;
      ;
    }
  }
}; /* jsfx_process */

}} // namespace artv::sstillwell
