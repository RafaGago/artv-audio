#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls

// From Reaper 6.27

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

struct _4x4 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::exciter;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  // definitions for environment function calls

  //----------------------------------------------------------------------------
  // stubs for JSFX special variables

  double jsfx_specialvar_get_srate() { return plugcontext->get_sample_rate(); }

  //----------------------------------------------------------------------------
public:
#if 0
  void set_slider1_slider (float v)
  {
    // Original slider line: slider1:0<0,100,0.1>Low Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct low_drive_tag {};

  void set (low_drive_tag, float v)
  {
    // Original slider line: slider1:0<0,100,0.1>Low Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider1) {
      return;
    }
    slider1 = v;
    slider();
  }

  static constexpr auto get_parameter (low_drive_tag)
  {
    // Original slider line: slider1:0<0,100,0.1>Low Drive (%)
    return float_param ("%", 0.0, 100.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider2_slider (float v)
  {
    // Original slider line: slider2:0<-12,12,0.1>Low Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct low_gain_tag {};

  void set (low_gain_tag, float v)
  {
    // Original slider line: slider2:0<-12,12,0.1>Low Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider2) {
      return;
    }
    slider2 = v;
    slider();
  }

  static constexpr auto get_parameter (low_gain_tag)
  {
    // Original slider line: slider2:0<-12,12,0.1>Low Gain (dB)
    return float_param ("dB", -12.0, 12.0, 0.0, 0.3);
  }

#endif
#if 0
  void set_slider3_slider (float v)
  {
    // Original slider line: slider3:0<0,100,0.1>Mid Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mid_drive_tag {};

  void set (mid_drive_tag, float v)
  {
    // Original slider line: slider3:0<0,100,0.1>Mid Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider3) {
      return;
    }
    slider3 = v;
    slider();
  }

  static constexpr auto get_parameter (mid_drive_tag)
  {
    // Original slider line: slider3:0<0,100,0.1>Mid Drive (%)
    return float_param ("%", 0.0, 100.0, 0.0, 0.3);
  }

#endif
#if 0
  void set_slider4_slider (float v)
  {
    // Original slider line: slider4:0<-12,12,0.1>Mid Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mid_gain_tag {};

  void set (mid_gain_tag, float v)
  {
    // Original slider line: slider4:0<-12,12,0.1>Mid Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider4) {
      return;
    }
    slider4 = v;
    slider();
  }

  static constexpr auto get_parameter (mid_gain_tag)
  {
    // Original slider line: slider4:0<-12,12,0.1>Mid Gain (dB)
    return float_param ("dB", -12.0, 12.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider5_slider (float v)
  {
    // Original slider line: slider5:0<0,100,0.1>High Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct high_drive_tag {};

  void set (high_drive_tag, float v)
  {
    // Original slider line: slider5:0<0,100,0.1>High Drive (%)
    // Range: min:0.0, max:100.0, default: 0.0, step: 0.1
    if (v == slider5) {
      return;
    }
    slider5 = v;
    slider();
  }

  static constexpr auto get_parameter (high_drive_tag)
  {
    // Original slider line: slider5:0<0,100,0.1>High Drive (%)
    return float_param ("", 0.0, 100.0, 0.0, 0.3);
  }

#endif
#if 0
  void set_slider6_slider (float v)
  {
    // Original slider line: slider6:0<-12,12,0.1>High Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct high_gain_tag {};

  void set (high_gain_tag, float v)
  {
    // Original slider line: slider6:0<-12,12,0.1>High Gain (dB)
    // Range: min:-12.0, max:12.0, default: 0.0, step: 0.1
    if (v == slider6) {
      return;
    }
    slider6 = v;
    slider();
  }

  static constexpr auto get_parameter (high_gain_tag)
  {
    // Original slider line: slider6:0<-12,12,0.1>High Gain (dB)
    return float_param ("", -12.0, 12.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_slider7_slider (float v)
  {
    // Original slider line: slider7:240<60,500,1>Low-Mid Crossover (Hz)
    // Range: min:60.0, max:500.0, default: 240.0, step: 1.0
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct low_mid_freq_tag {};

  void set (low_mid_freq_tag, float v)
  {
    // Original slider line: slider7:240<60,500,1>Low-Mid Crossover (Hz)
    // Range: min:60.0, max:500.0, default: 240.0, step: 1.0
    v = midi_note_to_hz (v);
    if (v == slider7) {
      return;
    }
    slider7 = v;
    slider();
  }

  static constexpr auto get_parameter (low_mid_freq_tag)
  {
    // Original slider line: slider7:240<60,500,1>Low-Mid Crossover (Hz)
    return frequency_parameter (60.0, 500.0, 240.0);
  }

#endif
#if 0
  void set_slider8_slider (float v)
  {
    // Original slider line: slider8:2400<510,10000,10>Mid-High Crossover (Hz)
    // Range: min:510.0, max:10000.0, default: 2400.0, step: 10.0
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct mid_high_freq_tag {};

  void set (mid_high_freq_tag, float v)
  {
    // Original slider line: slider8:2400<510,10000,10>Mid-High Crossover (Hz)
    // Range: min:510.0, max:10000.0, default: 2400.0, step: 10.0
    v = midi_note_to_hz (v);
    if (v == slider8) {
      return;
    }
    slider8 = v;
    slider();
  }

  static constexpr auto get_parameter (mid_high_freq_tag)
  {
    // Original slider line: slider8:2400<510,10000,10>Mid-High Crossover (Hz)
    return frequency_parameter (510.0, 10000.0, 2400.0);
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    low_drive_tag,
    low_gain_tag,
    mid_drive_tag,
    mid_gain_tag,
    high_drive_tag,
    high_gain_tag,
    low_mid_freq_tag,
    mid_high_freq_tag>;
  //----------------------------------------------------------------------------
private:
  // global/stateful variables for section "init"
  double db2log;
  double halfpi;
  double halfpiscaled;
  double log2db;
  double pi;
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
    db2log       = 0;
    halfpi       = 0;
    halfpiscaled = 0;
    log2db       = 0;
    pi           = 0;
    slider1      = 0;
    slider2      = 0;
    slider3      = 0;
    slider4      = 0;
    slider5      = 0;
    slider6      = 0;
    slider7      = 0;
    slider8      = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "slider"
  double ah;
  double al;
  double driveh2;
  double drivel2;
  double drivem2;
  double mixhg1;
  double mixhgd;
  double mixlg1;
  double mixlgd;
  double mixmg1;
  double mixmgd;
  //----------------------------------------------------------------------------
  void init_slider_variables()
  {
    ah      = 0;
    al      = 0;
    driveh2 = 0;
    drivel2 = 0;
    drivem2 = 0;
    mixhg1  = 0;
    mixhgd  = 0;
    mixlg1  = 0;
    mixlgd  = 0;
    mixmg1  = 0;
    mixmgd  = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double lfh;
  double lfl;
  double rfh;
  double rfl;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    lfh = 0;
    lfl = 0;
    rfh = 0;
    rfl = 0;
  }
  //----------------------------------------------------------------------------
  plugin_context* plugcontext = nullptr;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_slider_variables();
    init_block_variables();
    slider1      = 0.0;
    slider2      = 0.0;
    slider3      = 0.0;
    slider4      = 0.0;
    slider5      = 0.0;
    slider6      = 0.0;
    slider7      = 240.0;
    slider8      = 2400.0;
    log2db       = 8.6858896380650365530225783783321;
    db2log       = 0.11512925464970228420089957273422;
    pi           = 3.1415926535;
    halfpi       = pi / 2.;
    halfpiscaled = halfpi * 1.41254;
    slider();
  }
  //----------------------------------------------------------------------------
private:
  void slider()
  {
    double driveh  = 0.;
    double driveh1 = 0.;
    double drivel  = 0.;
    double drivel1 = 0.;
    double drivem  = 0.;
    double drivem1 = 0.;
    double gainh   = 0.;
    double gainl   = 0.;
    double gainm   = 0.;
    double mixh    = 0.;
    double mixh1   = 0.;
    double mixh2   = 0.;
    double mixhg   = 0.;
    double mixl    = 0.;
    double mixl1   = 0.;
    double mixl2   = 0.;
    double mixlg   = 0.;
    double mixm    = 0.;
    double mixm1   = 0.;
    double mixm2   = 0.;
    double mixmg   = 0.;
    mixl           = slider1 / 100.;
    mixm           = slider3 / 100.;
    mixh           = slider5 / 100.;
    drivel         = mixl;
    drivem         = mixm;
    driveh         = mixh;
    drivel1        = 1. / (1. - (drivel / 2.));
    drivem1        = 1. / (1. - (drivem / 2.));
    driveh1        = 1. / (1. - (driveh / 2.));
    drivel2        = drivel / 2.;
    drivem2        = drivem / 2.;
    driveh2        = driveh / 2.;
    al             = std::min (slider7, jsfx_specialvar_get_srate())
      / jsfx_specialvar_get_srate();
    ah = std::max (
      std::min (slider8, jsfx_specialvar_get_srate())
        / jsfx_specialvar_get_srate(),
      al);
    mixl1  = 1. - mixl;
    mixm1  = 1. - mixm;
    mixh1  = 1. - mixh;
    mixl2  = mixl / 2.;
    mixm2  = mixm / 2.;
    mixh2  = mixh / 2.;
    gainl  = std::exp (slider2 * db2log);
    gainm  = std::exp (slider4 * db2log);
    gainh  = std::exp (slider6 * db2log);
    mixlg  = mixl * gainl;
    mixmg  = mixm * gainm;
    mixhg  = mixh * gainh;
    mixlgd = mixl * gainl * drivel1;
    mixmgd = mixm * gainm * drivem1;
    mixhgd = mixh * gainh * driveh1;
    mixlg1 = mixl1 * gainl;
    mixmg1 = mixm1 * gainm;
    mixhg1 = mixh1 * gainh;
    ;
  }
  //----------------------------------------------------------------------------
public:
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint block_samples)
  {
    double dry0   = 0.;
    double dry0_h = 0.;
    double dry0_l = 0.;
    double dry0_m = 0.;
    double dry1   = 0.;
    double dry1_h = 0.;
    double dry1_l = 0.;
    double dry1_m = 0.;
    double high_l = 0.;
    double high_r = 0.;
    double lf1h   = 0.;
    double lf1l   = 0.;
    double low_l  = 0.;
    double low_r  = 0.;
    double mid_l  = 0.;
    double mid_r  = 0.;
    double rf1h   = 0.;
    double rf1l   = 0.;
    double wet0   = 0.;
    double wet0_h = 0.;
    double wet0_l = 0.;
    double wet0_m = 0.;
    double wet1   = 0.;
    double wet1_h = 0.;
    double wet1_l = 0.;
    double wet1_m = 0.;

    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = chnls[0][$$i];
      auto& spl1 = chnls[1][$$i];
      dry0       = spl0;
      dry1       = spl1;
      lf1h       = lfh;
      lfh        = dry0 + lfh - ah * lf1h;
      high_l     = dry0 - lfh * ah;
      lf1l       = lfl;
      lfl        = dry0 + lfl - al * lf1l;
      low_l      = lfl * al;
      mid_l      = dry0 - low_l - high_l;
      rf1h       = rfh;
      rfh        = dry1 + rfh - ah * rf1h;
      high_r     = dry1 - rfh * ah;
      rf1l       = rfl;
      rfl        = dry1 + rfl - al * rf1l;
      low_r      = rfl * al;
      mid_r      = dry1 - low_r - high_r;
      wet0_l     = mixlgd * low_l * (1. - std::abs (low_l * drivel2));
      wet0_m     = mixmgd * mid_l * (1. - std::abs (mid_l * drivem2));
      wet0_h     = mixhgd * high_l * (1. - std::abs (high_l * driveh2));
      wet0       = (wet0_l + wet0_m + wet0_h);
      dry0_l     = low_l * mixlg1;
      dry0_m     = mid_l * mixmg1;
      dry0_h     = high_l * mixhg1;
      dry0       = (dry0_l + dry0_m + dry0_h);
      wet1_l     = mixlgd * low_r * (1. - std::abs (low_r * drivel2));
      wet1_m     = mixmgd * mid_r * (1. - std::abs (mid_r * drivem2));
      wet1_h     = mixhgd * high_r * (1. - std::abs (high_r * driveh2));
      wet1       = (wet1_l + wet1_m + wet1_h);
      dry1_l     = low_r * mixlg1;
      dry1_m     = mid_r * mixmg1;
      dry1_h     = high_r * mixhg1;
      dry1       = (dry1_l + dry1_m + dry1_h);
      spl0       = dry0 + wet0;
      spl1       = dry1 + wet1;
      ;
    }
  }
}; /* jsfx_process */

}} // namespace artv::sstillwell
