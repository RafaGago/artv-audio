#pragma once
// Generated by jsfx2cpp.py. To be manually corrected.
// includes for environment function calls

// REMINDER: to generate this I had to comment the UI imports.
// From commit 7ecb9004b77fcaa43b4f98bcd158cc2c968f5814

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

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

namespace artv { namespace geraint_luff {

class sandwitch_amp {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::distortion;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  // definitions for environment function calls
private:
  static double eel2_eq (double lhs, double rhs)
  {
    return (double) (std::abs (lhs - rhs) < 0.00001);
  }
  static bool eel2_ne (double lhs, double rhs) { return !eel2_eq (lhs, rhs); }
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
  void set_limit_db_slider (float v)
  {
    // Original slider line: slider1:limit_db=-12<-60,0,0.1>-Level (dB)
    // Range: min:-60.0, max:0.0, default: -12.0, step: 0.1
    if (v == limit_db) {
      return;
    }
    limit_db = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct limit_db_tag {};

  void set (limit_db_tag, float v)
  {
    limit_db = v;
  }

  static constexpr auto get_parameter (limit_db_tag)
  {
    // Original slider line: slider1:limit_db=-12<-60,0,0.1>-Level (dB)
    return float_param ("dB", -60.0, 0.0, -12.0, 0.1);
  }

#endif
#if 0
  void set_asymmetry_slider (float v)
  {
    // Original slider line: slider2:asymmetry=0<-1,1,0.01>-Asymmetry
    // Range: min:-1.0, max:1.0, default: 0.0, step: 0.01
    if (v == asymmetry) {
      return;
    }
    asymmetry = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct asymmetry_tag {};

  void set (asymmetry_tag, float v)
  {
    asymmetry = v;
  }

  static constexpr auto get_parameter (asymmetry_tag)
  {
    // Original slider line: slider2:asymmetry=0<-1,1,0.01>-Asymmetry
    return float_param ("", -1.0, 1.0, 0.0, 0.01);
  }

#endif
#if 0
  void set_output_db_slider (float v)
  {
    // Original slider line: slider3:output_db=0<-12,20,0.1>-Gain (dB)
    // Range: min:-12.0, max:20.0, default: 0.0, step: 0.1
    if (v == output_db) {
      return;
    }
    output_db = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct output_db_tag {};

  void set (output_db_tag, float v)
  {
    output_db = v;
  }

  static constexpr auto get_parameter (output_db_tag)
  {
    // Original slider line: slider3:output_db=0<-12,20,0.1>-Gain (dB)
    return float_param ("dB", -12.0, 20.0, 0.0, 0.1);
  }

#endif
#if 0
  void set_distortion_width_slider (float v)
  {
    // Original slider line: slider4:distortion_width=0.5<0.1,2,0.01>-Width
    // Range: min:0.1, max:2.0, default: 0.5, step: 0.01
    if (v == distortion_width) {
      return;
    }
    distortion_width = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct distortion_width_tag {};

  void set (distortion_width_tag, float v)
  {
    distortion_width = v;
  }

  static constexpr auto get_parameter (distortion_width_tag)
  {
    // Original slider line: slider4:distortion_width=0.5<0.1,2,0.01>-Width
    return float_param ("", 0.1, 2.0, 0.5, 0.01);
  }

#endif
#if 0
  void set_filter_freq_slider (float v)
  {
    // Original slider line: slider5:filter_freq=900<50,4000,1>-Filter freq
    // Range: min:50.0, max:4000.0, default: 900.0, step: 1.0
    if (v == filter_freq) {
      return;
    }
    filter_freq = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct filter_freq_tag {};

  void set (filter_freq_tag, float v)
  {
    v           = midi_note_to_hz (v);
    filter_freq = v;
  }

  static constexpr auto get_parameter (filter_freq_tag)
  {
    // Original slider line: slider5:filter_freq=900<50,4000,1>-Filter freq
    return frequency_parameter (50.0, 4000.0, 900.0);
  }

#endif
#if 0
  void set_filter_gain_slider (float v)
  {
    // Original slider line: slider6:filter_gain=20<0,50,0.1>-Filter gain
    // Range: min:0.0, max:50.0, default: 20.0, step: 0.1
    if (v == filter_gain) {
      return;
    }
    filter_gain = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct filter_gain_tag {};

  void set (filter_gain_tag, float v)
  {
    filter_gain = v;
  }

  static constexpr auto get_parameter (filter_gain_tag)
  {
    // Original slider line: slider6:filter_gain=20<0,50,0.1>-Filter gain
    return float_param ("", 0.0, 50.0, 20.0, 0.1);
  }

#endif
#if 0
  void set_bandwidth_octaves_slider (float v)
  {
    // Original slider line: slider7:bandwidth_octaves=2<0.5,4,0.1>-Filter
    // bandwidth (octaves) Range: min:0.5, max:4.0, default: 2.0, step: 0.1
    if (v == bandwidth_octaves) {
      return;
    }
    bandwidth_octaves = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct bandwidth_octaves_tag {};

  void set (bandwidth_octaves_tag, float v)
  {
    bandwidth_octaves = v;
  }

  static constexpr auto get_parameter (bandwidth_octaves_tag)
  {
    // Original slider line: slider7:bandwidth_octaves=2<0.5,4,0.1>-Filter
    // bandwidth (octaves)
    return float_param ("Octaves", 0.5, 4.0, 2.0, 0.1);
  }

#endif
#if 0
  void set_secondary_gain_db_slider (float v)
  {
    // Original slider line: slider8:secondary_gain_db=0<-40,40,0.1>-Secondary
    // gain Range: min:-40.0, max:40.0, default: 0.0, step: 0.1
    if (v == secondary_gain_db) {
      return;
    }
    secondary_gain_db = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct secondary_gain_db_tag {};

  void set (secondary_gain_db_tag, float v)
  {
    secondary_gain_db = v;
  }

  static constexpr auto get_parameter (secondary_gain_db_tag)
  {
    // Original slider line: slider8:secondary_gain_db=0<-40,40,0.1>-Secondary
    // gain
    return float_param ("dB", -40.0, 40.0, 0.0, 0.1);
  }

#endif
#if 0
#else
  // Snippet for parameter boilerplate in the authors framework....
  using parameters = mp_list<
    limit_db_tag,
    asymmetry_tag,
    output_db_tag,
    distortion_width_tag,
    filter_freq_tag,
    filter_gain_tag,
    bandwidth_octaves_tag,
    secondary_gain_db_tag>;
#endif
private:
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double asymmetry;
  double asymmetry$step;
  double asymmetry$value;
  double bandwidth_octaves;
  double bandwidth_octaves$step;
  double bandwidth_octaves$value;
  double bufferlength;
  double bufferpos;
  double distortion_width;
  double distortion_width$step;
  double distortion_width$value;
  double filter_a0;
  double filter_a1;
  double filter_a2;
  double filter_b0;
  double filter_b1;
  double filter_b2;
  double filter_freq;
  double filter_freq$step;
  double filter_freq$value;
  double filter_gain;
  double filter_gain$step;
  double filter_gain$value;
  double freemem;
  double g_scale;
  double gain;
  double gfx_ext_retina;
  double inputbuffer;
  double limit_db;
  double limit_db$step;
  double limit_db$value;
  double limit_reference;
  double offset;
  double offset_post;
  double offset_tanh;
  double output_db;
  double output_db$step;
  double output_db$value;
  double outputbuffer;
  double secondary_gain;
  double secondary_gain_db;
  double secondary_gain_db$step;
  double secondary_gain_db$value;
  double width;
  double width_inv;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    asymmetry               = 0;
    asymmetry$step          = 0;
    asymmetry$value         = 0;
    bandwidth_octaves       = 0;
    bandwidth_octaves$step  = 0;
    bandwidth_octaves$value = 0;
    bufferlength            = 0;
    bufferpos               = 0;
    distortion_width        = 0;
    distortion_width$step   = 0;
    distortion_width$value  = 0;
    filter_a0               = 0;
    filter_a1               = 0;
    filter_a2               = 0;
    filter_b0               = 0;
    filter_b1               = 0;
    filter_b2               = 0;
    filter_freq             = 0;
    filter_freq$step        = 0;
    filter_freq$value       = 0;
    filter_gain             = 0;
    filter_gain$step        = 0;
    filter_gain$value       = 0;
    freemem                 = 0;
    g_scale                 = 0;
    gain                    = 0;
    gfx_ext_retina          = 0;
    inputbuffer             = 0;
    limit_db                = 0;
    limit_db$step           = 0;
    limit_db$value          = 0;
    limit_reference         = 0;
    offset                  = 0;
    offset_post             = 0;
    offset_tanh             = 0;
    output_db               = 0;
    output_db$step          = 0;
    output_db$value         = 0;
    outputbuffer            = 0;
    secondary_gain          = 0;
    secondary_gain_db       = 0;
    secondary_gain_db$step  = 0;
    secondary_gain_db$value = 0;
    width                   = 0;
    width_inv               = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double filter_left$x1;
  double filter_left$x2;
  double filter_left$y1;
  double filter_left$y2;
  double filter_left_inv$x1;
  double filter_left_inv$x2;
  double filter_left_inv$y1;
  double filter_left_inv$y2;
  double filter_right$x1;
  double filter_right$x2;
  double filter_right$y1;
  double filter_right$y2;
  double filter_right_inv$x1;
  double filter_right_inv$x2;
  double filter_right_inv$y1;
  double filter_right_inv$y2;
  double needs_update;
  double spl2;
  double spl3;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    filter_left$x1      = 0;
    filter_left$x2      = 0;
    filter_left$y1      = 0;
    filter_left$y2      = 0;
    filter_left_inv$x1  = 0;
    filter_left_inv$x2  = 0;
    filter_left_inv$y1  = 0;
    filter_left_inv$y2  = 0;
    filter_right$x1     = 0;
    filter_right$x2     = 0;
    filter_right$y1     = 0;
    filter_right$y2     = 0;
    filter_right_inv$x1 = 0;
    filter_right_inv$x2 = 0;
    filter_right_inv$y1 = 0;
    filter_right_inv$y2 = 0;
    needs_update        = 0;
    spl2                = 0;
    spl3                = 0;
  }
  //----------------------------------------------------------------------------
  plugin_context* plugcontext = nullptr;
  uint            _last_block_samples;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_block_variables();
    _last_block_samples = pc.get_max_block_samples();

    limit_db          = -12.0;
    asymmetry         = 0.0;
    output_db         = 0.0;
    distortion_width  = 0.5;
    filter_freq       = 900.0;
    filter_gain       = 20.0;
    bandwidth_octaves = 2.0;
    secondary_gain_db = 0.0;
    gfx_ext_retina    = 1.;
#if 0
    freemem           = ui_setup (0.) /*TODO: unknown call */;
#else
    freemem = 0;
#endif

    bufferlength = std::floor (0.1 * jsfx_specialvar_get_srate());
    inputbuffer  = freemem;
    outputbuffer = freemem + bufferlength;
    freemem += bufferlength * 2.;
    heap_reset (freemem);
    bufferpos = 0.;
    ;
    limit_reference = 1.;
    offset          = 0.;
    offset_post     = 0.;
    gain            = 1.;
    secondary_gain  = 1.;
    width_inv       = 1.;
    width           = width_inv;
    filter_a0       = 1.;
    filter_a1       = 0.;
    filter_a2       = 0.;
    filter_b0       = 1.;
    filter_b1       = 0.;
    filter_b2       = 0.;
    ;
    g_scale = [&] {
      if (g_scale) {
        return g_scale;
      }
      else {
        return 1.;
      }
    }();
    init$smoother_init (limit_db$value, limit_db, limit_db$step);
    init$smoother_init (asymmetry$value, asymmetry, asymmetry$step);
    init$smoother_init (output_db$value, output_db, output_db$step);
    init$smoother_init (
      distortion_width$value, distortion_width, distortion_width$step);
    init$smoother_init (filter_freq$value, filter_freq, filter_freq$step);
    init$smoother_init (filter_gain$value, filter_gain, filter_gain$step);
    init$smoother_init (
      bandwidth_octaves$value, bandwidth_octaves, bandwidth_octaves$step);
    init$smoother_init (
      secondary_gain_db$value, secondary_gain_db, secondary_gain_db$step);
    init$update();
    ;
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint block_samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    double avg       = 0.;
    double diff      = 0.;
    double left      = 0.;
    double right     = 0.;
    double smoothing = 0.;

    _last_block_samples = block_samples;

    smoothing
      = init$smoother_block_1 (limit_db, limit_db$value, limit_db$step)
      + init$smoother_block_1 (asymmetry, asymmetry$value, asymmetry$step)
      + init$smoother_block_1 (output_db, output_db$value, output_db$step)
      + init$smoother_block_1 (
          distortion_width, distortion_width$value, distortion_width$step)
      + init$smoother_block_1 (filter_freq, filter_freq$value, filter_freq$step)
      + init$smoother_block_1 (filter_gain, filter_gain$value, filter_gain$step)
      + init$smoother_block_1 (
          bandwidth_octaves, bandwidth_octaves$value, bandwidth_octaves$step)
      + init$smoother_block_1 (
          secondary_gain_db, secondary_gain_db$value, secondary_gain_db$step);
    if (needs_update) {
      needs_update = 0.;
      init$update();
    };
    for (int $$i = 0, $$end = block_samples; $$i < $$end; ++$$i) {
      auto& spl0 = outs[0][$$i];
      auto& spl1 = outs[1][$$i];
      spl0       = ins[0][$$i];
      spl1       = ins[1][$$i];

      spl0 += 1e-15;
      spl1 += 1e-15;
      if (smoothing) {
        init$smoother_sample (limit_db$value, limit_db$step);
        init$smoother_sample (asymmetry$value, asymmetry$step);
        init$smoother_sample (output_db$value, output_db$step);
        init$smoother_sample (distortion_width$value, distortion_width$step);
        init$smoother_sample (filter_freq$value, filter_freq$step);
        init$smoother_sample (filter_gain$value, filter_gain$step);
        init$smoother_sample (bandwidth_octaves$value, bandwidth_octaves$step);
        init$smoother_sample (secondary_gain_db$value, secondary_gain_db$step);
        init$update();
      }
      heap (inputbuffer + bufferpos) = (spl0 + spl1) * 0.5;
      left                           = spl0 + spl2 * secondary_gain;
      right                          = spl1 + spl3 * secondary_gain;
      left                           = block$filter (
        left, filter_left$x1, filter_left$x2, filter_left$y1, filter_left$y2);
      right = block$filter (
        right,
        filter_right$x1,
        filter_right$x2,
        filter_right$y1,
        filter_right$y2);
      avg   = (left + right) * 0.5;
      diff  = (right - left) * 0.5;
      left  = avg - diff * width_inv;
      right = avg + diff * width_inv;
      left  = block$distort (left);
      right = block$distort (right);
      avg   = (left + right) * 0.5;
      diff  = (right - left) * 0.5;
      left  = avg - diff * width;
      right = avg + diff * width;
      left  = block$inv_filter (
        left,
        filter_left_inv$x1,
        filter_left_inv$x2,
        filter_left_inv$y1,
        filter_left_inv$y2);
      right = block$inv_filter (
        right,
        filter_right_inv$x1,
        filter_right_inv$x2,
        filter_right_inv$y1,
        filter_right_inv$y2);
      spl0                            = gain * (left - spl2 * secondary_gain);
      spl1                            = gain * (right - spl3 * secondary_gain);
      heap (outputbuffer + bufferpos) = (spl0 + spl1) * 0.5;
      bufferpos += 1.;
      if (bufferpos >= bufferlength) {
        bufferpos = 0.;
      };
    }
  }
  // functions for section "init"
private:
  //----------------------------------------------------------------------------
  double init$sinh (double x)
  {
    return (std::exp (x) - std::exp (-x)) * 0.5;
  }
  //----------------------------------------------------------------------------
  double init$smoother_block_1 (double& this$, double& $value, double& $step)
  {
    return [&] {
      if (eel2_ne (this$, $value)) {
        $step = (this$ - $value) / (double) _last_block_samples;
        return 1.;
      }
      else {
        ($step = 0.);
        return $step;
      }
    }();
  }
  //----------------------------------------------------------------------------
  double init$smoother_init (double& $value, double& this$, double& $step)
  {
    $value = this$;
    $step  = 0.;
    return 1.;
  }
  //----------------------------------------------------------------------------
  double init$smoother_sample (double& $value, double& $step)
  {
    $value += $step;
    return $value;
  }
  //----------------------------------------------------------------------------
  double init$smoother_value (double& $value)
  {
    return $value;
  }
  //----------------------------------------------------------------------------
  double init$tanh (double x)
  {
    return [&] {
      if (x >= 20.) {
        return 1.;
      }
      else {
        return [&] {
          if (x <= -20.) {
            return -(1.);
          }
          else {
            return ((std::exp (2. * x) - 1.) / (std::exp (2. * x) + 1.));
          }
        }();
      }
    }();
  }
  //----------------------------------------------------------------------------
  double init$update()
  {
    double cos_w0 = 0.;
    double a      = 0.;
    double alpha  = 0.;
    double w0     = 0.;
    limit_reference
      = std::pow (10., init$smoother_value (limit_db$value) / 20.);
    offset      = init$smoother_value (asymmetry$value);
    offset_tanh = init$tanh (offset);
    gain        = std::pow (10., init$smoother_value (output_db$value) / 20.);
    secondary_gain
      = std::pow (10., init$smoother_value (secondary_gain_db$value) / 20.);
    width     = init$smoother_value (distortion_width$value);
    width_inv = 1. / width;
    w0        = 2. * 3.141592653589793 * init$smoother_value (filter_freq$value)
      / jsfx_specialvar_get_srate();
    cos_w0 = std::cos (w0);
    alpha  = std::sin (w0)
      * init$sinh (
              std::log (2.) / 2. * init$smoother_value (bandwidth_octaves$value)
              * w0 / std::sin (w0));
    a         = std::pow (10., init$smoother_value (filter_gain$value) / 40.);
    filter_a0 = 1. + alpha / a;
    filter_a1 = -2. * cos_w0;
    filter_a2 = 1. - alpha / a;
    filter_b0 = 1. + alpha * a;
    filter_b1 = -2. * cos_w0;
    filter_b2 = 1. - alpha * a;
    return filter_b2;
  }
  // functions for section "block"

  //----------------------------------------------------------------------------
  double block$distort (double v)
  {
    return limit_reference
      * (init$tanh (v / limit_reference - offset) + offset_tanh);
  }
  //----------------------------------------------------------------------------
  double block$filter (
    double  x,
    double& $x1,
    double& $x2,
    double& $y1,
    double& $y2)
  {
    double y = 0.;
    y = (filter_b0 * x + filter_b1 * $x1 + filter_b2 * $x2 - filter_a1 * $y1
         - filter_a2 * $y2)
      / filter_a0;
    $y2 = $y1;
    $y1 = y;
    $x2 = $x1;
    $x1 = x;
    return y;
  }
  //----------------------------------------------------------------------------
  double block$inv_filter (
    double  x,
    double& $x1,
    double& $x2,
    double& $y1,
    double& $y2)
  {
    double y = 0.;
    y = (filter_a0 * x + filter_a1 * $x1 + filter_a2 * $x2 - filter_b1 * $y1
         - filter_b2 * $y2)
      / filter_b0;
    $x2 = $x1;
    $x1 = x;
    $y2 = $y1;
    $y1 = y;
    return y;
  }
}; /* jsfx_process */
}} // namespace artv::geraint_luff
