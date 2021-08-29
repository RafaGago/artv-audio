#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"

// Generated by jsfx2cpp.py. To be manually corrected.

// REMINDER: to generate this I had to comment the UI imports.
// From commit 7ecb9004b77fcaa43b4f98bcd158cc2c968f5814

// includes for environment function calls
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace geraint_luff {

class echo_cycles {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::delay;
  //----------------------------------------------------------------------------
private:
  // definitions for environment function calls
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
  void jsfx_memset (size_t idx, int val, size_t size)
  {
    std::memset (&heapmem[idx], val, size * sizeof heapmem[0]);
  }
  //----------------------------------------------------------------------------
  // stubs for JSFX special variables
  double jsfx_specialvar_get_beat_position()
  {
    return plugcontext->get_play_state().quarter_note_position;
  }

  double jsfx_specialvar_get_play_state()
  {
    return plugcontext->get_play_state().is_playing ? 1. : 0.;
  }

  double jsfx_specialvar_get_srate() { return plugcontext->get_sample_rate(); }

  double jsfx_specialvar_get_samplesblock() { return samples_block; }

  double jsfx_specialvar_get_tempo()
  {
    return plugcontext->get_play_state().bpm;
  }
  //----------------------------------------------------------------------------
public:
#if 0
  void set_input_width_slider (float v)
  {
    // Original slider line: slider1:input_width=0.5<0,1>-Input width
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    if (v == input_width) {
      return;
    }
    input_width = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct input_width_tag {};

  void set (input_width_tag, float v)
  {
    // Original slider line: slider1:input_width=0.5<0,1>-Input width
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    v /= 100.;
#if 0
    if (v == input_width) {
      return;
    }
#endif
    input_width = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (input_width_tag)
  {
    // Original slider line: slider1:input_width=0.5<0,1>-Input width
    return float_param ("%", 0.0, 100.0, 50, 0.1);
  }
#endif
#if 0
  void set_output_variation_slider (float v)
  {
    // Original slider line: slider2:output_variation=0.5<0,1>-Output variation
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    if (v == output_variation) {
      return;
    }
    output_variation = v;
#if 0
    slider();
#endif

  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct output_variation_tag {};

  void set (output_variation_tag, float v)
  {
    // Original slider line: slider2:output_variation=0.5<0,1>-Output variation
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    v /= 100.;
#if 0
    if (v == output_variation) {
      return;
    }
#endif
    output_variation = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (output_variation_tag)
  {
    // Original slider line: slider2:output_variation=0.5<0,1>-Output variation
    return float_param ("%", 0.0, 100.0, 50, 0.1);
  }
#endif
#if 0
  void set_delay_ms_slider (float v)
  {
    // Original slider line: slider3:delay_ms=0<0,1000>-Delay (ms)
    // Range: min:0.0, max:1000.0, default: 0.0, step: None
    if (v == delay_ms) {
      return;
    }
    delay_ms = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct delay_ms_tag {};

  void set (delay_ms_tag, float v)
  {
// Original slider line: slider3:delay_ms=0<0,1000>-Delay (ms)
// Range: min:0.0, max:1000.0, default: 0.0, step: None
#if 0
    if (v == delay_ms) {
      return;
    }
#endif
    delay_ms = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (delay_ms_tag)
  {
    // Original slider line: slider3:delay_ms=0<0,1000>-Delay (ms)
    return float_param ("ms", 0.0, 1000.0, 0.0, 0.1, 0.4);
  }

#endif
#if 0
  void set_delay_beats_slider (float v)
  {
    // Original slider line: slider4:delay_beats=0.75<0,4,0.25>-Delay (beats)
    // Range: min:0.0, max:4.0, default: 0.75, step: 0.25
    if (v == delay_beats) {
      return;
    }
    delay_beats = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct delay_beats_tag {};

  void set (delay_beats_tag, float v)
  {
// Original slider line: slider4:delay_beats=0.75<0,4,0.25>-Delay (beats)
// Range: min:0.0, max:4.0, default: 0.75, step: 0.25
#if 0
    if (v == delay_beats) {
      return;
    }
#endif
    delay_beats = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (delay_beats_tag)
  {
    // Original slider line: slider4:delay_beats=0.75<0,4,0.25>-Delay (beats)
    return float_param ("beats", 0.0, 4.0, 0.75, 0.25);
  }
#endif
#if 0
  void set_feedback_ratio_slider (float v)
  {
    // Original slider line: slider5:feedback_ratio=0.5<0,0.99>-Feedback
    // Range: min:0.0, max:0.99, default: 0.5, step: None
    if (v == feedback_ratio) {
      return;
    }
    feedback_ratio = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct feedback_ratio_tag {};

  void set (feedback_ratio_tag, float v)
  {
    // Original slider line: slider5:feedback_ratio=0.5<0,0.99>-Feedback
    // Range: min:0.0, max:0.99, default: 0.5, step: None
    v /= 101.01010101;
#if 0
    if (v == feedback_ratio) {
      return;
    }
#endif
    feedback_ratio = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (feedback_ratio_tag)
  {
    // Original slider line: slider5:feedback_ratio=0.5<0,0.99>-Feedback
    return float_param ("%", 0.0, 100, 50, 0.01);
  }
#endif
#if 0
  void set_feedback_rotation_slider (float v)
  {
    // Original slider line:
    // slider6:feedback_rotation=2.323<0,6.283185307179586>-Feedback rotation
    // Range: min:0.0, max:6.283185307179586, default: 2.323, step: None
    if (v == feedback_rotation) {
      return;
    }
    feedback_rotation = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct feedback_rotation_tag {};

  void set (feedback_rotation_tag, float v)
  {
    v *= (6.283185307179586 / 360.);
// Original slider line:
// slider6:feedback_rotation=2.323<0,6.283185307179586>-Feedback rotation
// Range: min:0.0, max:6.283185307179586, default: 2.323, step: None
#if 0
    if (v == feedback_rotation) {
      return;
    }
#endif
    feedback_rotation = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (feedback_rotation_tag)
  {
    // Original slider line:
    // slider6:feedback_rotation=2.323<0,6.283185307179586>-Feedback rotation
    return float_param ("deg", 0.0, 360, 90, 0.1);
  }

#endif
#if 0
  void set_input_rotation_initial_slider (float v)
  {
    // Original slider line:
    // slider7:input_rotation_initial=2.776<0,6.283185307179586>-Input rotation
    // Range: min:0.0, max:6.283185307179586, default: 2.776, step: None
    if (v == input_rotation_initial) {
      return;
    }
    input_rotation_initial = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct input_rotation_initial_tag {};

  void set (input_rotation_initial_tag, float v)
  {
    // Original slider line:
    // slider7:input_rotation_initial=2.776<0,6.283185307179586>-Input rotation
    // Range: min:0.0, max:6.283185307179586, default: 2.776, step: None
    v *= (6.283185307179586 / 360.);
#if 0
    if (v == input_rotation_initial) {
      return;
    }
#endif
    input_rotation_initial = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (input_rotation_initial_tag)
  {
    // Original slider line:
    // slider7:input_rotation_initial=2.776<0,6.283185307179586>-Input rotation
    return float_param ("deg", 0.0, 360, 2., 0.1);
  }
#endif
#if 0
  void set_output_dry_slider (float v)
  {
    // Original slider line: slider8:output_dry=1<0,1>-Output Dry
    // Range: min:0.0, max:1.0, default: 1.0, step: None
    if (v == output_dry) {
      return;
    }
    output_dry = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct output_dry_tag {};

  void set (output_dry_tag, float v)
  {
    // Original slider line: slider8:output_dry=1<0,1>-Output Dry
    // Range: min:0.0, max:1.0, default: 1.0, step: None
    if (v > -50.) {
      v = db_to_gain (v);
    }
    else {
      v = 0.;
    }

#if 0
    if (v == output_dry) {
      return;
    }
#endif
    output_dry = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (output_dry_tag)
  {
    // Original slider line: slider8:output_dry=1<0,1>-Output Dry
    return float_param ("dB", -50, 0.0, -50.0, 0.001, 1.3);
  }

#endif
#if 0
  void set_output_wet_slider (float v)
  {
    // Original slider line: slider9:output_wet=0.5<0,1>-Output Wet
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    if (v == output_wet) {
      return;
    }
    output_wet = v;
#if 0
    slider();
#endif

  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct output_wet_tag {};

  void set (output_wet_tag, float v)
  {
    // Original slider line: slider9:output_wet=0.5<0,1>-Output Wet
    // Range: min:0.0, max:1.0, default: 0.5, step: None
    v          = db_to_gain (v);
#if 0
    if (v == output_wet) {
      return;
    }
#endif
    output_wet = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (output_wet_tag)
  {
    // Original slider line: slider9:output_wet=0.5<0,1>-Output Wet
    return float_param ("dB", -50, 0.0, -3., 0.001, 1.3);
  }

#endif
#if 0
  void set_filter_freq_slider (float v)
  {
    // Original slider line: slider10:filter_freq=2000<50,18000>-Filter Hz
    // Range: min:50.0, max:18000.0, default: 2000.0, step: None
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
// Original slider line: slider10:filter_freq=2000<50 ,18000>-Filter Hz
// Range: min:50.0, max:18000.0, default: 2000.0, step: None
#if 0
    if (v == filter_freq) {
      return;
    }
#endif
    v           = midi_note_to_hz (v);
    filter_freq = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (filter_freq_tag)
  {
    // Original slider line: slider10:filter_freq=2000<50,18000>-Filter Hz
    return frequency_parameter (50.0, 18000.0, 2000.0);
  }

#endif
#if 0
  void set_filter_db_slider (float v)
  {
    // Original slider line: slider11:filter_db=2<0,24>-Filter dB (skirt)
    // Range: min:0.0, max:24.0, default: 2.0, step: None
    if (v == filter_db) {
      return;
    }
    filter_db = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct filter_db_tag {};

  void set (filter_db_tag, float v)
  {
// Original slider line: slider11:filter_db=2<0,24>-Filter dB (skirt)
// Range: min:0.0, max:24.0, default: 2.0, step: None
#if 0
    if (v == filter_db) {
      return;
    }
#endif
    filter_db = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (filter_db_tag)
  {
    // Original slider line: slider11:filter_db=2<0,24>-Filter dB (skirt)
    return float_param ("dB", 0.0, 24.0, 2.0, 0.1);
  }

#endif
#if 0
  void set_filter_bandwidth_slider (float v)
  {
    // Original slider line: slider12:filter_bandwidth=2<0.1,5>-Filter bandwidth
    // Range: min:0.1, max:5.0, default: 2.0, step: None
    if (v == filter_bandwidth) {
      return;
    }
    filter_bandwidth = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct filter_bandwidth_tag {};

  void set (filter_bandwidth_tag, float v)
  {
// Original slider line: slider12:filter_bandwidth=2<0.1,5>-Filter bandwidth
// Range: min:0.1, max:5.0, default: 2.0, step: None
#if 0
    if (v == filter_bandwidth) {
      return;
    }
#endif
    filter_bandwidth = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (filter_bandwidth_tag)
  {
    // Original slider line: slider12:filter_bandwidth=2<0.1,5>-Filter bandwidth
    return float_param ("", 0.1, 5.0, 2.0, 0.05);
  }

#endif
#if 0
  void set_rotation_mode_slider (float v)
  {
    // Original slider line: slider13:rotation_mode=0<0,1,1{fixed,LFO}>-Rotation
    // mode Range: min:0.0, max:1.0, default: 0.0, step: 1.0
    if (v == rotation_mode) {
      return;
    }
    rotation_mode = v;
    slider();
  }
#else
  // Snippet for parameter boilerplate in the authors framework....
  struct rotation_mode_tag {};

  void set (rotation_mode_tag, float v)
  {
// Original slider line: slider13:rotation_mode=0<0,1,1{fixed,LFO}>-Rotation
// mode Range: min:0.0, max:1.0, default: 0.0, step: 1.0
#if 0
    if (v == rotation_mode) {
      return;
    }
#endif
    rotation_mode = v;
#if 0
    slider();
#endif
  }

  static constexpr auto get_parameter (rotation_mode_tag)
  {
    // Original slider line: slider13:rotation_mode=0<0,1,1{fixed,LFO}>-Rotation
    // mode
    return choice_param (0, make_cstr_array ("Fixed", "LFO"));
  }

#endif
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    input_width_tag,
    output_variation_tag,
    delay_ms_tag,
    delay_beats_tag,
    feedback_ratio_tag,
    feedback_rotation_tag,
    input_rotation_initial_tag,
    output_dry_tag,
    output_wet_tag,
    filter_freq_tag,
    filter_db_tag,
    filter_bandwidth_tag,
    rotation_mode_tag>;

private:
  //----------------------------------------------------------------------------
  // global/stateful variables for section "init"
  double buffer_length;
  double cos_w0;
  double delay_beats;
  double delay_buffer_left;
  double delay_buffer_mono;
  double delay_buffer_right;
  double delay_gain;
  double delay_gain_step;
  double delay_ms;
  double feedback_ratio;
  double feedback_ratio$this$step;
  double feedback_ratio$this$value;
  double feedback_rotation;
  double feedback_rotation$this$step;
  double feedback_rotation$this$value;
  double filter_bandwidth;
  double filter_bandwidth$this$step;
  double filter_bandwidth$this$value;
  double filter_db;
  double filter_db$this$step;
  double filter_db$this$value;
  double filter_freq;
  double filter_freq$this$step;
  double filter_freq$this$value;
  double filter_gain;
  double filter_left$this$a1;
  double filter_left$this$a2;
  double filter_left$this$b0;
  double filter_left$this$b1;
  double filter_left$this$b2;
  double filter_left$this$x1;
  double filter_left$this$x2;
  double filter_left$this$y1;
  double filter_left$this$y2;
  double filter_mono$this$a1;
  double filter_mono$this$a2;
  double filter_mono$this$b0;
  double filter_mono$this$b1;
  double filter_mono$this$b2;
  double filter_mono$this$x1;
  double filter_mono$this$x2;
  double filter_mono$this$y1;
  double filter_mono$this$y2;
  double filter_right$this$a1;
  double filter_right$this$a2;
  double filter_right$this$b0;
  double filter_right$this$b1;
  double filter_right$this$b2;
  double filter_right$this$x1;
  double filter_right$this$x2;
  double filter_right$this$y1;
  double filter_right$this$y2;
  double freemem;
  double gfx_ext_retina;
  double input_cos;
  double input_rot_cos;
  double input_rot_sin;
  double input_rotation;
  double input_rotation$this$step;
  double input_rotation$this$value;
  double input_rotation_initial;
  double input_sin;
  double input_width;
  double max_delay;
  double max_delay_samples;
  double output_dry;
  double output_dry$this$step;
  double output_dry$this$value;
  double output_variation;
  double output_wet;
  double output_wet$this$step;
  double output_wet$this$value;
  double rotation_mode;
  double rotation_offset;
  double rotation_speed;
  double vis_sample;
  //----------------------------------------------------------------------------
  void init_init_variables()
  {
    buffer_length                = 0;
    cos_w0                       = 0;
    delay_beats                  = 0;
    delay_buffer_left            = 0;
    delay_buffer_mono            = 0;
    delay_buffer_right           = 0;
    delay_gain                   = 0;
    delay_gain_step              = 0;
    delay_ms                     = 0;
    feedback_ratio               = 0;
    feedback_ratio$this$step     = 0;
    feedback_ratio$this$value    = 0;
    feedback_rotation            = 0;
    feedback_rotation$this$step  = 0;
    feedback_rotation$this$value = 0;
    filter_bandwidth             = 0;
    filter_bandwidth$this$step   = 0;
    filter_bandwidth$this$value  = 0;
    filter_db                    = 0;
    filter_db$this$step          = 0;
    filter_db$this$value         = 0;
    filter_freq                  = 0;
    filter_freq$this$step        = 0;
    filter_freq$this$value       = 0;
    filter_gain                  = 0;
    filter_left$this$a1          = 0;
    filter_left$this$a2          = 0;
    filter_left$this$b0          = 0;
    filter_left$this$b1          = 0;
    filter_left$this$b2          = 0;
    filter_left$this$x1          = 0;
    filter_left$this$x2          = 0;
    filter_left$this$y1          = 0;
    filter_left$this$y2          = 0;
    filter_mono$this$a1          = 0;
    filter_mono$this$a2          = 0;
    filter_mono$this$b0          = 0;
    filter_mono$this$b1          = 0;
    filter_mono$this$b2          = 0;
    filter_mono$this$x1          = 0;
    filter_mono$this$x2          = 0;
    filter_mono$this$y1          = 0;
    filter_mono$this$y2          = 0;
    filter_right$this$a1         = 0;
    filter_right$this$a2         = 0;
    filter_right$this$b0         = 0;
    filter_right$this$b1         = 0;
    filter_right$this$b2         = 0;
    filter_right$this$x1         = 0;
    filter_right$this$x2         = 0;
    filter_right$this$y1         = 0;
    filter_right$this$y2         = 0;
    freemem                      = 0;
    gfx_ext_retina               = 0;
    input_cos                    = 0;
    input_rot_cos                = 0;
    input_rot_sin                = 0;
    input_rotation               = 0;
    input_rotation$this$step     = 0;
    input_rotation$this$value    = 0;
    input_rotation_initial       = 0;
    input_sin                    = 0;
    input_width                  = 0;
    max_delay                    = 0;
    max_delay_samples            = 0;
    output_dry                   = 0;
    output_dry$this$step         = 0;
    output_dry$this$value        = 0;
    output_variation             = 0;
    output_wet                   = 0;
    output_wet$this$step         = 0;
    output_wet$this$value        = 0;
    rotation_mode                = 0;
    rotation_offset              = 0;
    rotation_speed               = 0;
    vis_sample                   = 0;
  }
  //----------------------------------------------------------------------------
  // global/stateful variables for section "block"
  double buffer_read_index;
  double delay_samples;
  double delay_seconds;
  double feedback_cos;
  double feedback_sin;
  double input_ll;
  double input_lr;
  double input_mono;
  double input_rl;
  double input_rr;
  double output_off_channel;
  double output_on_channel;
  //----------------------------------------------------------------------------
  void init_block_variables()
  {
    buffer_read_index  = 0;
    delay_samples      = 0;
    delay_seconds      = 0;
    feedback_cos       = 0;
    feedback_sin       = 0;
    input_ll           = 0;
    input_lr           = 0;
    input_mono         = 0;
    input_rl           = 0;
    input_rr           = 0;
    output_off_channel = 0;
    output_on_channel  = 0;
  }
  plugin_context* plugcontext = nullptr;
  //----------------------------------------------------------------------------
public:
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    init_init_variables();
    init_block_variables();
#if 0
    input_width            = 0.5;
    output_variation       = 0.5;
    delay_ms               = 0.0;
    delay_beats            = 0.75;
    feedback_ratio         = 0.5;
    feedback_rotation      = 2.323;
    input_rotation_initial = 2.776;
    output_dry             = 1.0;
    output_wet             = 0.5;
    filter_freq            = 2000.0;
    filter_db              = 2.0;
    filter_bandwidth       = 2.0;
    rotation_mode          = 0.0;
#else
    // I changed the defaults and ranges
    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
#endif

    max_delay         = 8.;
    max_delay_samples = jsfx_specialvar_get_srate() * max_delay + 100.;
    gfx_ext_retina    = 1.;
    freemem           = 0.;
#if 0
    freemem                = ui_setup (freemem) /*TODO: unknown call */;
#endif
    buffer_length = max_delay_samples;
    (delay_buffer_mono = freemem);
    freemem = delay_buffer_mono + buffer_length;
    (delay_buffer_left = freemem);
    freemem = delay_buffer_left + buffer_length;
    (delay_buffer_right = freemem);
    freemem = delay_buffer_right + buffer_length;
    heap_reset (freemem);
    jsfx_memset (delay_buffer_mono, 0., buffer_length);
    jsfx_memset (delay_buffer_left, 0., buffer_length);
    jsfx_memset (delay_buffer_right, 0., buffer_length);
    vis_sample      = 0.;
    delay_gain      = 0.;
    delay_gain_step = 0.;
    init$filter_init (
      0.,
      filter_mono$this$y2,
      filter_mono$this$y1,
      filter_mono$this$x2,
      filter_mono$this$x1,
      filter_mono$this$b0,
      filter_mono$this$a2,
      filter_mono$this$a1,
      filter_mono$this$b2,
      filter_mono$this$b1);
    init$filter_init (
      0.,
      filter_left$this$y2,
      filter_left$this$y1,
      filter_left$this$x2,
      filter_left$this$x1,
      filter_left$this$b0,
      filter_left$this$a2,
      filter_left$this$a1,
      filter_left$this$b2,
      filter_left$this$b1);
    init$filter_init (
      0.,
      filter_right$this$y2,
      filter_right$this$y1,
      filter_right$this$x2,
      filter_right$this$x1,
      filter_right$this$b0,
      filter_right$this$a2,
      filter_right$this$a1,
      filter_right$this$b2,
      filter_right$this$b1);
    init$smoother_init (
      output_dry$this$value, output_dry, output_dry$this$step);
    init$smoother_init (
      output_wet$this$value, output_wet, output_wet$this$step);
    init$smoother_init (
      feedback_ratio$this$value, feedback_ratio, feedback_ratio$this$step);
    init$smoother_init (
      feedback_rotation$this$value,
      feedback_rotation,
      feedback_rotation$this$step);
    init$smoother_init (
      input_rotation$this$value, input_rotation, input_rotation$this$step);
    init$smoother_init (
      filter_freq$this$value, filter_freq, filter_freq$this$step);
    init$smoother_init (filter_db$this$value, filter_db, filter_db$this$step);
    init$smoother_init (
      filter_bandwidth$this$value,
      filter_bandwidth,
      filter_bandwidth$this$step);
    input_cos      = 1.;
    input_sin      = 0.;
    rotation_speed = 0.;
    ;
    init$recalculate();
    init$recalculate_input_rotation();
    rotation_offset = [&] {
      if (eel2_eq (rotation_mode, 1.)) {
        return jsfx_specialvar_get_beat_position() * rotation_speed;
      }
      else {
        return 0.;
      }
    }();
    ;
  }
  //----------------------------------------------------------------------------
  uint samples_block = 512;
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    samples_block               = samples;
    double buffer_write_index   = 0.;
    double delay_samples_target = 0.;
    double delayed_left         = 0.;
    double delayed_mono         = 0.;
    double delayed_right        = 0.;
    double left                 = 0.;
    double mono                 = 0.;
    double right                = 0.;
    double smoothing            = 0.;

    (rotation_offset = 0.);
    (rotation_offset += rotation_speed
       * (jsfx_specialvar_get_tempo() / 60. * jsfx_specialvar_get_samplesblock()
          / jsfx_specialvar_get_srate()));
    if (eel2_eq (rotation_mode, 1.)) {
      if (
        eel2_eq (jsfx_specialvar_get_play_state(), 1.)
        || eel2_eq (jsfx_specialvar_get_play_state(), 5.)) {
        (rotation_offset
         = jsfx_specialvar_get_beat_position() * rotation_speed);
      }
      else {
        rotation_offset;
      }
    }
    else {
      rotation_offset;
    }
    rotation_offset -= std::floor (rotation_offset);
    input_rotation
      = input_rotation_initial + 6.283185307179586 * rotation_offset;
    smoothing
      = init$smoother_block_1 (
          output_dry, output_dry$this$value, output_dry$this$step)
      + init$smoother_block_1 (
          output_wet, output_wet$this$value, output_wet$this$step)
      + init$smoother_block_1 (
          feedback_ratio, feedback_ratio$this$value, feedback_ratio$this$step)
      + init$smoother_block_2 (
          6.283185307179586,
          feedback_rotation,
          feedback_rotation$this$value,
          feedback_rotation$this$step)
      + init$smoother_block_2 (
          6.283185307179586,
          input_rotation,
          input_rotation$this$value,
          input_rotation$this$step)
      + init$smoother_block_1 (
          filter_freq, filter_freq$this$value, filter_freq$this$step)
      + init$smoother_block_1 (
          filter_db, filter_db$this$value, filter_db$this$step)
      + init$smoother_block_1 (
          filter_bandwidth,
          filter_bandwidth$this$value,
          filter_bandwidth$this$step);
    init$recalculate();
    init$recalculate_input_rotation();
    delay_seconds
      = delay_ms * 0.001 + delay_beats * 60. / jsfx_specialvar_get_tempo();
    delay_samples_target = std::max (
      1., std::floor (delay_seconds * jsfx_specialvar_get_srate() + 0.5));
    (delay_gain_step = (1. - delay_gain) / jsfx_specialvar_get_samplesblock());
    if (eel2_ne (delay_samples_target, delay_samples)) {
      if (delay_gain <= 0.01) {
        delay_samples = delay_samples_target;
        delay_gain_step
          = (1. - delay_gain) / jsfx_specialvar_get_samplesblock();
      }
      else {
        (delay_gain_step = -delay_gain / jsfx_specialvar_get_samplesblock());
      }
    }
    else {
      delay_gain_step;
    };
    for (int $$i = 0, $$end = samples; $$i < $$end; ++$$i) {
      T& spl0 = chnls[0][$$i];
      T& spl1 = chnls[1][$$i];
      if (smoothing) {
        init$smoother_sample (output_dry$this$value, output_dry$this$step);
        init$smoother_sample (output_wet$this$value, output_wet$this$step);
        init$smoother_sample (
          feedback_ratio$this$value, feedback_ratio$this$step);
        init$smoother_sample (
          feedback_rotation$this$value, feedback_rotation$this$step);
        init$smoother_sample (
          input_rotation$this$value, input_rotation$this$step);
        init$smoother_sample (filter_freq$this$value, filter_freq$this$step);
        init$smoother_sample (filter_db$this$value, filter_db$this$step);
        init$smoother_sample (
          filter_bandwidth$this$value, filter_bandwidth$this$step);
        init$recalculate();
        init$recalculate_input_rotation();
      }
      else {
        if (eel2_eq (rotation_mode, 1.)) {
          init$smoother_sample (
            input_rotation$this$value, input_rotation$this$step);
          init$recalculate_input_rotation();
        }
      }
      vis_sample += 1.;
      delay_gain += delay_gain_step;
      mono         = (spl0 + spl1) * input_mono;
      left         = input_ll * spl0 + input_lr * spl1;
      right        = input_rl * spl0 + input_rr * spl1;
      delayed_mono = init$filter_sample (
        heap (delay_buffer_mono + buffer_read_index),
        filter_mono$this$b0,
        filter_mono$this$b1,
        filter_mono$this$x1,
        filter_mono$this$b2,
        filter_mono$this$x2,
        filter_mono$this$a1,
        filter_mono$this$y1,
        filter_mono$this$a2,
        filter_mono$this$y2);
      delayed_left = init$filter_sample (
        heap (delay_buffer_left + buffer_read_index),
        filter_left$this$b0,
        filter_left$this$b1,
        filter_left$this$x1,
        filter_left$this$b2,
        filter_left$this$x2,
        filter_left$this$a1,
        filter_left$this$y1,
        filter_left$this$a2,
        filter_left$this$y2);
      delayed_right = init$filter_sample (
        heap (delay_buffer_right + buffer_read_index),
        filter_right$this$b0,
        filter_right$this$b1,
        filter_right$this$x1,
        filter_right$this$b2,
        filter_right$this$x2,
        filter_right$this$a1,
        filter_right$this$y1,
        filter_right$this$a2,
        filter_right$this$y2);
      mono += init$smoother_value (feedback_ratio$this$value) * delayed_mono;
      left += feedback_cos * delayed_left - feedback_sin * delayed_right;
      right += feedback_sin * delayed_left + feedback_cos * delayed_right;
      buffer_write_index = buffer_read_index + delay_samples;
      while (buffer_write_index >= buffer_length) {
        buffer_write_index -= buffer_length;
      }
      heap (delay_buffer_mono + buffer_write_index)  = mono * delay_gain;
      heap (delay_buffer_left + buffer_write_index)  = left * delay_gain;
      heap (delay_buffer_right + buffer_write_index) = right * delay_gain;
      spl0 = init$smoother_value (output_dry$this$value) * spl0
        + init$smoother_value (output_wet$this$value)
          * (delayed_mono + output_on_channel * delayed_left + output_off_channel * delayed_right);
      spl1 = init$smoother_value (output_dry$this$value) * spl1
        + init$smoother_value (output_wet$this$value)
          * (delayed_mono + output_on_channel * delayed_right + output_off_channel * delayed_left);
      buffer_read_index += 1.;
      if (buffer_read_index >= buffer_length) {
        (buffer_read_index = 0.);
      };
    }
  }

private:
  // functions for section "init"
  //----------------------------------------------------------------------------
  double init$filter_bandpass_2 (
    double  normfreq,
    double  bw,
    double& this$a1,
    double& this$a2,
    double& this$b0,
    double& this$b1,
    double& this$b2)
  {
    double inv_a0 = 0.;
    double alpha  = 0.;
    double b1     = 0.;
    double w0     = 0.;
    double sin_w0 = 0.;
    double q      = 0.;
    w0
      = 2. * 3.141592653589793 * std::max (0.000001, std::min (0.49, normfreq));
    sin_w0 = std::sin (w0);
    cos_w0 = std::cos (w0);
    (q = -bw);
    if (bw > 0.) {
      (q = 0.5 / init$sinh (std::log (2.) * 0.5 * bw * w0 / sin_w0));
    }
    else {
      q;
    }
    alpha   = sin_w0 * 0.5 / q;
    inv_a0  = 1. / (1. + alpha);
    this$a1 = -2. * cos_w0 * inv_a0;
    this$a2 = (1. - alpha) * inv_a0;
    this$b0 = alpha * inv_a0;
    this$b1 = 0.;
    this$b2 = -this$b0;
    return this$b2;
  }
  //----------------------------------------------------------------------------
  double init$filter_bandpass_3 (
    double  normfreq,
    double  bw,
    double  gain,
    double& this$b0,
    double& this$b1,
    double& this$a1,
    double& this$b2,
    double& this$a2)
  {
    double inv_a0 = 0.;
    double alpha  = 0.;
    double b1     = 0.;
    double w0     = 0.;
    double sin_w0 = 0.;
    double q      = 0.;
    init$filter_bandpass_2 (
      normfreq, bw, this$a1, this$a2, this$b0, this$b1, this$b2);
    this$b0 += (1. - this$b0) * gain;
    this$b1 += (this$a1 - this$b1) * gain;
    this$b2 += (this$a2 - this$b2) * gain;
    return this$b2;
  }
  //----------------------------------------------------------------------------
  double init$filter_init (
    double  steady_value,
    double& this$y2,
    double& this$y1,
    double& this$x2,
    double& this$x1,
    double& this$b0,
    double& this$a2,
    double& this$a1,
    double& this$b2,
    double& this$b1)
  {
    this$y2 = steady_value;
    this$y1 = this$y2;
    this$x2 = this$y1;
    this$x1 = this$x2;
    this$b0 = 1.;
    this$a2 = 0.;
    this$a1 = this$a2;
    this$b2 = this$a1;
    this$b1 = this$b2;
    return this$b1;
  }
  //----------------------------------------------------------------------------
  double init$filter_sample (
    double  x0,
    double& this$b0,
    double& this$b1,
    double& this$x1,
    double& this$b2,
    double& this$x2,
    double& this$a1,
    double& this$y1,
    double& this$a2,
    double& this$y2)
  {
    double y0 = 0.;
    y0        = this$b0 * x0 + this$b1 * this$x1 + this$b2 * this$x2
      - this$a1 * this$y1 - this$a2 * this$y2;
    this$x2 = this$x1;
    this$x1 = x0;
    this$y2 = this$y1;
    this$y1 = y0;
    return this$y1;
  }
  //----------------------------------------------------------------------------
  double init$recalculate()
  {
    input_mono = std::sqrt (0.5);
    input_cos  = std::cos ((input_width + 1.) * 3.141592653589793 / 4.);
    input_sin  = std::sin ((input_width + 1.) * 3.141592653589793 / 4.);
    output_on_channel
      = output_variation + (1. - output_variation) / std::sqrt (2.);
    output_off_channel = (output_variation - 1.) / std::sqrt (2.);
    delay_seconds
      = delay_ms * 0.001 + delay_beats * 60. / jsfx_specialvar_get_tempo();
    rotation_speed = - [&]{ if ( feedback_rotation > 3.141592653589793 ) { return feedback_rotation - 2. * 3.141592653589793 ; } else { return feedback_rotation ; } }()  / ( std::sqrt ( 2. ) + 2. * delay_beats ) ;
    filter_gain
      = std::pow (10., -init$smoother_value (filter_db$this$value) / 20.);
    feedback_cos = init$smoother_value (feedback_ratio$this$value)
      * std::cos (init$smoother_value (feedback_rotation$this$value));
    feedback_sin = init$smoother_value (feedback_ratio$this$value)
      * std::sin (init$smoother_value (feedback_rotation$this$value));
    init$filter_bandpass_3 (
      init$smoother_value (filter_freq$this$value)
        / jsfx_specialvar_get_srate(),
      init$smoother_value (filter_bandwidth$this$value),
      filter_gain,
      filter_mono$this$b0,
      filter_mono$this$b1,
      filter_mono$this$a1,
      filter_mono$this$b2,
      filter_mono$this$a2);
    init$filter_bandpass_3 (
      init$smoother_value (filter_freq$this$value)
        / jsfx_specialvar_get_srate(),
      init$smoother_value (filter_bandwidth$this$value),
      filter_gain,
      filter_left$this$b0,
      filter_left$this$b1,
      filter_left$this$a1,
      filter_left$this$b2,
      filter_left$this$a2);
    return init$filter_bandpass_3 (
      init$smoother_value (filter_freq$this$value)
        / jsfx_specialvar_get_srate(),
      init$smoother_value (filter_bandwidth$this$value),
      filter_gain,
      filter_right$this$b0,
      filter_right$this$b1,
      filter_right$this$a1,
      filter_right$this$b2,
      filter_right$this$a2);
  }
  //----------------------------------------------------------------------------
  double init$recalculate_input_rotation()
  {
    input_rot_cos = std::cos (init$smoother_value (input_rotation$this$value));
    input_rot_sin = std::sin (init$smoother_value (input_rotation$this$value));
    input_ll      = input_rot_cos * input_sin - input_rot_sin * input_cos;
    input_lr      = input_rot_cos * input_cos - input_rot_sin * input_sin;
    input_rl      = input_rot_sin * input_sin + input_rot_cos * input_cos;
    input_rr      = input_rot_sin * input_cos + input_rot_cos * input_sin;
    return input_rr;
  }
  //----------------------------------------------------------------------------
  double init$sinh (double x)
  {
    double e_x = 0.;
    return [&] {
      if (x > 20.) {
        return (std::exp (x) * 0.5);
      }
      else {
        return [&] {
          if (x < -20.) {
            return (-0.5 * std::exp (-x));
          }
          else {
            e_x = std::exp (x);
            return (e_x - 1. / e_x) * 0.5;
          }
        }();
      }
    }();
  }
  //----------------------------------------------------------------------------
  double init$smoother_block_1 (
    double& this$,
    double& this$value,
    double& this$step)
  {
    return [&] {
      if (eel2_ne (this$, this$value)) {
        this$step = (this$ - this$value) / jsfx_specialvar_get_samplesblock();
        return 1.;
      }
      else {
        (this$step = 0.);
        return this$step;
      }
    }();
  }
  //----------------------------------------------------------------------------
  double init$smoother_block_2 (
    double  modulo,
    double& this$,
    double& this$value,
    double& this$step)
  {
    double half = 0.;
    return [&] {
      if (eel2_ne (this$, this$value)) {
        half = modulo * 0.5;
        while (this$ - this$value > half) {
          this$value += modulo;
        }
        while (this$ - this$value < -half) {
          this$value -= modulo;
        }
        this$step = (this$ - this$value) / jsfx_specialvar_get_samplesblock();
        return 1.;
      }
      else {
        (this$step = 0.);
        return this$step;
      }
    }();
  }
  //----------------------------------------------------------------------------
  double init$smoother_init (
    double& this$value,
    double& this$,
    double& this$step)
  {
    this$value = this$;
    this$step  = 0.;
    return 1.;
  }
  //----------------------------------------------------------------------------
  double init$smoother_sample (double& this$value, double& this$step)
  {
    this$value += this$step;
    return this$value;
  }
  //----------------------------------------------------------------------------
  double init$smoother_value (double& this$value) { return this$value; }
}; /* jsfx_process */
}} // namespace artv::geraint_luff

#pragma GCC diagnostic pop
