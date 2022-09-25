#pragma once

#define ARTV_FDNST8_USE_BROKEN_DELAY_TIME_BEHAVIOR 1

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <stdint.h>
#include <utility>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/dsp/third_party/airwindows/all_consoles.hpp"
#include "artv-common/dsp/third_party/airwindows/busscolors4.hpp"

#include "artv-common/dsp/third_party/chokehold/consolidator.hpp"
#include "artv-common/dsp/third_party/chokehold/gate_expander.hpp"
#include "artv-common/dsp/third_party/chokehold/signal_crusher.hpp"
#include "artv-common/dsp/third_party/chokehold/track_comp.hpp"

#include "artv-common/dsp/third_party/chow/phaser.hpp"

#include "artv-common/dsp/third_party/dragonfly/early_reflections.hpp"
#include "artv-common/dsp/third_party/dragonfly/hall.hpp"
#include "artv-common/dsp/third_party/dragonfly/plate.hpp"
#include "artv-common/dsp/third_party/dragonfly/room.hpp"

#include "artv-common/dsp/third_party/geraintluff/atlantis_reverb.hpp"
#include "artv-common/dsp/third_party/geraintluff/echo_cycles.hpp"
#include "artv-common/dsp/third_party/geraintluff/ripple.hpp"
#include "artv-common/dsp/third_party/geraintluff/sandwitch_amp.hpp"
#include "artv-common/dsp/third_party/geraintluff/spring_box.hpp"

#include "artv-common/dsp/third_party/liteon/nonlinear.hpp"
#include "artv-common/dsp/third_party/liteon/sonic_enhancer.hpp"
#include "artv-common/dsp/third_party/liteon/stereotilt.hpp"
#if 0 // sounds broken
#include "artv-common/dsp/third_party/liteon/tubeharmonics.hpp"
#endif

#include "artv-common/dsp/third_party/ljkb/luftikus.hpp"

#include "artv-common/dsp/third_party/rubberband/rubberband.hpp"

#include "artv-common/dsp/third_party/saike/stereo_bub3.hpp"
#if 0
#include "artv-common/dsp/third_party/saike/tanh_aa.hpp"
#endif
#include "artv-common/dsp/third_party/saike/transience.hpp"

#include "artv-common/dsp/third_party/smashed_transistors/ze_big_chorus3.hpp"
#include "artv-common/dsp/third_party/smashed_transistors/ze_little_scanner_chorus.hpp"

#include "artv-common/dsp/third_party/sonic_anomaly/bass_professor1.hpp"
#include "artv-common/dsp/third_party/sonic_anomaly/slax.hpp"
#if 0
#include "artv-common/dsp/third_party/sonic_anomaly/vola2.hpp"
#endif

#include "artv-common/dsp/third_party/soundtouch/soundtouch.hpp"

#include "artv-common/dsp/third_party/sstillwell/1175.hpp"
#include "artv-common/dsp/third_party/sstillwell/4x4.hpp"
#include "artv-common/dsp/third_party/sstillwell/eventhorizon2.hpp"
#include "artv-common/dsp/third_party/sstillwell/fairlychildish.hpp"
#include "artv-common/dsp/third_party/sstillwell/hugebooty.hpp"
#include "artv-common/dsp/third_party/sstillwell/majortom.hpp"
#include "artv-common/dsp/third_party/sstillwell/rbj1073.hpp"

#include "artv-common/dsp/third_party/tal/reverb2.hpp"

#include "artv-common/dsp/third_party/shabtronic/fdn_verb.hpp"

#include "artv-common/dsp/third_party/witti/bbd_echo_stereo.hpp"

#include "artv-common/dsp/own/classes/oversampled.hpp"

#include "artv-common/dsp/own/fx/diffuse_delay.hpp"
#include "artv-common/dsp/own/fx/eq4x.hpp"
#include "artv-common/dsp/own/fx/filter2x.hpp"
#include "artv-common/dsp/own/fx/lin_eq4x.hpp"
#include "artv-common/dsp/own/fx/mod.hpp"
#include "artv-common/dsp/own/fx/phaser.hpp"
#include "artv-common/dsp/own/fx/pitch_shifter.hpp"
#include "artv-common/dsp/own/fx/reverb.hpp"
#include "artv-common/dsp/own/fx/transient_gate.hpp"
#include "artv-common/dsp/own/fx/waveshaper.hpp"

#include "artv-common/dsp/own/fx/sound_delay.hpp"
#if 0
#include "artv-common/dsp/own/fx/polyphase-fir-tester.hpp"
#endif

#include "mix-maxtrix/crossover.hpp"
#include "mix-maxtrix/lin_crossover.hpp"
#include "mix-maxtrix/wonky_crossover.hpp"

#include "artv-common/dsp/own/classes/mix.hpp"
#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/gui_util.hpp"
#include "artv-common/juce/look_and_feel.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/parameters.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "mix-maxtrix/buffers.hpp"
#include "mix-maxtrix/order_and_buffering.hpp"

#define VERSION_INT VERSION_GET (VERSION_MAJOR, VERSION_MINOR, VERSION_REV)

namespace artv { namespace parameters {

// parameter definitions in one place ------------------------------------------
constexpr std::size_t n_stereo_busses = 8;

parameter_cpp_class_define (
  in_selection,
  n_stereo_busses,
  param_common(),
  toggle_buttons_param (
    0,
    make_cstr_array ("1", "2 ", "3", "4", "5", "6 ", "7", "8"),
    16),
  void // defaulted
);

parameter_cpp_class_define (
  out_selection,
  n_stereo_busses,
  param_common(),
  toggle_buttons_param (
    0,
    make_cstr_array ("1", "2 ", "3", "4", "5", "6 ", "7", "8"),
    16),
  void // defaulted
);

parameter_cpp_class_define (
  mixer_sends,
  1,
  param_common(),
  toggle_buttons_param (
    0,
    make_cstr_array (
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">",
      ">"),
    16),
  void // defaulted
);

parameter_cpp_class_define (
  volume,
  n_stereo_busses,
  param_common ("Gain"),
  float_param ("dB", -50.f, 12.f, 0.f, 0.01f, 1.3f),
  slider_ext);

parameter_cpp_class_define (
  global_volume,
  1,
  param_common ("Global Gain"),
  float_param ("dB", -20.f, 20.f, 0.f, 0.01f),
  slider_ext);

struct routing_mode {
  enum {
    full,
    blocks_of_4,
    blocks_of_2,
  };
};
//------------------------------------------------------------------------------
static constexpr uint n_parallel_buses (uint rm)
{
  return parameters::n_stereo_busses >> rm;
}
//------------------------------------------------------------------------------

parameter_cpp_class_define (
  routing,
  1,
  param_common ("Bus Routing"),
  choice_param (
    0,
    make_cstr_array (
      "[1,2,3,4,5,6,7,8]",
      "[1,2,3,4]>[5,6,7,8]",
      "[1,2]>[3,4]>[5,6]>[7,8]"),
    10,
    false),
  combobox_ext);

parameter_cpp_class_define (
  fx_mix,
  n_stereo_busses,
  param_common ("Mix"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::dry_wet_ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wet_pan,
  n_stereo_busses,
  param_common ("Pan"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::wet_pan_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wet_balance,
  n_stereo_busses,
  param_common ("M/S"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::wet_ms_ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  dry_pan,
  n_stereo_busses,
  param_common ("Pan"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::dry_pan_tag {}),
  slider_ext);

parameter_cpp_class_define (
  pan,
  n_stereo_busses,
  param_common ("GPan"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::pan_tag {}),
  slider_ext);

parameter_cpp_class_define (
  dry_balance,
  n_stereo_busses,
  param_common ("M/S"),
  dry_wet_mixer::get_parameter (dry_wet_mixer::dry_ms_ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mute_solo,
  n_stereo_busses,
  param_common(),
  toggle_buttons_param (0, make_cstr_array ("Mute", "Solo"), 4),
  void // defaulted
);

// clang-format off
static constexpr auto channel_modifs_names = make_cstr_array (
  "P", "L ", "R", "S");
// clang-format on

static constexpr uint n_channel_modifiers = channel_modifs_names.size();

parameter_cpp_class_define (
  channel_modifs,
  n_stereo_busses,
  param_common(),
  toggle_buttons_param (0, channel_modifs_names, 16),
  void // defaulted
);
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  busscolors4_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<airwindows::busscolors4>(),
    declptr<airwindows::busscolors4::drive_tag>()),
  airwindows::busscolors4::get_parameter (
    airwindows::busscolors4::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  busscolors4_color,
  n_stereo_busses,
  param_common (
    "Color",
    declptr<airwindows::busscolors4>(),
    declptr<airwindows::busscolors4::color_tag>()),
  airwindows::busscolors4::get_parameter (
    airwindows::busscolors4::color_tag {}),
  slider_ext);

parameter_cpp_class_define (
  busscolors4_dry_wet,
  n_stereo_busses,
  param_common (
    "Dry/Wet",
    declptr<airwindows::busscolors4>(),
    declptr<airwindows::busscolors4::dry_wet_tag>()),
  airwindows::busscolors4::get_parameter (
    airwindows::busscolors4::dry_wet_tag {}),
  slider_ext);

using busscolors4_params
  = mp_list<busscolors4_drive, busscolors4_color, busscolors4_dry_wet>;
//------------------------------------------------------------------------------
static constexpr uint crossover_n_bands = 4;
static constexpr uint n_crossovers      = crossover_n_bands - 1;

using mm_lr_crossv = mixmaxtrix_lr_crossover<crossover_n_bands>;

parameter_cpp_class_define (
  lr_crossv_band1_frequency,
  1,
  param_common (
    "1.Freq",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band1_freq_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band2_frequency,
  1,
  param_common (
    "2.Freq",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band2_freq_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band3_frequency,
  1,
  param_common (
    "3.Freq",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band3_freq_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band3_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band1_mode,
  1,
  param_common (
    "1.Mode",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band1_mode_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band1_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band2_mode,
  1,
  param_common (
    "2.Mode",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band2_mode_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band2_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band3_mode,
  1,
  param_common (
    "3.Mode",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band3_mode_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band3_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band1_diff,
  1,
  param_common (
    "1.LR-Diff",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band1_diff_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band1_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band2_diff,
  1,
  param_common (
    "2.LR-Diff",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band2_diff_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band2_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band3_diff,
  1,
  param_common (
    "3.LR-Diff",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band3_diff_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band3_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band1_out,
  1,
  param_common (
    "1.Send",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band1_out_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band1_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band2_out,
  1,
  param_common (
    "2.Send",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band2_out_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band2_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lr_crossv_band3_out,
  1,
  param_common (
    "3.Send",
    declptr<mm_lr_crossv>(),
    declptr<mm_lr_crossv::band3_out_tag>()),
  mm_lr_crossv::get_parameter (mm_lr_crossv::band3_out_tag {}),
  slider_ext);

using lr_crossv_params = mp_list<
  lr_crossv_band1_frequency,
  lr_crossv_band1_diff,
  lr_crossv_band1_mode,
  lr_crossv_band1_out,
  lr_crossv_band2_frequency,
  lr_crossv_band2_diff,
  lr_crossv_band2_mode,
  lr_crossv_band2_out,
  lr_crossv_band3_frequency,
  lr_crossv_band3_diff,
  lr_crossv_band3_mode,
  lr_crossv_band3_out>;
//------------------------------------------------------------------------------
using mm_wonky_crossv = mixmaxtrix_wonky_crossover<crossover_n_bands>;

parameter_cpp_class_define (
  wonky_crossv_band1_frequency,
  1,
  param_common (
    "1.Freq",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band1_freq_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band2_frequency,
  1,
  param_common (
    "2.Freq",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band2_freq_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band3_frequency,
  1,
  param_common (
    "3.Freq",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band3_freq_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band3_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band1_mode,
  1,
  param_common (
    "1.Mode",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band1_mode_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band1_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band2_mode,
  1,
  param_common (
    "2.Mode",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band2_mode_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band2_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band3_mode,
  1,
  param_common (
    "3.Mode",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band3_mode_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band3_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band1_diff,
  1,
  param_common (
    "1.LR-Diff",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band1_diff_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band1_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band2_diff,
  1,
  param_common (
    "2.LR-Diff",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band2_diff_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band2_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band3_diff,
  1,
  param_common (
    "3.LR-Diff",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band3_diff_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band3_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band1_out,
  1,
  param_common (
    "1.Send",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band1_out_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band1_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band2_out,
  1,
  param_common (
    "2.Send",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band2_out_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band2_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wonky_crossv_band3_out,
  1,
  param_common (
    "3.Send",
    declptr<mm_wonky_crossv>(),
    declptr<mm_wonky_crossv::band3_out_tag>()),
  mm_wonky_crossv::get_parameter (mm_wonky_crossv::band3_out_tag {}),
  slider_ext);

using wonky_crossv_params = mp_list<
  wonky_crossv_band1_frequency,
  wonky_crossv_band1_diff,
  wonky_crossv_band1_mode,
  wonky_crossv_band1_out,
  wonky_crossv_band2_frequency,
  wonky_crossv_band2_diff,
  wonky_crossv_band2_mode,
  wonky_crossv_band2_out,
  wonky_crossv_band3_frequency,
  wonky_crossv_band3_diff,
  wonky_crossv_band3_mode,
  wonky_crossv_band3_out>;

//------------------------------------------------------------------------------
using mm_lin_iir_crossv = mixmaxtrix_linphase_iir_crossover<crossover_n_bands>;

parameter_cpp_class_define (
  lin_iir_crossv_band1_frequency,
  1,
  param_common (
    "1.Freq",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band1_freq_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band2_frequency,
  1,
  param_common (
    "2.Freq",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band2_freq_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band3_frequency,
  1,
  param_common (
    "3.Freq",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band3_freq_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band3_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band1_mode,
  1,
  param_common (
    "1.Mode",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band1_mode_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band1_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band2_mode,
  1,
  param_common (
    "2.Mode",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band2_mode_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band2_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band3_mode,
  1,
  param_common (
    "3.Mode",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band3_mode_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band3_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band1_diff,
  1,
  param_common (
    "1.LR-Diff",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band1_diff_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band1_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band2_diff,
  1,
  param_common (
    "2.LR-Diff",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band2_diff_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band2_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band3_diff,
  1,
  param_common (
    "3.LR-Diff",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band3_diff_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band3_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band1_out,
  1,
  param_common (
    "1.Send",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band1_out_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band1_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band2_out,
  1,
  param_common (
    "2.Send",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band2_out_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band2_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band3_out,
  1,
  param_common (
    "3.Send",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band3_out_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band3_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band1_quality,
  1,
  param_common (
    "1.Quality",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band1_quality_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band1_quality_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band2_quality,
  1,
  param_common (
    "2.Quality",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band2_quality_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band2_quality_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_iir_crossv_band3_quality,
  1,
  param_common (
    "3.Quality",
    declptr<mm_lin_iir_crossv>(),
    declptr<mm_lin_iir_crossv::band3_quality_tag>()),
  mm_lin_iir_crossv::get_parameter (mm_lin_iir_crossv::band3_quality_tag {}),
  slider_ext);

using lin_iir_crossv_params = mp_list<
  lin_iir_crossv_band1_frequency,
  lin_iir_crossv_band1_diff,
  lin_iir_crossv_band1_mode,
  lin_iir_crossv_band1_out,
  lin_iir_crossv_band2_frequency,
  lin_iir_crossv_band2_diff,
  lin_iir_crossv_band2_mode,
  lin_iir_crossv_band2_out,
  lin_iir_crossv_band3_frequency,
  lin_iir_crossv_band3_diff,
  lin_iir_crossv_band3_mode,
  lin_iir_crossv_band3_out,
  lin_iir_crossv_band1_quality,
  lin_iir_crossv_band2_quality,
  lin_iir_crossv_band3_quality>;
//------------------------------------------------------------------------------
static constexpr uint console_n_sends = 1 + (crossover_n_bands - 1);
static constexpr uint console_n_elems = n_stereo_busses + console_n_sends;

using consoles_dsp_type = airwindows::all_consoles<console_n_elems>;

parameter_cpp_class_define (
  consoles_type,
  n_stereo_busses,
  param_common (
    "Type",
    declptr<consoles_dsp_type>(),
    declptr<consoles_dsp_type::type_tag>()),
  consoles_dsp_type::get_parameter (consoles_dsp_type::type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consoles_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<consoles_dsp_type>(),
    declptr<consoles_dsp_type::drive_tag>()),
  consoles_dsp_type::get_parameter (consoles_dsp_type::drive_tag {}),
  slider_ext);

using consoles_params = mp_list<consoles_type, consoles_drive>;
//------------------------------------------------------------------------------
#if 0
parameter_cpp_class_define (
  density2_density,
  n_stereo_busses,
  param_common (
    "Density",
    declptr<airwindows::density2>(),
    declptr<airwindows::density2::density_tag>()),
  airwindows::density2::get_parameter (airwindows::density2::density_tag {}),
  slider_ext);

parameter_cpp_class_define (
  density2_highpass,
  n_stereo_busses,
  param_common (
    "Highpass",
    declptr<airwindows::density2>(),
    declptr<airwindows::density2::highpass_tag>()),
  airwindows::density2::get_parameter (airwindows::density2::highpass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  density2_dry_wet,
  n_stereo_busses,
  param_common (
    "Dry/Wet",
    declptr<airwindows::density2>(),
    declptr<airwindows::density2::dry_wet_tag>()),
  airwindows::density2::get_parameter (airwindows::density2::dry_wet_tag {}),
  slider_ext);

using density2_params = mp_list<density2_density, density2_highpass>;
#endif
//------------------------------------------------------------------------------
#if 0
parameter_cpp_class_define (
  focus_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<airwindows::focus>(),
    declptr<airwindows::focus::drive_tag>()),
  airwindows::focus::get_parameter (airwindows::focus::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  focus_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<airwindows::focus>(),
    declptr<airwindows::focus::mode_tag>()),
  airwindows::focus::get_parameter (airwindows::focus::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  focus_focus,
  n_stereo_busses,
  param_common (
    "Focus",
    declptr<airwindows::focus>(),
    declptr<airwindows::focus::focus_tag>()),
  airwindows::focus::get_parameter (airwindows::focus::focus_tag {}),
  slider_ext);

parameter_cpp_class_define (
  focus_dry_wet,
  n_stereo_busses,
  param_common (
    "Dry/Wet",
    declptr<airwindows::focus>(),
    declptr<airwindows::focus::dry_wet_tag>()),
  airwindows::focus::get_parameter (airwindows::focus::dry_wet_tag {}),
  slider_ext);

using focus_params
  = mp_list<focus_drive, focus_mode, focus_focus, focus_dry_wet>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  stereo_bub3_delay,
  n_stereo_busses,
  param_common (
    "Delay",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::delay_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::delay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_lowcut,
  n_stereo_busses,
  param_common (
    "Crossovr",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::lowcut_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::lowcut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_side_level,
  n_stereo_busses,
  param_common (
    "Old Side",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::sidelevel_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::sidelevel_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_strength,
  n_stereo_busses,
  param_common (
    "Strength",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::strength_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::strength_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_vibratoamount,
  n_stereo_busses,
  param_common (
    "Vib Amt",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::vibratoamount_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::vibratoamount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_vibratospeed,
  n_stereo_busses,
  param_common (
    "Vib Speed",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::vibratospeed_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::vibratospeed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_saturation,
  n_stereo_busses,
  param_common (
    "Sat",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::saturation_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::saturation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_bub3_sidehp,
  n_stereo_busses,
  param_common (
    "Side HP",
    declptr<saike::stereo_bub3>(),
    declptr<saike::stereo_bub3::sidehp_tag>()),
  saike::stereo_bub3::get_parameter (saike::stereo_bub3::sidehp_tag {}),
  slider_ext);

using stereo_bub3_params = mp_list<
  stereo_bub3_delay,
  stereo_bub3_strength,
  stereo_bub3_lowcut,
  stereo_bub3_side_level,
  stereo_bub3_vibratospeed,
  stereo_bub3_vibratoamount,
  stereo_bub3_saturation,
  stereo_bub3_sidehp>;
//------------------------------------------------------------------------------
#if 0
parameter_cpp_class_define (
  smooth_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<saike::smooth>(),
    declptr<saike::smooth::drive_tag>()),
  saike::smooth::get_parameter (saike::smooth::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  smooth_slew,
  n_stereo_busses,
  param_common (
    "Slew",
    declptr<saike::smooth>(),
    declptr<saike::smooth::slew_tag>()),
  saike::smooth::get_parameter (saike::smooth::slew_tag {}),
  slider_ext);

parameter_cpp_class_define (
  smooth_ceiling,
  n_stereo_busses,
  param_common (
    "Ceiling",
    declptr<saike::smooth>(),
    declptr<saike::smooth::ceiling_tag>()),
  saike::smooth::get_parameter (saike::smooth::ceiling_tag {}),
  slider_ext);

parameter_cpp_class_define (
  smooth_warmth,
  n_stereo_busses,
  param_common (
    "Warmth",
    declptr<saike::smooth>(),
    declptr<saike::smooth::warmth_tag>()),
  saike::smooth::get_parameter (saike::smooth::warmth_tag {}),
  slider_ext);

using smooth_params
  = mp_list<smooth_drive, smooth_slew, smooth_ceiling, smooth_warmth>;
#endif
//------------------------------------------------------------------------------
#if 0
parameter_cpp_class_define (
  tanh_aa_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::gain_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_ceiling,
  n_stereo_busses,
  param_common (
    "Ceiling",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::ceiling_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::ceiling_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_c_alias,
  n_stereo_busses,
  param_common (
    "AA Mode",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::c_alias_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::c_alias_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_fix_dc,
  n_stereo_busses,
  param_common (
    "Fix DC",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::fix_dc_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::fix_dc_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_smooth_time,
  n_stereo_busses,
  param_common (
    "Intertia",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::smooth_time_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::smooth_time_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_shelf,
  n_stereo_busses,
  param_common (
    "HF Correct",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::shelf_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::shelf_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_func,
  n_stereo_busses,
  param_common (
    "Function",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::func_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::func_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tanh_aa_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<saike::tanh_aa>(),
    declptr<saike::tanh_aa::oversampling_tag>()),
  saike::tanh_aa::get_parameter (saike::tanh_aa::oversampling_tag {}),
  slider_ext);

using tanh_aa_params = mp11::mp_list<
  tanh_aa_func,
  tanh_aa_gain,
  tanh_aa_ceiling,
  tanh_aa_smooth_time,
  tanh_aa_c_alias,
  tanh_aa_oversampling,
  tanh_aa_fix_dc,
  tanh_aa_shelf>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  bbe_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<liteon::bbe>(),
    declptr<liteon::bbe::drive_tag>()),
  liteon::bbe::get_parameter (liteon::bbe::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbe_low_contour,
  n_stereo_busses,
  param_common (
    "L Contour",
    declptr<liteon::bbe>(),
    declptr<liteon::bbe::low_contour_tag>()),
  liteon::bbe::get_parameter (liteon::bbe::low_contour_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbe_process,
  n_stereo_busses,
  param_common (
    "H Contour",
    declptr<liteon::bbe>(),
    declptr<liteon::bbe::process_tag>()),
  liteon::bbe::get_parameter (liteon::bbe::process_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbe_cv,
  n_stereo_busses,
  param_common (
    "Peak Det.",
    declptr<liteon::bbe>(),
    declptr<liteon::bbe::cv_tag>()),
  liteon::bbe::get_parameter (liteon::bbe::cv_tag {}),
  slider_ext);

using bbe_params = mp_list<bbe_drive, bbe_low_contour, bbe_process, bbe_cv>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  nonlinear_saturation,
  n_stereo_busses,
  param_common (
    "Sat",
    declptr<liteon::nonlinear>(),
    declptr<liteon::nonlinear::saturation_tag>()),
  liteon::nonlinear::get_parameter (liteon::nonlinear::saturation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  nonlinear_fluctuation,
  n_stereo_busses,
  param_common (
    "Fluctu",
    declptr<liteon::nonlinear>(),
    declptr<liteon::nonlinear::fluctuation_tag>()),
  liteon::nonlinear::get_parameter (liteon::nonlinear::fluctuation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  nonlinear_noise_floor,
  n_stereo_busses,
  param_common (
    "Noise Flr",
    declptr<liteon::nonlinear>(),
    declptr<liteon::nonlinear::noise_floor_tag>()),
  liteon::nonlinear::get_parameter (liteon::nonlinear::noise_floor_tag {}),
  slider_ext);

parameter_cpp_class_define (
  nonlinear_drive,
  n_stereo_busses,
  param_common (
    "Output",
    declptr<liteon::nonlinear>(),
    declptr<liteon::nonlinear::output_tag>()),
  liteon::nonlinear::get_parameter (liteon::nonlinear::output_tag {}),
  slider_ext);

using nonlinear_params = mp11::
  mp_list<nonlinear_saturation, nonlinear_fluctuation, nonlinear_noise_floor>;
//------------------------------------------------------------------------------
#if ADD_2ND_TIER_FX
parameter_cpp_class_define (
  pseudostereo_amount,
  n_stereo_busses,
  param_common (
    "Amount",
    declptr<liteon::pseudostereo>(),
    declptr<liteon::pseudostereo::amount_tag>()),
  liteon::pseudostereo::get_parameter (liteon::pseudostereo::amount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  pseudostereo_delay,
  n_stereo_busses,
  param_common (
    "Delay",
    declptr<liteon::pseudostereo>(),
    declptr<liteon::pseudostereo::delay_tag>()),
  liteon::pseudostereo::get_parameter (liteon::pseudostereo::delay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  pseudostereo_balance,
  n_stereo_busses,
  param_common (
    "Balance",
    declptr<liteon::pseudostereo>(),
    declptr<liteon::pseudostereo::balance_tag>()),
  liteon::pseudostereo::get_parameter (liteon::pseudostereo::balance_tag {}),
  slider_ext);

using pseudostereo_params = mp11::
  mp_list<pseudostereo_amount, pseudostereo_delay, pseudostereo_balance>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  stereo_tilt_frequency,
  n_stereo_busses,
  param_common (
    "Freq",
    declptr<liteon::stereo_tilt>(),
    declptr<liteon::stereo_tilt::frequency_tag>()),
  liteon::stereo_tilt::get_parameter (liteon::stereo_tilt::frequency_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_tilt_tilt,
  n_stereo_busses,
  param_common (
    "Tilt",
    declptr<liteon::stereo_tilt>(),
    declptr<liteon::stereo_tilt::tilt_tag>()),
  liteon::stereo_tilt::get_parameter (liteon::stereo_tilt::tilt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo_tilt_balance,
  n_stereo_busses,
  param_common (
    "Balance",
    declptr<liteon::stereo_tilt>(),
    declptr<liteon::stereo_tilt::balance_tag>()),
  liteon::stereo_tilt::get_parameter (liteon::stereo_tilt::balance_tag {}),
  slider_ext);

using stereo_tilt_params
  = mp_list<stereo_tilt_frequency, stereo_tilt_tilt, stereo_tilt_balance>;
//------------------------------------------------------------------------------
#if 0 // Sounds broken
parameter_cpp_class_define (
  tube_harmonics_even,
  n_stereo_busses,
  param_common (
    "Even H",
    declptr<liteon::tube_harmonics>(),
    declptr<liteon::tube_harmonics::even_tag>()),
  liteon::tube_harmonics::get_parameter (liteon::tube_harmonics::even_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tube_harmonics_odd,
  n_stereo_busses,
  param_common (
    "Odd H",
    declptr<liteon::tube_harmonics>(),
    declptr<liteon::tube_harmonics::odd_tag>()),
  liteon::tube_harmonics::get_parameter (liteon::tube_harmonics::odd_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tube_harmonics_fluctuation,
  n_stereo_busses,
  param_common (
    "Fluctuation",
    declptr<liteon::tube_harmonics>(),
    declptr<liteon::tube_harmonics::fluctuation_tag>()),
  liteon::tube_harmonics::get_parameter (
    liteon::tube_harmonics::fluctuation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tube_harmonics_ts_input,
  n_stereo_busses,
  param_common (
    "Ts In",
    declptr<liteon::tube_harmonics>(),
    declptr<liteon::tube_harmonics::ts_input_tag>()),
  liteon::tube_harmonics::get_parameter (
    liteon::tube_harmonics::ts_input_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tube_harmonics_ts_output,
  n_stereo_busses,
  param_common (
    "Ts Out",
    declptr<liteon::tube_harmonics>(),
    declptr<liteon::tube_harmonics::ts_output_tag>()),
  liteon::tube_harmonics::get_parameter (
    liteon::tube_harmonics::ts_output_tag {}),
  slider_ext);

using tube_harmonics_params = mp_list<
  tube_harmonics_even,
  tube_harmonics_odd,
  tube_harmonics_ts_input,
  tube_harmonics_ts_output,
  tube_harmonics_fluctuation>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  transience_mode,
  n_stereo_busses,
  param_common (
    "Detector",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::mode_tag>()),
  saike::transience::get_parameter (saike::transience::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_attack,
  n_stereo_busses,
  param_common (
    "Attk T",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::sattack_tag>()),
  saike::transience::get_parameter (saike::transience::sattack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_attack_amt,
  n_stereo_busses,
  param_common (
    "Attk Amt",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::strength_tag>()),
  saike::transience::get_parameter (saike::transience::strength_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_decay,
  n_stereo_busses,
  param_common (
    "Dec Time",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::sdecay_tag>()),
  saike::transience::get_parameter (saike::transience::sdecay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_decay_amt,
  n_stereo_busses,
  param_common (
    "Dec Amt",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::strength2_tag>()),
  saike::transience::get_parameter (saike::transience::strength2_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_gainsmoothing,
  n_stereo_busses,
  param_common (
    "G Smooth",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::gainsmoothing_tag>()),
  saike::transience::get_parameter (saike::transience::gainsmoothing_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_gaincompensation,
  n_stereo_busses,
  param_common (
    "G Compen",
    declptr<updownsampled<saike::transience>>(),
    declptr<saike::transience::compensate_tag>()),
  saike::transience::get_parameter (saike::transience::compensate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transience_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<saike::transience>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<saike::transience>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using transience_params = mp_list<
  transience_attack,
  transience_attack_amt,
  transience_decay,
  transience_decay_amt,
  transience_mode,
  transience_gainsmoothing,
  transience_gaincompensation,
  transience_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  slax_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<sonic_anomaly::slax::gain_tag>()),
  sonic_anomaly::slax::get_parameter (sonic_anomaly::slax::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  slax_peak,
  n_stereo_busses,
  param_common (
    "Peak",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<sonic_anomaly::slax::peak_tag>()),
  sonic_anomaly::slax::get_parameter (sonic_anomaly::slax::peak_tag {}),
  slider_ext);

parameter_cpp_class_define (
  slax_emphasis,
  n_stereo_busses,
  param_common (
    "Empha",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<sonic_anomaly::slax::emphasis_tag>()),
  sonic_anomaly::slax::get_parameter (sonic_anomaly::slax::emphasis_tag {}),
  slider_ext);

parameter_cpp_class_define (
  slax_ratio,
  n_stereo_busses,
  param_common (
    "Ratio",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<sonic_anomaly::slax::comp_lim_tag>()),
  sonic_anomaly::slax::get_parameter (sonic_anomaly::slax::comp_lim_tag {}),
  slider_ext);

parameter_cpp_class_define (
  slax_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<sonic_anomaly::slax::mode_tag>()),
  sonic_anomaly::slax::get_parameter (sonic_anomaly::slax::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  slax_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<sonic_anomaly::slax>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<sonic_anomaly::slax>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using slax_params = mp_list<
  slax_gain,
  slax_peak,
  slax_ratio,
  slax_mode,
  slax_emphasis,
  slax_oversampling>;
//------------------------------------------------------------------------------
#if 0
// Has latency...
parameter_cpp_class_define (
  vola2_sc_filter,
  n_stereo_busses,
  param_common (
    "SC Filter",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::sc_filter_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::sc_filter_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_lf,
  n_stereo_busses,
  param_common (
    "LF Cut",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::lf_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::lf_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_attack,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::attack_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::attack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_push_down,
  n_stereo_busses,
  param_common (
    "Push Down",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::push_down_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::push_down_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_pull_up,
  n_stereo_busses,
  param_common (
    "Pull Up",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::pull_up_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::pull_up_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_recovery,
  n_stereo_busses,
  param_common (
    "Recovery",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::recovery_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::recovery_tag {}),
  slider_ext);

parameter_cpp_class_define (
  vola2_emphasis,
  n_stereo_busses,
  param_common (
    "Empha",
    declptr<sonic_anomaly::vola2>(),
    declptr<sonic_anomaly::vola2::emphasis_tag>()),
  sonic_anomaly::vola2::get_parameter (sonic_anomaly::vola2::emphasis_tag {}),
  slider_ext);

using vola2_params = mp_list<
  vola2_sc_filter,
  vola2_lf,
  vola2_attack,
  vola2_push_down,
  vola2_pull_up,
  vola2_recovery,
  vola2_emphasis>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  bass_professor_amount,
  n_stereo_busses,
  param_common (
    "Amount",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::amount_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::amount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_presence,
  n_stereo_busses,
  param_common (
    "Presence",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::presence_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::presence_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_dirt,
  n_stereo_busses,
  param_common (
    "Dirt",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::dirt_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::dirt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_bass,
  n_stereo_busses,
  param_common (
    "Bass",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::bass_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::bass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_middle,
  n_stereo_busses,
  param_common (
    "Middle",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::middle_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::middle_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_treble,
  n_stereo_busses,
  param_common (
    "Treble",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::treble_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::treble_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_depth,
  n_stereo_busses,
  param_common (
    "Depth",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::depth_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bass_professor_lf_cut,
  n_stereo_busses,
  param_common (
    "LF Cut",
    declptr<sonic_anomaly::bass_professor>(),
    declptr<sonic_anomaly::bass_professor::lf_cut_tag>()),
  sonic_anomaly::bass_professor::get_parameter (
    sonic_anomaly::bass_professor::lf_cut_tag {}),
  slider_ext);

using bass_professor_params = mp_list<
  bass_professor_amount,
  bass_professor_presence,
  bass_professor_bass,
  bass_professor_dirt,
  bass_professor_middle,
  bass_professor_depth,
  bass_professor_treble,
  bass_professor_lf_cut>;
//------------------------------------------------------------------------------
#if ADD_2ND_TIER_FX
parameter_cpp_class_define (
  transpire_sensitivity,
  n_stereo_busses,
  param_common (
    "Sensitivity",
    declptr<sonic_anomaly::transpire>(),
    declptr<sonic_anomaly::transpire::sensitivity_tag>()),
  sonic_anomaly::transpire::get_parameter (
    sonic_anomaly::transpire::sensitivity_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transpire_attack,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<sonic_anomaly::transpire>(),
    declptr<sonic_anomaly::transpire::attack_tag>()),
  sonic_anomaly::transpire::get_parameter (
    sonic_anomaly::transpire::attack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transpire_sustain,
  n_stereo_busses,
  param_common (
    "Sustain",
    declptr<sonic_anomaly::transpire>(),
    declptr<sonic_anomaly::transpire::sustain_tag>()),
  sonic_anomaly::transpire::get_parameter (
    sonic_anomaly::transpire::sustain_tag {}),
  slider_ext);

using transpire_params
  = mp_list<transpire_sensitivity, transpire_attack, transpire_sustain>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  atlantis_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::dry_db_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::dry_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_wet,
  n_stereo_busses,
  param_common (
    "Wet",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::wet_db_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::wet_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_decay,
  n_stereo_busses,
  param_common (
    "Decay",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::decay_seconds_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::decay_seconds_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_window,
  n_stereo_busses,
  param_common (
    "Window",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::window_ms_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::window_ms_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_shimmer,
  n_stereo_busses,
  param_common (
    "Shimmer",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::shimmer_factor_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::shimmer_factor_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_shimmer_tone,
  n_stereo_busses,
  param_common (
    "Shm 5ths",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::shimmer_tone_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::shimmer_tone_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_detune_shift,
  n_stereo_busses,
  param_common (
    "Dtune Bias",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::detune_shift_bias_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::detune_shift_bias_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_detune_spd,
  n_stereo_busses,
  param_common (
    "Dtune Spd",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::detune_cents_per_second_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::detune_cents_per_second_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_compressor_threshold,
  n_stereo_busses,
  param_common (
    "Comp Thres",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::compressor_threshold_db_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::compressor_threshold_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_compressor_ratio,
  n_stereo_busses,
  param_common (
    "Comp Ratio",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::compressor_ratio_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::compressor_ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_damping,
  n_stereo_busses,
  param_common (
    "Damp Amt",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::damping_strength_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::damping_strength_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_damping_low,
  n_stereo_busses,
  param_common (
    "Damp L",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::low_damping_hz_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::low_damping_hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_damping_high,
  n_stereo_busses,
  param_common (
    "Damp H",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::high_damping_hz_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::high_damping_hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::ducking_speed_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  atlantis_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<geraint_luff::atlantis_reverb>(),
    declptr<geraint_luff::atlantis_reverb::ducking_threshold_tag>()),
  geraint_luff::atlantis_reverb::get_parameter (
    geraint_luff::atlantis_reverb::ducking_threshold_tag {}),
  slider_ext);

using atlantis_reverb_params = mp_list<
  // atlantis_dry,
  // atlantis_wet,
  atlantis_decay,
  atlantis_window,
  atlantis_shimmer,
  atlantis_shimmer_tone,

  atlantis_detune_shift,
  atlantis_detune_spd,
  atlantis_compressor_threshold,
  atlantis_compressor_ratio,

  atlantis_ducking_threshold,
  atlantis_ducking_speed,

  atlantis_damping,
  atlantis_damping_low,
  atlantis_damping_high>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  echo_cycles_input_width,
  n_stereo_busses,
  param_common (
    "In Width",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::input_width_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::input_width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_output_variation,
  n_stereo_busses,
  param_common (
    "Out Var",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::output_variation_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::output_variation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_delay_ms,
  n_stereo_busses,
  param_common (
    "Delay",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::delay_ms_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::delay_ms_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_delay_beats,
  n_stereo_busses,
  param_common (
    "Del Sync",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::delay_beats_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::delay_beats_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_feedback_ratio,
  n_stereo_busses,
  param_common (
    "Fb Ratio",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::feedback_ratio_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::feedback_ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_feedback_rotation,
  n_stereo_busses,
  param_common (
    "Fb Rot",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::feedback_rotation_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::feedback_rotation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_input_rotation_initial,
  n_stereo_busses,
  param_common (
    "In Rot",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::input_rotation_initial_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::input_rotation_initial_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_output_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::output_dry_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::output_dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_output_wet,
  n_stereo_busses,
  param_common (
    "Wet",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::output_wet_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::output_wet_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_filter_freq,
  n_stereo_busses,
  param_common (
    "Filt Frq",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::filter_freq_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::filter_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_filter_db,
  n_stereo_busses,
  param_common (
    "Filt Drv",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::filter_db_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::filter_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_filter_bandwidth,
  n_stereo_busses,
  param_common (
    "Filt BW",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::filter_bandwidth_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::filter_bandwidth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_rotation_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::rotation_mode_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::rotation_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::ducking_speed_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  echo_cycles_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<geraint_luff::echo_cycles>(),
    declptr<geraint_luff::echo_cycles::ducking_threshold_tag>()),
  geraint_luff::echo_cycles::get_parameter (
    geraint_luff::echo_cycles::ducking_threshold_tag {}),
  slider_ext);
using echo_cycles_params = mp_list<
  // echo_cycles_output_dry,
  // echo_cycles_output_wet,
  echo_cycles_delay_ms,
  echo_cycles_delay_beats,
  echo_cycles_feedback_ratio,
  echo_cycles_feedback_rotation,

  echo_cycles_input_rotation_initial,
  echo_cycles_rotation_mode,
  echo_cycles_input_width,
  echo_cycles_output_variation,

  echo_cycles_ducking_threshold,
  echo_cycles_ducking_speed,
  echo_cycles_filter_freq,
  echo_cycles_filter_db,
  echo_cycles_filter_bandwidth>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  ripple_control_band_count,
  n_stereo_busses,
  param_common (
    "Bands",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::control_band_count_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::control_band_count_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_low_freq,
  n_stereo_busses,
  param_common (
    "Low Freq",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::low_freq_tag>()),
  geraint_luff::ripple::get_parameter (geraint_luff::ripple::low_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_high_freq,
  n_stereo_busses,
  param_common (
    "High Freq",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::high_freq_tag>()),
  geraint_luff::ripple::get_parameter (geraint_luff::ripple::high_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_phase_offset,
  n_stereo_busses,
  param_common (
    "Ph Offset",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::cycle_phase_offset_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::cycle_phase_offset_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_phase_stereo_offset,
  n_stereo_busses,
  param_common (
    "Ph Stereo",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::cycle_phase_stereo_offset_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::cycle_phase_stereo_offset_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_phase_lfo_hz,
  n_stereo_busses,
  param_common (
    "LFO",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::cycle_phase_lfo_hz_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::cycle_phase_lfo_hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_phase_invbeat,
  n_stereo_busses,
  param_common (
    "LFO Sync",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::cycle_phase_invbeat_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::cycle_phase_invbeat_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_control_filter_mode,
  n_stereo_busses,
  param_common (
    "Filt Type",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::control_filter_mode_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::control_filter_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_filter_db,
  n_stereo_busses,
  param_common (
    "Filt Strength",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::filter_db_tag>()),
  geraint_luff::ripple::get_parameter (geraint_luff::ripple::filter_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_filter_width_factor,
  n_stereo_busses,
  param_common (
    "Filt Width",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::filter_width_factor_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::filter_width_factor_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_output_gain_db,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<geraint_luff::ripple::output_gain_db_tag>()),
  geraint_luff::ripple::get_parameter (
    geraint_luff::ripple::output_gain_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ripple_oversampling,
  n_stereo_busses,
  param_common (
    "Upsample",
    declptr<upsampled<geraint_luff::ripple>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<geraint_luff::ripple>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using ripple_params = mp_list<
  ripple_control_band_count,
  ripple_control_filter_mode,
  ripple_low_freq,
  ripple_high_freq,
  ripple_phase_offset,
  ripple_phase_stereo_offset,
  ripple_phase_lfo_hz,
  ripple_phase_invbeat,
  ripple_filter_db,
  ripple_filter_width_factor,
  ripple_output_gain_db,
  ripple_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  spring_box_density,
  n_stereo_busses,
  param_common (
    "Density",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::density_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::density_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_density_sync,
  n_stereo_busses,
  param_common (
    "Dens Sync",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::density_sync_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::density_sync_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::feedback_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_late_bias,
  n_stereo_busses,
  param_common (
    "Late Bias",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::late_bias_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::late_bias_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_shape,
  n_stereo_busses,
  param_common (
    "Shape",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::shape_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::shape_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_filter,
  n_stereo_busses,
  param_common (
    "Filter",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::filter_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::filter_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_detune,
  n_stereo_busses,
  param_common (
    "Detune",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::detune_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::detune_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_speed,
  n_stereo_busses,
  param_common (
    "Ch Speed",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::speed_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_chorus_alignment,
  n_stereo_busses,
  param_common (
    "Ch Align",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::chorus_alignment_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::chorus_alignment_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::ducking_speed_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  spring_box_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<geraint_luff::spring_box>(),
    declptr<geraint_luff::spring_box::ducking_threshold_tag>()),
  geraint_luff::spring_box::get_parameter (
    geraint_luff::spring_box::ducking_threshold_tag {}),
  slider_ext);

using spring_box_params = mp_list<
  spring_box_density,
  spring_box_density_sync,
  spring_box_feedback,
  spring_box_late_bias,

  spring_box_shape,
  spring_box_filter,
  spring_box_ducking_threshold,
  spring_box_ducking_speed,

  spring_box_detune,
  spring_box_speed,
  spring_box_chorus_alignment>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  sandwitch_amp_limit_db,
  n_stereo_busses,
  param_common (
    "Limit",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::limit_db_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::limit_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_asymmetry,
  n_stereo_busses,
  param_common (
    "Asym",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::asymmetry_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::asymmetry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_output_db,
  n_stereo_busses,
  param_common (
    "Output",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::output_db_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::output_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_distortion_width,
  n_stereo_busses,
  param_common (
    "Dist Wd",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::distortion_width_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::distortion_width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_filter_freq,
  n_stereo_busses,
  param_common (
    "Filter F",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::filter_freq_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::filter_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_filter_gain,
  n_stereo_busses,
  param_common (
    "Filter G",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::filter_gain_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::filter_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_bandwidth_octaves,
  n_stereo_busses,
  param_common (
    "BW",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::bandwidth_octaves_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::bandwidth_octaves_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sandwitch_amp_secondary_gain_db,
  n_stereo_busses,
  param_common (
    "Sec Gain",
    declptr<geraint_luff::sandwitch_amp>(),
    declptr<geraint_luff::sandwitch_amp::secondary_gain_db_tag>()),
  geraint_luff::sandwitch_amp::get_parameter (
    geraint_luff::sandwitch_amp::secondary_gain_db_tag {}),
  slider_ext);

using sandwitch_amp_params = mp_list<
  sandwitch_amp_limit_db,
  sandwitch_amp_output_db,
  sandwitch_amp_asymmetry,
  sandwitch_amp_distortion_width,
  sandwitch_amp_filter_freq,
  sandwitch_amp_filter_gain,
  sandwitch_amp_bandwidth_octaves,
  sandwitch_amp_secondary_gain_db>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  bbd_echo_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::feedback_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_lfo_depth,
  n_stereo_busses,
  param_common (
    "Lfo Amt",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::lfo_depth_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::lfo_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_lfo_speed,
  n_stereo_busses,
  param_common (
    "Lfo Spd",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::lfo_speed_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::lfo_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_stages,
  n_stereo_busses,
  param_common (
    "Stages",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::stages_tag>()),
  witti::bbd_echo_stereo::get_parameter (witti::bbd_echo_stereo::stages_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_delay,
  n_stereo_busses,
  param_common (
    "Time",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::delay_tag>()),
  witti::bbd_echo_stereo::get_parameter (witti::bbd_echo_stereo::delay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_delay_sync_l,
  n_stereo_busses,
  param_common (
    "SynTime L",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::delay_sync_l_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::delay_sync_l_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_delay_sync_r,
  n_stereo_busses,
  param_common (
    "SynTime R",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::delay_sync_r_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::delay_sync_r_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_hp_filter,
  n_stereo_busses,
  param_common (
    "HP Freq",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::hp_filter_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::hp_filter_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_hp_res,
  n_stereo_busses,
  param_common (
    "HP Res",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::hp_res_tag>()),
  witti::bbd_echo_stereo::get_parameter (witti::bbd_echo_stereo::hp_res_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_lp_filter,
  n_stereo_busses,
  param_common (
    "LP Freq",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::lp_filter_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::lp_filter_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_lp_res,
  n_stereo_busses,
  param_common (
    "LP Res",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::lp_res_tag>()),
  witti::bbd_echo_stereo::get_parameter (witti::bbd_echo_stereo::lp_res_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_clock_offset,
  n_stereo_busses,
  param_common (
    "Clk Offset",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::clock_offset_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::clock_offset_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_clock_scale,
  n_stereo_busses,
  param_common (
    "Clk Scale",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::clock_scale_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::clock_scale_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_clock_curve,
  n_stereo_busses,
  param_common (
    "Clk Curve",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::clock_curve_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::clock_curve_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_hiss,
  n_stereo_busses,
  param_common (
    "Hiss",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::hiss_tag>()),
  witti::bbd_echo_stereo::get_parameter (witti::bbd_echo_stereo::hiss_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_age,
  n_stereo_busses,
  param_common (
    "LR-Diff",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::decalibration_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::decalibration_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::ducking_speed_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  bbd_echo_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<witti::bbd_echo_stereo>(),
    declptr<witti::bbd_echo_stereo::ducking_threshold_tag>()),
  witti::bbd_echo_stereo::get_parameter (
    witti::bbd_echo_stereo::ducking_threshold_tag {}),
  slider_ext);

using bbd_echo_params = mp_list<
// bbd_echo_dry_wet,
#if 0
  bbd_echo_out,
#endif
  bbd_echo_stages,
  bbd_echo_feedback,
  bbd_echo_delay_sync_l,
  bbd_echo_delay_sync_r,

  bbd_echo_delay,
  bbd_echo_clock_scale,
  bbd_echo_clock_curve,
  bbd_echo_clock_offset,

  bbd_echo_lfo_depth,
  bbd_echo_lfo_speed,
  bbd_echo_hp_filter,
  bbd_echo_hp_res,

  bbd_echo_lp_filter,
  bbd_echo_lp_res,
  bbd_echo_hiss,
  bbd_echo_age,

  bbd_echo_ducking_threshold,
  bbd_echo_ducking_speed>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  df_er_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::dry_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_er_wet,
  n_stereo_busses,
  param_common (
    "Wet",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::wet_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::wet_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_er_size,
  n_stereo_busses,
  param_common (
    "Size",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::size_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::size_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_er_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::width_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_er_low_cut,
  n_stereo_busses,
  param_common (
    "Low Cut",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::low_cut_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::low_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_er_high_cut,
  n_stereo_busses,
  param_common (
    "High Cut",
    declptr<dragonfly::early_reflections>(),
    declptr<dragonfly::early_reflections::high_cut_tag>()),
  dragonfly::early_reflections::get_parameter (
    dragonfly::early_reflections::high_cut_tag {}),
  slider_ext);

using df_er_params = mp_list<
  // df_er_dry,
  // df_er_wet,
  df_er_size,
  df_er_width,
  df_er_low_cut,
  df_er_high_cut>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  df_plate_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::dry_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_wet,
  n_stereo_busses,
  param_common (
    "Wet",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::wet_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::wet_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::width_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_low_cut,
  n_stereo_busses,
  param_common (
    "Low Cut",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::low_cut_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::low_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_high_cut,
  n_stereo_busses,
  param_common (
    "High Cut",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::high_cut_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::high_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_algorithm,
  n_stereo_busses,
  param_common (
    "Algo",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::algorithm_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::algorithm_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_predelay,
  n_stereo_busses,
  param_common (
    "Predelay",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::predelay_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_decay,
  n_stereo_busses,
  param_common (
    "Decay",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::decay_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_early_damp,
  n_stereo_busses,
  param_common (
    "Damp",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::early_damp_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::early_damp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::ducking_speed_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_plate_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<dragonfly::plate>(),
    declptr<dragonfly::plate::ducking_threshold_tag>()),
  dragonfly::plate::get_parameter (dragonfly::plate::ducking_threshold_tag {}),
  slider_ext);

using df_plate_params = mp_list<
  // df_plate_dry,
  // df_plate_wet,
  df_plate_algorithm,
  df_plate_predelay,
  df_plate_decay,
  df_plate_early_damp,
  df_plate_low_cut,
  df_plate_high_cut,
  df_plate_ducking_threshold,
  df_plate_ducking_speed,
  df_plate_width>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  df_hall_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::dry_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_early,
  n_stereo_busses,
  param_common (
    "Early",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::early_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::early_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_early_send,
  n_stereo_busses,
  param_common (
    "Early2Late",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::early_send_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::early_send_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_late,
  n_stereo_busses,
  param_common (
    "Late",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::late_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::late_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_size,
  n_stereo_busses,
  param_common (
    "Size",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::size_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::size_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::width_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_predelay,
  n_stereo_busses,
  param_common (
    "Predelay",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::predelay_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_decay,
  n_stereo_busses,
  param_common (
    "Decay",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::decay_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_diffuse,
  n_stereo_busses,
  param_common (
    "Diffuse",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::diffuse_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::diffuse_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_low_cut,
  n_stereo_busses,
  param_common (
    "Low Cut",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::low_cut_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::low_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_high_cut,
  n_stereo_busses,
  param_common (
    "High Cut",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::high_cut_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::high_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_low_x_over,
  n_stereo_busses,
  param_common (
    "Low XOver",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::low_x_over_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::low_x_over_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_high_x_over,
  n_stereo_busses,
  param_common (
    "High XOver",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::high_x_over_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::high_x_over_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_low_mult,
  n_stereo_busses,
  param_common (
    "Low Mult",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::low_mult_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::low_mult_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_high_mult,
  n_stereo_busses,
  param_common (
    "High Mult",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::high_mult_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::high_mult_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_spin,
  n_stereo_busses,
  param_common (
    "Spin",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::spin_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::spin_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_wander,
  n_stereo_busses,
  param_common (
    "Wander",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::wander_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::wander_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_modulation,
  n_stereo_busses,
  param_common (
    "Mod",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::modulation_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::modulation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::ducking_speed_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_hall_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<dragonfly::hall>(),
    declptr<dragonfly::hall::ducking_threshold_tag>()),
  dragonfly::hall::get_parameter (dragonfly::hall::ducking_threshold_tag {}),
  slider_ext);

using df_hall_params = mp_list<
  // df_hall_dry,
  df_hall_early,
  df_hall_early_send,
  df_hall_late,
  df_hall_predelay,

  df_hall_decay,
  df_hall_size,
  df_hall_diffuse,
  df_hall_spin,

  df_hall_wander,
  df_hall_modulation,
  df_hall_low_cut,
  df_hall_high_cut,
  df_hall_low_x_over,
  df_hall_high_x_over,
  df_hall_low_mult,
  df_hall_high_mult,

  df_hall_ducking_threshold,
  df_hall_ducking_speed,
  df_hall_width>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  df_room_dry,
  n_stereo_busses,
  param_common (
    "Dry",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::dry_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_early,
  n_stereo_busses,
  param_common (
    "Early",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::early_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::early_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_early_send,
  n_stereo_busses,
  param_common (
    "Early2Late",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::early_send_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::early_send_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_late,
  n_stereo_busses,
  param_common (
    "Late",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::late_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::late_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_size,
  n_stereo_busses,
  param_common (
    "Size",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::size_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::size_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::width_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_predelay,
  n_stereo_busses,
  param_common (
    "Predelay",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::predelay_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_decay,
  n_stereo_busses,
  param_common (
    "Decay",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::decay_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_diffuse,
  n_stereo_busses,
  param_common (
    "Diffuse",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::diffuse_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::diffuse_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_spin,
  n_stereo_busses,
  param_common (
    "Spin",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::spin_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::spin_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_wander,
  n_stereo_busses,
  param_common (
    "Wander",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::wander_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::wander_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_low_cut,
  n_stereo_busses,
  param_common (
    "Low Cut",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::low_cut_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::low_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_high_cut,
  n_stereo_busses,
  param_common (
    "High Cut",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::high_cut_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::high_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_low_boost,
  n_stereo_busses,
  param_common (
    "Low Boost",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::low_boost_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::low_boost_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_low_boost_freq,
  n_stereo_busses,
  param_common (
    "Low B Freq.",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::low_boost_freq_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::low_boost_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_early_damp,
  n_stereo_busses,
  param_common (
    "Early Damp",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::early_damp_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::early_damp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_late_damp,
  n_stereo_busses,
  param_common (
    "Late Damp",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::late_damp_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::late_damp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::ducking_speed_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  df_room_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<dragonfly::room>(),
    declptr<dragonfly::room::ducking_threshold_tag>()),
  dragonfly::room::get_parameter (dragonfly::room::ducking_threshold_tag {}),
  slider_ext);

using df_room_params = mp_list<
  // df_room_dry,
  df_room_early,
  df_room_early_send,
  df_room_late,
  df_room_predelay,

  df_room_decay,
  df_room_size,
  df_room_early_damp,
  df_room_late_damp,

  df_room_low_cut,
  df_room_high_cut,
  df_room_low_boost_freq,
  df_room_low_boost,

  df_room_diffuse,
  df_room_spin,
  df_room_wander,
  df_room_width,

  df_room_ducking_threshold,
  df_room_ducking_speed>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  _4x4_low_drive,
  n_stereo_busses,
  param_common (
    "L< Drive",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::low_drive_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::low_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_mid_drive,
  n_stereo_busses,
  param_common (
    "M Drive",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::mid_drive_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::mid_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_high_drive,
  n_stereo_busses,
  param_common (
    "H Drive",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::high_drive_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::high_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_low_gain,
  n_stereo_busses,
  param_common (
    "L Gain",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::low_gain_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::low_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_mid_gain,
  n_stereo_busses,
  param_common (
    "M Gain",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::mid_gain_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::mid_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_high_gain,
  n_stereo_busses,
  param_common (
    "H Gain",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::high_gain_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::high_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_low_mid_freq,
  n_stereo_busses,
  param_common (
    "L/M Freq",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::low_mid_freq_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::low_mid_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _4x4_mid_high_freq,
  n_stereo_busses,
  param_common (
    "M/H Freq",
    declptr<sstillwell::_4x4>(),
    declptr<sstillwell::_4x4::mid_high_freq_tag>()),
  sstillwell::_4x4::get_parameter (sstillwell::_4x4::mid_high_freq_tag {}),
  slider_ext);

using _4x4_params = mp_list<
  _4x4_low_mid_freq,
  _4x4_mid_high_freq,
  _4x4_low_drive,
  _4x4_low_gain,
  _4x4_mid_drive,
  _4x4_mid_gain,
  _4x4_high_drive,
  _4x4_high_gain>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  _1175_threshold,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<sstillwell::_1175::threshold_tag>()),
  sstillwell::_1175::get_parameter (sstillwell::_1175::threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _1175_ratio,
  n_stereo_busses,
  param_common (
    "Ratio",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<sstillwell::_1175::ratio_tag>()),
  sstillwell::_1175::get_parameter (sstillwell::_1175::ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _1175_attack,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<sstillwell::_1175::attack_tag>()),
  sstillwell::_1175::get_parameter (sstillwell::_1175::attack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _1175_release,
  n_stereo_busses,
  param_common (
    "Release",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<sstillwell::_1175::release_tag>()),
  sstillwell::_1175::get_parameter (sstillwell::_1175::release_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _1175_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<sstillwell::_1175::gain_tag>()),
  sstillwell::_1175::get_parameter (sstillwell::_1175::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  _1175_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<sstillwell::_1175>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<sstillwell::_1175>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using _1175_params = mp_list<
  _1175_threshold,
  _1175_ratio,
  _1175_attack,
  _1175_release,
  _1175_gain,
  _1175_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  fairly_childish_threshold,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::threshold_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_bias,
  n_stereo_busses,
  param_common (
    "Bias",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::bias_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::bias_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_agc_range,
  n_stereo_busses,
  param_common (
    "AGC Range",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::agc_range_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::agc_range_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_makeup,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::makeup_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::makeup_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_time_constant,
  n_stereo_busses,
  param_common (
    "Time K",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::time_constant_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::time_constant_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_rms_window,
  n_stereo_busses,
  param_common (
    "RMS Win",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<sstillwell::fairly_childish::rms_window_tag>()),
  sstillwell::fairly_childish::get_parameter (
    sstillwell::fairly_childish::rms_window_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fairly_childish_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<sstillwell::fairly_childish>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<sstillwell::fairly_childish>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using fairly_childish_params = mp_list<
  fairly_childish_threshold,
  fairly_childish_bias,
  fairly_childish_time_constant,
  fairly_childish_rms_window,
  fairly_childish_makeup,
  fairly_childish_agc_range,
  fairly_childish_oversampling>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  huge_booty_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<sstillwell::huge_booty>(),
    declptr<sstillwell::huge_booty::drive_tag>()),
  sstillwell::huge_booty::get_parameter (sstillwell::huge_booty::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  huge_booty_frequency,
  n_stereo_busses,
  param_common (
    "Freq",
    declptr<sstillwell::huge_booty>(),
    declptr<sstillwell::huge_booty::frequency_tag>()),
  sstillwell::huge_booty::get_parameter (
    sstillwell::huge_booty::frequency_tag {}),
  slider_ext);

parameter_cpp_class_define (
  huge_booty_mix,
  n_stereo_busses,
  param_common (
    "Mix",
    declptr<sstillwell::huge_booty>(),
    declptr<sstillwell::huge_booty::mix_tag>()),
  sstillwell::huge_booty::get_parameter (sstillwell::huge_booty::mix_tag {}),
  slider_ext);

using hugebooty_params = mp_list<huge_booty_drive, huge_booty_frequency>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  major_tom_threshold,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::threshold_tag>()),
  sstillwell::major_tom::get_parameter (
    sstillwell::major_tom::threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_ratio,
  n_stereo_busses,
  param_common (
    "Ratio",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::ratio_tag>()),
  sstillwell::major_tom::get_parameter (sstillwell::major_tom::ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::gain_tag>()),
  sstillwell::major_tom::get_parameter (sstillwell::major_tom::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_knee,
  n_stereo_busses,
  param_common (
    "Knee",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::knee_tag>()),
  sstillwell::major_tom::get_parameter (sstillwell::major_tom::knee_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_agc,
  n_stereo_busses,
  param_common (
    "AGC",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::agc_tag>()),
  sstillwell::major_tom::get_parameter (sstillwell::major_tom::agc_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_detection,
  n_stereo_busses,
  param_common (
    "Detector",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::detection_tag>()),
  sstillwell::major_tom::get_parameter (
    sstillwell::major_tom::detection_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_detection_src,
  n_stereo_busses,
  param_common (
    "Detection",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<sstillwell::major_tom::detection_src_tag>()),
  sstillwell::major_tom::get_parameter (
    sstillwell::major_tom::detection_src_tag {}),
  slider_ext);

parameter_cpp_class_define (
  major_tom_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<sstillwell::major_tom>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<sstillwell::major_tom>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using major_tom_params = mp_list<
  major_tom_threshold,
  major_tom_ratio,
  major_tom_gain,
  major_tom_knee,
  major_tom_agc,
  major_tom_detection,
  major_tom_detection_src,
  major_tom_oversampling>;

//------------------------------------------------------------------------------
#if 0
parameter_cpp_class_define (
  master_tom_threshold,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::threshold_tag>()),
  sstillwell::master_tom::get_parameter (
    sstillwell::master_tom::threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_ratio,
  n_stereo_busses,
  param_common (
    "Ratio",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::ratio_tag>()),
  sstillwell::master_tom::get_parameter (sstillwell::master_tom::ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::gain_tag>()),
  sstillwell::master_tom::get_parameter (sstillwell::master_tom::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_knee,
  n_stereo_busses,
  param_common (
    "Knee",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::knee_tag>()),
  sstillwell::master_tom::get_parameter (sstillwell::master_tom::knee_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_agc,
  n_stereo_busses,
  param_common (
    "AGC",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::agc_tag>()),
  sstillwell::master_tom::get_parameter (sstillwell::master_tom::agc_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_detection,
  n_stereo_busses,
  param_common (
    "Detector",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::detection_tag>()),
  sstillwell::master_tom::get_parameter (
    sstillwell::master_tom::detection_tag {}),
  slider_ext);

parameter_cpp_class_define (
  master_tom_detection_src,
  n_stereo_busses,
  param_common (
    "Detection",
    declptr<sstillwell::master_tom>(),
    declptr<sstillwell::master_tom::detection_src_tag>()),
  sstillwell::master_tom::get_parameter (
    sstillwell::master_tom::detection_src_tag {}),
  slider_ext);

using master_tom_params = mp_list<
  master_tom_threshold,
  master_tom_ratio,
  master_tom_gain,
  master_tom_knee,
  master_tom_agc,
  master_tom_detection,
  master_tom_detection_src>;
#endif
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  event_horizon_2_threshold,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<sstillwell::event_horizon_2>>(),
    declptr<sstillwell::event_horizon_2::threshold_tag>()),
  sstillwell::event_horizon_2::get_parameter (
    sstillwell::event_horizon_2::threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  event_horizon_2_ceiling,
  n_stereo_busses,
  param_common (
    "Ceiling",
    declptr<updownsampled<sstillwell::event_horizon_2>>(),
    declptr<sstillwell::event_horizon_2::ceiling_tag>()),
  sstillwell::event_horizon_2::get_parameter (
    sstillwell::event_horizon_2::ceiling_tag {}),
  slider_ext);

parameter_cpp_class_define (
  event_horizon_2_release,
  n_stereo_busses,
  param_common (
    "Release",
    declptr<updownsampled<sstillwell::event_horizon_2>>(),
    declptr<sstillwell::event_horizon_2::release_tag>()),
  sstillwell::event_horizon_2::get_parameter (
    sstillwell::event_horizon_2::release_tag {}),
  slider_ext);

parameter_cpp_class_define (
  event_horizon_2_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<sstillwell::event_horizon_2>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<sstillwell::event_horizon_2>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using event_horizon_2_params = mp11::mp_list<
  event_horizon_2_threshold,
  event_horizon_2_ceiling,
  event_horizon_2_release,
  event_horizon_2_oversampling>;
//----------------------------------------------------------------------------
parameter_cpp_class_define (
  rbj1073_hpf,
  n_stereo_busses,
  param_common (
    "HPF",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::hpf_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::hpf_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rbj1073_low_shelf,
  n_stereo_busses,
  param_common (
    "L Shelf",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::low_shelf_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::low_shelf_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rbj1073_low_gain,
  n_stereo_busses,
  param_common (
    "L Gain",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::low_gain_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::low_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rbj1073_mid_freq,
  n_stereo_busses,
  param_common (
    "M Freq",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::mid_freq_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::mid_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rbj1073_mid_gain,
  n_stereo_busses,
  param_common (
    "M Gain",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::mid_gain_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::mid_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rbj1073_high_12k_gain,
  n_stereo_busses,
  param_common (
    "12kHz Gain",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::high_12k_gain_tag>()),
  sstillwell::rbj1073::get_parameter (
    sstillwell::rbj1073::high_12k_gain_tag {}),
  slider_ext);

// unused, post FX
parameter_cpp_class_define (
  rbj1073_gain,
  n_stereo_busses,
  param_common (
    "12kHz Gain",
    declptr<sstillwell::rbj1073>(),
    declptr<sstillwell::rbj1073::gain_tag>()),
  sstillwell::rbj1073::get_parameter (sstillwell::rbj1073::gain_tag {}),
  slider_ext);

using rbj1073_params = mp_list<
  rbj1073_low_shelf,
  rbj1073_low_gain,
  rbj1073_mid_freq,
  rbj1073_mid_gain,
  rbj1073_high_12k_gain,
  rbj1073_hpf>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  consolidator_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::operation_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::operation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_ingain,
  n_stereo_busses,
  param_common (
    "In Gain",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::dbgain_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::dbgain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_sidechain_freq,
  n_stereo_busses,
  param_common (
    "SC HP",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::scfreq_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::scfreq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_trim,
  n_stereo_busses,
  param_common (
    "Out Gain",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::dbtrim_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::dbtrim_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_st_mode,
  n_stereo_busses,
  param_common (
    "St Mode",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::midside_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::midside_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_channel_link,
  n_stereo_busses,
  param_common (
    "St Mode",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<chokehold::consolidator::linkamount_tag>()),
  chokehold::consolidator::get_parameter (
    chokehold::consolidator::linkamount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  consolidator_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<chokehold::consolidator>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<chokehold::consolidator>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using consolidator_params = mp_list<
  consolidator_mode,
  consolidator_ingain,
  consolidator_sidechain_freq,
  consolidator_trim,
  consolidator_channel_link,
  consolidator_st_mode,
  consolidator_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  gate_expander_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::operation_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::operation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_gatethresh,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::gatethresh_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::gatethresh_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_range,
  n_stereo_busses,
  param_common (
    "Range",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::gaterange_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::gaterange_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_hysteresis,
  n_stereo_busses,
  param_common (
    "Hyst",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::gatehyst_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::gatehyst_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_attack,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::gateattack_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::gateattack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_release,
  n_stereo_busses,
  param_common (
    "Release",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::gaterelease_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::gaterelease_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_sidechain_freq,
  n_stereo_busses,
  param_common (
    "SC Freq",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::scfreq_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::scfreq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_channel_link,
  n_stereo_busses,
  param_common (
    "Chnl Link",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::linkamount_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::linkamount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_channel_mode,
  n_stereo_busses,
  param_common (
    "Chnl Mod",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<chokehold::gate_expander::routing_tag>()),
  chokehold::gate_expander::get_parameter (
    chokehold::gate_expander::routing_tag {}),
  slider_ext);

parameter_cpp_class_define (
  gate_expander_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<chokehold::gate_expander>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<chokehold::gate_expander>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using gate_expander_params = mp_list<
  gate_expander_gatethresh,
  gate_expander_hysteresis,
  gate_expander_attack,
  gate_expander_release,
  gate_expander_mode,
  gate_expander_range,
  gate_expander_sidechain_freq,
  gate_expander_channel_link,
  gate_expander_channel_mode,
  gate_expander_oversampling>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  track_comp_dbgain,
  n_stereo_busses,
  param_common (
    "In Gain",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::dbgain_tag>()),
  chokehold::track_comp::get_parameter (chokehold::track_comp::dbgain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compfeedbk,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compfeedbk_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::compfeedbk_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compwindow,
  n_stereo_busses,
  param_common (
    "FB RMS",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compwindow_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::compwindow_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compthresh,
  n_stereo_busses,
  param_common (
    "Thres",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compthresh_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::compthresh_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compratio,
  n_stereo_busses,
  param_common (
    "Ratio",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compratio_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::compratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compattack,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compattack_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::compattack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_comprelease,
  n_stereo_busses,
  param_common (
    "Release",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::comprelease_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::comprelease_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_compknee,
  n_stereo_busses,
  param_common (
    "Knee",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::compknee_tag>()),
  chokehold::track_comp::get_parameter (chokehold::track_comp::compknee_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_comprange,
  n_stereo_busses,
  param_common (
    "Range",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::comprange_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::comprange_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_linkamount,
  n_stereo_busses,
  param_common (
    "St Link",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::linkamount_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::linkamount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_scfreq,
  n_stereo_busses,
  param_common (
    "Sidech HP",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::scfreq_tag>()),
  chokehold::track_comp::get_parameter (chokehold::track_comp::scfreq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_autogain,
  n_stereo_busses,
  param_common (
    "AGC",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::autogain_tag>()),
  chokehold::track_comp::get_parameter (chokehold::track_comp::autogain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_saturation,
  n_stereo_busses,
  param_common (
    "Sat",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::saturation_tag>()),
  chokehold::track_comp::get_parameter (
    chokehold::track_comp::saturation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_dbtrim,
  n_stereo_busses,
  param_common (
    "Out Gain",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<chokehold::track_comp::dbtrim_tag>()),
  chokehold::track_comp::get_parameter (chokehold::track_comp::dbtrim_tag {}),
  slider_ext);

parameter_cpp_class_define (
  track_comp_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<chokehold::track_comp>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<chokehold::track_comp>::get_parameter (
    oversampled_amount_tag {}),
  slider_ext);

using track_comp_params = mp11::mp_list<
  track_comp_compthresh,
  track_comp_compratio,
  track_comp_compattack,
  track_comp_comprelease,
  track_comp_compknee,
  track_comp_scfreq,
  track_comp_compfeedbk,
  track_comp_compwindow,
  track_comp_dbgain,
  track_comp_dbtrim,
  track_comp_comprange,
  track_comp_saturation,
  track_comp_linkamount,
  track_comp_autogain,
  track_comp_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  signal_crusher_down,
  n_stereo_busses,
  param_common (
    "DownSmpl",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::down_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::down_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_dnfilt,
  n_stereo_busses,
  param_common (
    "DS Mode",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::dnfilt_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::dnfilt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_up,
  n_stereo_busses,
  param_common (
    "UpSmpl",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::up_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::up_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_upfilt,
  n_stereo_busses,
  param_common (
    "US Mode",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::upfilt_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::upfilt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_ratio,
  n_stereo_busses,
  param_common (
    "SRate Div",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::ratio_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::ratio_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_bits,
  n_stereo_busses,
  param_common (
    "Bits",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::bits_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::bits_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_dither,
  n_stereo_busses,
  param_common (
    "Dither",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::dither_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::dither_tag {}),
  slider_ext);

parameter_cpp_class_define (
  signal_crusher_blank,
  n_stereo_busses,
  param_common (
    "Blanking",
    declptr<chokehold::signal_crusher>(),
    declptr<chokehold::signal_crusher::blank_tag>()),
  chokehold::signal_crusher::get_parameter (
    chokehold::signal_crusher::blank_tag {}),
  slider_ext);

using signal_crusher_params = mp_list<
  signal_crusher_down,
  signal_crusher_dnfilt,
  signal_crusher_up,
  signal_crusher_upfilt,
  signal_crusher_ratio,
  signal_crusher_bits,
  signal_crusher_dither,
  signal_crusher_blank>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  sound_delay_main,
  n_stereo_busses,
  param_common (
    "Main",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_ms_tag>()),
  sound_delay::get_parameter (sound_delay::delay_ms_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_main_samples,
  n_stereo_busses,
  param_common (
    "Main",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_samples_tag>()),
  sound_delay::get_parameter (sound_delay::delay_samples_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_l,
  n_stereo_busses,
  param_common (
    "Left",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_ms_l_offset_tag>()),
  sound_delay::get_parameter (sound_delay::delay_ms_l_offset_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_r,
  n_stereo_busses,
  param_common (
    "Right",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_ms_r_offset_tag>()),
  sound_delay::get_parameter (sound_delay::delay_ms_r_offset_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_l_sync,
  n_stereo_busses,
  param_common (
    "L Sync",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_beats_l_tag>()),
  sound_delay::get_parameter (sound_delay::delay_beats_l_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_r_sync,
  n_stereo_busses,
  param_common (
    "R Sync",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_beats_r_tag>()),
  sound_delay::get_parameter (sound_delay::delay_beats_r_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_mid_ms,
  n_stereo_busses,
  param_common (
    "Mid",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_mid_ms_tag>()),
  sound_delay::get_parameter (sound_delay::delay_mid_ms_tag {}),
  slider_ext);

parameter_cpp_class_define (
  sound_delay_side_ms,
  n_stereo_busses,
  param_common (
    "Side",
    declptr<sound_delay>(),
    declptr<sound_delay::delay_side_ms_tag>()),
  sound_delay::get_parameter (sound_delay::delay_side_ms_tag {}),
  slider_ext);

using sound_delay_params = mp_list<
  sound_delay_main,
  sound_delay_main_samples,
  sound_delay_l,
  sound_delay_r,
  sound_delay_l_sync,
  sound_delay_r_sync,
  sound_delay_mid_ms,
  sound_delay_side_ms>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  zebigchorus3_algo,
  n_stereo_busses,
  param_common (
    "Algo",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_algo_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_algo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_delay,
  n_stereo_busses,
  param_common (
    "Delay",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_delay_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_delay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_dmod,
  n_stereo_busses,
  param_common (
    "Delay Mod",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_dmod_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_dmod_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_dispstatic,
  n_stereo_busses,
  param_common (
    "Disp Static",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_dispstatic_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_dispstatic_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_dispmod,
  n_stereo_busses,
  param_common (
    "Disp Mod",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_dispmod_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_dispmod_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_rate,
  n_stereo_busses,
  param_common (
    "Rate",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_rate_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_rate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_ratedisp,
  n_stereo_busses,
  param_common (
    "Rate Disp",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_ratedisp_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_ratedisp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_grate,
  n_stereo_busses,
  param_common (
    "G Rate",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_grate_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_rate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_gratedisp,
  n_stereo_busses,
  param_common (
    "G R Disper",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_gratedisp_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_ratedisp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_g0,
  n_stereo_busses,
  param_common (
    "G0",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_g0_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_g0_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_g1,
  n_stereo_busses,
  param_common (
    "G1",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_g1_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_g1_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_gain_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zebigchorus3_drywet,
  n_stereo_busses,
  param_common (
    "Dry/Wet",
    declptr<smashed_transistors::ze_big_chorus3>(),
    declptr<smashed_transistors::ze_big_chorus3::sl_drywet_tag>()),
  smashed_transistors::ze_big_chorus3::get_parameter (
    smashed_transistors::ze_big_chorus3::sl_drywet_tag {}),
  slider_ext);

using zebigchorus3_params = mp_list<
  zebigchorus3_algo,
  zebigchorus3_gain,
  zebigchorus3_delay,
  zebigchorus3_dmod,
  zebigchorus3_dispstatic,
  zebigchorus3_dispmod,
  zebigchorus3_rate,
  zebigchorus3_ratedisp,
  zebigchorus3_grate,
  zebigchorus3_gratedisp,
  zebigchorus3_g1,
  zebigchorus3_g0>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  zelittlechorus_lfoa,
  n_stereo_busses,
  param_common (
    "LFO A",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_lfoa_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_lfoa_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_lfob,
  n_stereo_busses,
  param_common (
    "LFO B",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_lfob_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_lfob_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_lfomix,
  n_stereo_busses,
  param_common (
    "LFO AB Mix",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_lfom_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_lfom_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_depth,
  n_stereo_busses,
  param_common (
    "Mod Depth",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_depth_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_f1,
  n_stereo_busses,
  param_common (
    "F1",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_f1_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_f1_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_f0,
  n_stereo_busses,
  param_common (
    "F0",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_f0_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_f0_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_r,
  n_stereo_busses,
  param_common (
    "R",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_r_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_r_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_lfof,
  n_stereo_busses,
  param_common (
    "LFO F Rate",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_lfof_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_lfof_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_fb,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_fb_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_fb_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_fbtype,
  n_stereo_busses,
  param_common (
    "Fb Type",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_fbtype_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_fbtype_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_fblp,
  n_stereo_busses,
  param_common (
    "Fb Lowpass",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_fblp_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_fblp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_type,
  n_stereo_busses,
  param_common (
    "Type",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_type_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  zelittlechorus_selfpm,
  n_stereo_busses,
  param_common (
    "Self PM",
    declptr<smashed_transistors::ze_little_scanner_chorus>(),
    declptr<smashed_transistors::ze_little_scanner_chorus::sl_selfpm_tag>()),
  smashed_transistors::ze_little_scanner_chorus::get_parameter (
    smashed_transistors::ze_little_scanner_chorus::sl_selfpm_tag {}),
  slider_ext);

using zelittlechorus_params = mp_list<
  zelittlechorus_type,
  zelittlechorus_fbtype,
  zelittlechorus_fb,
  zelittlechorus_fblp,
  zelittlechorus_f1,
  zelittlechorus_f0,
  zelittlechorus_r,
  zelittlechorus_lfof,
  zelittlechorus_lfoa,
  zelittlechorus_lfob,
  zelittlechorus_lfomix,
  zelittlechorus_depth,
  zelittlechorus_selfpm>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  chow_phaser_lfo_depth,
  n_stereo_busses,
  param_common (
    "Depth",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::lfo_depth_tag>()),
  chow::phaser::get_parameter (chow::phaser::lfo_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_lfo_freq,
  n_stereo_busses,
  param_common (
    "Freq",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::lfo_freq_tag>()),
  chow::phaser::get_parameter (chow::phaser::lfo_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::feedback_tag>()),
  chow::phaser::get_parameter (chow::phaser::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_modulation,
  n_stereo_busses,
  param_common (
    "Mod",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::modulation_tag>()),
  chow::phaser::get_parameter (chow::phaser::modulation_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_stages,
  n_stereo_busses,
  param_common (
    "Stages",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::stages_tag>()),
  chow::phaser::get_parameter (chow::phaser::stages_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_freq_mult,
  n_stereo_busses,
  param_common (
    "Freq Mult",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::freq_mult_tag>()),
  chow::phaser::get_parameter (chow::phaser::freq_mult_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_skew,
  n_stereo_busses,
  param_common (
    "Skew",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::skew_tag>()),
  chow::phaser::get_parameter (chow::phaser::skew_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_src_channel,
  n_stereo_busses,
  param_common (
    "Chnl In",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::src_channel_tag>()),
  chow::phaser::get_parameter (chow::phaser::src_channel_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_d1,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::d1_tag>()),
  chow::phaser::get_parameter (chow::phaser::d1_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_d2,
  n_stereo_busses,
  param_common (
    "Trash",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::d2_tag>()),
  chow::phaser::get_parameter (chow::phaser::d2_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_d3,
  n_stereo_busses,
  param_common (
    "Dirt",
    declptr<updownsampled<chow::phaser>>(),
    declptr<chow::phaser::d3_tag>()),
  chow::phaser::get_parameter (chow::phaser::d3_tag {}),
  slider_ext);

parameter_cpp_class_define (
  chow_phaser_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<chow::phaser>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<chow::phaser>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using chow_phaser_params = mp_list<
  chow_phaser_lfo_depth,
  chow_phaser_lfo_freq,
  chow_phaser_feedback,
  chow_phaser_modulation,
  chow_phaser_stages,
  chow_phaser_freq_mult,
  chow_phaser_skew,
  chow_phaser_src_channel,
  chow_phaser_d1,
  chow_phaser_d2,
  chow_phaser_d3,
  chow_phaser_oversampling>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  tal_reverb2_decay,
  n_stereo_busses,
  param_common (
    "Decay",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::decay_time_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::decay_time_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_predelay,
  n_stereo_busses,
  param_common (
    "Predelay",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::predelay_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_predelay_sync,
  n_stereo_busses,
  param_common (
    "PreDly Syn",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::predelay_sync_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::predelay_sync_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_lowshelf_frequency,
  n_stereo_busses,
  param_common (
    "LCut Frq",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::lowshelf_frequency_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::lowshelf_frequency_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_peak_frequency,
  n_stereo_busses,
  param_common (
    "Notch Frq",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::peak_frequency_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::peak_frequency_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_highshelf_frequency,
  n_stereo_busses,
  param_common (
    "HiCut Frq",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::highshelf_frequency_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::highshelf_frequency_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_lowshelf_gain,
  n_stereo_busses,
  param_common (
    "LCut Amt",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::lowshelf_gain_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::lowshelf_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_peak_gain,
  n_stereo_busses,
  param_common (
    "Notch Amt",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::peak_gain_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::peak_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_highshelf_gain,
  n_stereo_busses,
  param_common (
    "HiCut Amt",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::highshelf_gain_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::highshelf_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_stereo_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::stereo_width_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::stereo_width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::ducking_speed_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  tal_reverb2_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<tal::reverb2>(),
    declptr<tal::reverb2::ducking_threshold_tag>()),
  tal::reverb2::get_parameter (tal::reverb2::ducking_threshold_tag {}),
  slider_ext);

using tal_reverb2_params = mp_list<
  tal_reverb2_decay,
  tal_reverb2_stereo_width,
  tal_reverb2_lowshelf_frequency,
  tal_reverb2_lowshelf_gain,
  tal_reverb2_peak_frequency,
  tal_reverb2_peak_gain,
  tal_reverb2_highshelf_frequency,
  tal_reverb2_highshelf_gain,
  tal_reverb2_predelay,
  tal_reverb2_predelay_sync,
  tal_reverb2_ducking_threshold,
  tal_reverb2_ducking_speed>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  luftikus_gain_10hz,
  n_stereo_busses,
  param_common (
    "10Hz",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_10hz_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_10hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_gain_40hz,
  n_stereo_busses,
  param_common (
    "40Hz",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_40hz_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_40hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_gain_160hz,
  n_stereo_busses,
  param_common (
    "160Hz",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_160hz_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_160hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_gain_640hz,
  n_stereo_busses,
  param_common (
    "640Hz",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_640hz_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_640hz_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_gain_2k5,
  n_stereo_busses,
  param_common (
    "2.5kHz",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_2k5_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_2k5_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_gain_hi,
  n_stereo_busses,
  param_common (
    "HiShelf",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::gain_hi_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::gain_hi_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_freq_hi,
  n_stereo_busses,
  param_common (
    "HiShelf F",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::freq_hi_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::freq_hi_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::mode_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  luftikus_keep_gain,
  n_stereo_busses,
  param_common (
    "Keep Gain",
    declptr<ljkb::luftikus>(),
    declptr<ljkb::luftikus::keep_gain_tag>()),
  ljkb::luftikus::get_parameter (ljkb::luftikus::keep_gain_tag {}),
  slider_ext);

using luftikus_params = mp_list<
  luftikus_gain_10hz,
  luftikus_gain_40hz,
  luftikus_gain_160hz,
  luftikus_gain_640hz,
  luftikus_gain_2k5,
  luftikus_gain_hi,
  luftikus_freq_hi,
  luftikus_mode,
  luftikus_keep_gain>;

//------------------------------------------------------------------------------
parameter_cpp_class_define (
  fdnverb_time,
  n_stereo_busses,
  param_common (
    "Time",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::time_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::time_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::feedback_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_density,
  n_stereo_busses,
  param_common (
    "Density",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::density_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::density_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_cascade,
  n_stereo_busses,
  param_common (
    "Cascade",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::cascade_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::cascade_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_mod_rate,
  n_stereo_busses,
  param_common (
    "Mod Rate",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::mod_rate_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::mod_rate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_mod_depth,
  n_stereo_busses,
  param_common (
    "Mod Depth",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::mod_depth_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::mod_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_stereo,
  n_stereo_busses,
  param_common (
    "Stereo",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::stereo_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_smooth,
  n_stereo_busses,
  param_common (
    "Smooth",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::smooth_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::smooth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_hipass,
  n_stereo_busses,
  param_common (
    "Hipass",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::hipass_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::hipass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_lopass,
  n_stereo_busses,
  param_common (
    "Lopass",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::lopass_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::lopass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_pre_shift,
  n_stereo_busses,
  param_common (
    "Pre Shift",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::pre_shift_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::pre_shift_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_post_shift,
  n_stereo_busses,
  param_common (
    "Post Shift",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::post_shift_tag>()),
  shabtronic::fdn_verb::get_parameter (shabtronic::fdn_verb::post_shift_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::ducking_speed_tag>()),
  shabtronic::fdn_verb::get_parameter (
    shabtronic::fdn_verb::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  fdnverb_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<shabtronic::fdn_verb>(),
    declptr<shabtronic::fdn_verb::ducking_threshold_tag>()),
  shabtronic::fdn_verb::get_parameter (
    shabtronic::fdn_verb::ducking_threshold_tag {}),
  slider_ext);

using fdnverb_params = mp_list<
  fdnverb_time,
  fdnverb_feedback,
  fdnverb_density,
  fdnverb_cascade,

  fdnverb_mod_rate,
  fdnverb_mod_depth,
  fdnverb_pre_shift,
  fdnverb_post_shift,

  fdnverb_stereo,
  fdnverb_smooth,
  fdnverb_hipass,
  fdnverb_lopass,

  fdnverb_ducking_threshold,
  fdnverb_ducking_speed>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  myphaser_stages,
  n_stereo_busses,
  param_common (
    "N Stages",
    declptr<upsampled<phaser>>(),
    declptr<phaser::stages_tag>()),
  phaser::get_parameter (phaser::stages_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_rate,
  n_stereo_busses,
  param_common (
    "LFO Freq",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_rate_tag>()),
  phaser::get_parameter (phaser::lfo_rate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_rate_sync,
  n_stereo_busses,
  param_common (
    "LFO Sync",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_rate_sync_tag>()),
  phaser::get_parameter (phaser::lfo_rate_sync_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_depth,
  n_stereo_busses,
  param_common (
    "LFO Amt",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_depth_tag>()),
  phaser::get_parameter (phaser::lfo_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_start_phase,
  n_stereo_busses,
  param_common (
    "LFO Phase",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_start_phase_tag>()),
  phaser::get_parameter (phaser::lfo_start_phase_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_stereo,
  n_stereo_busses,
  param_common (
    "LFO Ster",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_stereo_tag>()),
  phaser::get_parameter (phaser::lfo_stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_lfo_wave,
  n_stereo_busses,
  param_common (
    "LFO Wave",
    declptr<upsampled<phaser>>(),
    declptr<phaser::lfo_wave_tag>()),
  phaser::get_parameter (phaser::lfo_wave_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_low_freq,
  n_stereo_busses,
  param_common (
    "Low Freq/T",
    declptr<upsampled<phaser>>(),
    declptr<phaser::low_freq_tag>()),
  phaser::get_parameter (phaser::low_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_high_freq,
  n_stereo_busses,
  param_common (
    "High Freq/T",
    declptr<upsampled<phaser>>(),
    declptr<phaser::high_freq_tag>()),
  phaser::get_parameter (phaser::high_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<upsampled<phaser>>(),
    declptr<phaser::feedback_tag>()),
  phaser::get_parameter (phaser::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_feedback_hp,
  n_stereo_busses,
  param_common (
    "Feedbk HP",
    declptr<upsampled<phaser>>(),
    declptr<phaser::feedback_hp_tag>()),
  phaser::get_parameter (phaser::feedback_hp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_feedback_lp,
  n_stereo_busses,
  param_common (
    "Feedbk LP",
    declptr<upsampled<phaser>>(),
    declptr<phaser::feedback_lp_tag>()),
  phaser::get_parameter (phaser::feedback_lp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_feedback_sat,
  n_stereo_busses,
  param_common (
    "Feedbk Sat",
    declptr<upsampled<phaser>>(),
    declptr<phaser::feedback_sat_tag>()),
  phaser::get_parameter (phaser::feedback_sat_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_q,
  n_stereo_busses,
  param_common (
    "Q/Special",
    declptr<upsampled<phaser>>(),
    declptr<phaser::q_tag>()),
  phaser::get_parameter (phaser::q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_stages_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<upsampled<phaser>>(),
    declptr<phaser::stages_mode_tag>()),
  phaser::get_parameter (phaser::stages_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_topology,
  n_stereo_busses,
  param_common (
    "Topology",
    declptr<upsampled<phaser>>(),
    declptr<phaser::topology_tag>()),
  phaser::get_parameter (phaser::topology_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_parallel_mix,
  n_stereo_busses,
  param_common (
    "Parallel Mix",
    declptr<upsampled<phaser>>(),
    declptr<phaser::parallel_mix_tag>()),
  phaser::get_parameter (phaser::parallel_mix_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_delay_feedback,
  n_stereo_busses,
  param_common (
    "Flan Fback",
    declptr<upsampled<phaser>>(),
    declptr<phaser::delay_feedback_tag>()),
  phaser::get_parameter (phaser::delay_feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_delay_time,
  n_stereo_busses,
  param_common (
    "Flan Del",
    declptr<upsampled<phaser>>(),
    declptr<phaser::delay_time_tag>()),
  phaser::get_parameter (phaser::delay_time_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_delay_lfo,
  n_stereo_busses,
  param_common (
    "Flan LFO",
    declptr<upsampled<phaser>>(),
    declptr<phaser::delay_lfo_tag>()),
  phaser::get_parameter (phaser::delay_lfo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  myphaser_oversampling,
  n_stereo_busses,
  param_common (
    "Upsample",
    declptr<upsampled<phaser>>(),
    declptr<oversampled_amount_tag>()),
  upsampled<phaser>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using myphaser_params = mp_list<
  myphaser_lfo_rate,
  myphaser_lfo_rate_sync,
  myphaser_lfo_depth,
  myphaser_stages,
  myphaser_low_freq,
  myphaser_high_freq,
  myphaser_q,
  myphaser_topology,

  myphaser_feedback,
  myphaser_feedback_sat,
  myphaser_feedback_hp,
  myphaser_feedback_lp,
  myphaser_delay_feedback,
  myphaser_delay_time,
  myphaser_delay_lfo,
  myphaser_parallel_mix,

  myphaser_stages_mode,
  myphaser_lfo_wave,
  myphaser_lfo_stereo,
  myphaser_lfo_start_phase,
  myphaser_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  eq4x_band1_type,
  n_stereo_busses,
  param_common (
    "1.Type",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band1_type_tag>()),
  eq4x::get_parameter (eq4x::band1_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band1_freq,
  n_stereo_busses,
  param_common (
    "1.Freq",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band1_freq_tag>()),
  eq4x::get_parameter (eq4x::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band1_q,
  n_stereo_busses,
  param_common (
    "1.Q/Order",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band1_q_tag>()),
  eq4x::get_parameter (eq4x::band1_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band1_gain,
  n_stereo_busses,
  param_common (
    "1.Gain",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band1_gain_tag>()),
  eq4x::get_parameter (eq4x::band1_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band2_type,
  n_stereo_busses,
  param_common (
    "2.Type",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band2_type_tag>()),
  eq4x::get_parameter (eq4x::band2_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band2_freq,
  n_stereo_busses,
  param_common (
    "2.Freq",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band2_freq_tag>()),
  eq4x::get_parameter (eq4x::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band2_q,
  n_stereo_busses,
  param_common (
    "2.Q/Order",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band2_q_tag>()),
  eq4x::get_parameter (eq4x::band2_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band2_gain,
  n_stereo_busses,
  param_common (
    "2.Gain",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band2_gain_tag>()),
  eq4x::get_parameter (eq4x::band2_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band3_type,
  n_stereo_busses,
  param_common (
    "3.Type",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band3_type_tag>()),
  eq4x::get_parameter (eq4x::band3_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band3_freq,
  n_stereo_busses,
  param_common (
    "3.Freq",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band3_freq_tag>()),
  eq4x::get_parameter (eq4x::band3_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band3_q,
  n_stereo_busses,
  param_common (
    "3.Q/Order",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band3_q_tag>()),
  eq4x::get_parameter (eq4x::band3_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band3_gain,
  n_stereo_busses,
  param_common (
    "3.Gain",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band3_gain_tag>()),
  eq4x::get_parameter (eq4x::band3_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_type,
  n_stereo_busses,
  param_common (
    "4.Type",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band4_type_tag>()),
  eq4x::get_parameter (eq4x::band4_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_freq,
  n_stereo_busses,
  param_common (
    "4.Freq",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band4_freq_tag>()),
  eq4x::get_parameter (eq4x::band4_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_q,
  n_stereo_busses,
  param_common (
    "4.Q/Order",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band4_q_tag>()),
  eq4x::get_parameter (eq4x::band4_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_gain,
  n_stereo_busses,
  param_common (
    "4.Gain",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band4_gain_tag>()),
  eq4x::get_parameter (eq4x::band4_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band1_tolerance,
  n_stereo_busses,
  param_common (
    "1.LR-Diff",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band1_diff_tag>()),
  eq4x::get_parameter (eq4x::band1_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band2_tolerance,
  n_stereo_busses,
  param_common (
    "2.LR-Diff",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band2_diff_tag>()),
  eq4x::get_parameter (eq4x::band2_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band3_tolerance,
  n_stereo_busses,
  param_common (
    "3.LR-Diff",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band3_diff_tag>()),
  eq4x::get_parameter (eq4x::band3_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_tolerance,
  n_stereo_busses,
  param_common (
    "4.LR-Diff",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::band4_diff_tag>()),
  eq4x::get_parameter (eq4x::band4_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_band4_topology,
  n_stereo_busses,
  param_common (
    "Topology",
    declptr<upsampled<eq4x>>(),
    declptr<eq4x::topology_tag>()),
  eq4x::get_parameter (eq4x::topology_tag {}),
  slider_ext);

parameter_cpp_class_define (
  eq4x_upsampling,
  n_stereo_busses,
  param_common (
    "Upsample",
    declptr<upsampled<eq4x>>(),
    declptr<oversampled_amount_tag>()),
  upsampled<eq4x>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using eq4x_params = mp_list<
  eq4x_band1_type,
  eq4x_band1_gain,
  eq4x_band1_freq,
  eq4x_band1_q,
  eq4x_band2_type,
  eq4x_band2_gain,
  eq4x_band2_freq,
  eq4x_band2_q,
  eq4x_band3_type,
  eq4x_band3_gain,
  eq4x_band3_freq,
  eq4x_band3_q,
  eq4x_band4_type,
  eq4x_band4_gain,
  eq4x_band4_freq,
  eq4x_band4_q,
  eq4x_band1_tolerance,
  eq4x_band2_tolerance,
  eq4x_band3_tolerance,
  eq4x_band4_tolerance,
  eq4x_band4_topology,
  eq4x_upsampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  lin_eq4x_band1_type,
  n_stereo_busses,
  param_common (
    "1.Type",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band1_type_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band1_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band1_freq,
  n_stereo_busses,
  param_common (
    "1.Freq",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band1_freq_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band1_q,
  n_stereo_busses,
  param_common ("1.Q", declptr<lin_eq4x>(), declptr<lin_eq4x::band1_q_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band1_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band1_gain,
  n_stereo_busses,
  param_common (
    "1.Gain",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band1_gain_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band1_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band2_type,
  n_stereo_busses,
  param_common (
    "2.Type",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band2_type_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band2_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band2_freq,
  n_stereo_busses,
  param_common (
    "2.Freq",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band2_freq_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band2_q,
  n_stereo_busses,
  param_common ("2.Q", declptr<lin_eq4x>(), declptr<lin_eq4x::band2_q_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band2_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band2_gain,
  n_stereo_busses,
  param_common (
    "2.Gain",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band2_gain_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band2_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band3_type,
  n_stereo_busses,
  param_common (
    "3.Type",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band3_type_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band3_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band3_freq,
  n_stereo_busses,
  param_common (
    "3.Freq",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band3_freq_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band3_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band3_q,
  n_stereo_busses,
  param_common ("3.Q", declptr<lin_eq4x>(), declptr<lin_eq4x::band3_q_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band3_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band3_gain,
  n_stereo_busses,
  param_common (
    "3.Gain",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band3_gain_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band3_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band4_type,
  n_stereo_busses,
  param_common (
    "4.Type",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band4_type_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band4_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band4_freq,
  n_stereo_busses,
  param_common (
    "4.Freq",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band4_freq_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band4_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band4_q,
  n_stereo_busses,
  param_common ("4.Q", declptr<lin_eq4x>(), declptr<lin_eq4x::band4_q_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band4_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band4_gain,
  n_stereo_busses,
  param_common (
    "4.Gain",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band4_gain_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band4_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band1_tolerance,
  n_stereo_busses,
  param_common (
    "1.LR-Diff",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band1_diff_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band1_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band2_tolerance,
  n_stereo_busses,
  param_common (
    "2.LR-Diff",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band2_diff_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band2_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band3_tolerance,
  n_stereo_busses,
  param_common (
    "3.LR-Diff",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band3_diff_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band3_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_band4_tolerance,
  n_stereo_busses,
  param_common (
    "4.LR-Diff",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::band4_diff_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::band4_diff_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lin_eq4x_quality,
  n_stereo_busses,
  param_common (
    "Quality",
    declptr<lin_eq4x>(),
    declptr<lin_eq4x::quality_tag>()),
  lin_eq4x::get_parameter (lin_eq4x::quality_tag {}),
  slider_ext);

using lin_eq4x_params = mp_list<
  lin_eq4x_band1_type,
  lin_eq4x_band1_gain,
  lin_eq4x_band1_freq,
  lin_eq4x_band1_q,
  lin_eq4x_band2_type,
  lin_eq4x_band2_gain,
  lin_eq4x_band2_freq,
  lin_eq4x_band2_q,
  lin_eq4x_band3_type,
  lin_eq4x_band3_gain,
  lin_eq4x_band3_freq,
  lin_eq4x_band3_q,
  lin_eq4x_band4_type,
  lin_eq4x_band4_gain,
  lin_eq4x_band4_freq,
  lin_eq4x_band4_q,
  lin_eq4x_band1_tolerance,
  lin_eq4x_band2_tolerance,
  lin_eq4x_band3_tolerance,
  lin_eq4x_band4_tolerance,
  lin_eq4x_quality>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  filter2x_band1_type,
  n_stereo_busses,
  param_common (
    "1.Type",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_type_tag>()),
  filter2x::get_parameter (filter2x::band1_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_topology,
  n_stereo_busses,
  param_common (
    "1.Topo",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_topology_tag>()),
  filter2x::get_parameter (filter2x::band1_topology_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_freq,
  n_stereo_busses,
  param_common (
    "1.Freq",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_freq_tag>()),
  filter2x::get_parameter (filter2x::band1_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_reso,
  n_stereo_busses,
  param_common (
    "1.Reso",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_reso_tag>()),
  filter2x::get_parameter (filter2x::band1_reso_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_drive,
  n_stereo_busses,
  param_common (
    "1.Drive",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_drive_tag>()),
  filter2x::get_parameter (filter2x::band1_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_gain,
  n_stereo_busses,
  param_common (
    "1.Gain",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_gain_tag>()),
  filter2x::get_parameter (filter2x::band1_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_feedback,
  n_stereo_busses,
  param_common (
    "1.Feedbk",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_feedback_tag>()),
  filter2x::get_parameter (filter2x::band1_feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_tolerance,
  n_stereo_busses,
  param_common (
    "1.LR-Diff",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_tolerance_tag>()),
  filter2x::get_parameter (filter2x::band1_tolerance_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_type,
  n_stereo_busses,
  param_common (
    "2.Type",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_type_tag>()),
  filter2x::get_parameter (filter2x::band2_type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_topology,
  n_stereo_busses,
  param_common (
    "2.Topo",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_topology_tag>()),
  filter2x::get_parameter (filter2x::band2_topology_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_freq,
  n_stereo_busses,
  param_common (
    "2.Freq",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_freq_tag>()),
  filter2x::get_parameter (filter2x::band2_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_reso,
  n_stereo_busses,
  param_common (
    "2.Reso",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_reso_tag>()),
  filter2x::get_parameter (filter2x::band2_reso_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_drive,
  n_stereo_busses,
  param_common (
    "2.Drive",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_drive_tag>()),
  filter2x::get_parameter (filter2x::band2_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_gain,
  n_stereo_busses,
  param_common (
    "2.Gain",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_gain_tag>()),
  filter2x::get_parameter (filter2x::band2_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_feedback,
  n_stereo_busses,
  param_common (
    "2.Feedbk",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_feedback_tag>()),
  filter2x::get_parameter (filter2x::band2_feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_tolerance,
  n_stereo_busses,
  param_common (
    "2.LR-Diff",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_tolerance_tag>()),
  filter2x::get_parameter (filter2x::band2_tolerance_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_ef_to_freq,
  n_stereo_busses,
  param_common (
    "1.EF:Freq",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_ef_to_freq_tag>()),
  filter2x::get_parameter (filter2x::band1_ef_to_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_ef_to_freq,
  n_stereo_busses,
  param_common (
    "2.EF:Freq",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_ef_to_freq_tag>()),
  filter2x::get_parameter (filter2x::band2_ef_to_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band1_ef_to_reso,
  n_stereo_busses,
  param_common (
    "1.EF:reso",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band1_ef_to_reso_tag>()),
  filter2x::get_parameter (filter2x::band1_ef_to_reso_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_band2_ef_to_reso,
  n_stereo_busses,
  param_common (
    "2.EF:reso",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::band2_ef_to_reso_tag>()),
  filter2x::get_parameter (filter2x::band2_ef_to_reso_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_envfollow_sensitivity,
  n_stereo_busses,
  param_common (
    "EF Sens",
    declptr<updownsampled<filter2x>>(),
    declptr<filter2x::envfollow_sensitivity_tag>()),
  filter2x::get_parameter (filter2x::envfollow_sensitivity_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_envfollow_attack,
  n_stereo_busses,
  param_common (
    "EF Attack",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::envfollow_attack_tag>()),
  filter2x::get_parameter (filter2x::envfollow_attack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_envfollow_release,
  n_stereo_busses,
  param_common (
    "EF Rel",
    declptr<oversampled<filter2x>>(),
    declptr<filter2x::envfollow_release_tag>()),
  filter2x::get_parameter (filter2x::envfollow_release_tag {}),
  slider_ext);

parameter_cpp_class_define (
  filter2x_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<oversampled<filter2x>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<filter2x>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using filter2x_params = mp_list<
  filter2x_band1_type,
  filter2x_band1_topology,
  filter2x_band1_freq,
  filter2x_band1_reso,
  filter2x_band1_drive,
  filter2x_band1_gain,
  filter2x_band1_feedback,
  filter2x_band1_tolerance,

  filter2x_envfollow_attack,
  filter2x_envfollow_release,
  filter2x_envfollow_sensitivity,
  filter2x_oversampling,
  filter2x_band1_ef_to_freq,
  filter2x_band2_ef_to_freq,
  filter2x_band1_ef_to_reso,
  filter2x_band2_ef_to_reso,

  filter2x_band2_type,
  filter2x_band2_topology,
  filter2x_band2_freq,
  filter2x_band2_reso,
  filter2x_band2_drive,
  filter2x_band2_gain,
  filter2x_band2_feedback,
  filter2x_band2_tolerance>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  waveshaper_band_gain,
  n_stereo_busses,
  param_common (
    "Sat Vol",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::saturated_out_tag>()),
  waveshaper::get_parameter (waveshaper::saturated_out_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_trim,
  n_stereo_busses,
  param_common (
    "In Trim",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::trim_tag>()),
  waveshaper::get_parameter (waveshaper::trim_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_drive,
  n_stereo_busses,
  param_common (
    "Drive",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::drive_tag>()),
  waveshaper::get_parameter (waveshaper::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_drive_balance,
  n_stereo_busses,
  param_common (
    "LR-Diff",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::drive_balance_tag>()),
  waveshaper::get_parameter (waveshaper::drive_balance_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_type,
  n_stereo_busses,
  param_common (
    "Shape",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::type_tag>()),
  waveshaper::get_parameter (waveshaper::type_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::mode_tag>()),
  waveshaper::get_parameter (waveshaper::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_band_lo,
  n_stereo_busses,
  param_common (
    "Sat Lo",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::lo_cut_tag>()),
  waveshaper::get_parameter (waveshaper::lo_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_band_hi,
  n_stereo_busses,
  param_common (
    "Sat Hi",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::hi_cut_tag>()),
  waveshaper::get_parameter (waveshaper::hi_cut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_lo_mode,
  n_stereo_busses,
  param_common (
    "Sat LoCfg",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::lo_mode_tag>()),
  waveshaper::get_parameter (waveshaper::lo_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_hi_mode,
  n_stereo_busses,
  param_common (
    "Sat HiCfg",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::hi_mode_tag>()),
  waveshaper::get_parameter (waveshaper::hi_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_emphasis_amount,
  n_stereo_busses,
  param_common (
    "Tone Amt",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::emphasis_amount_tag>()),
  waveshaper::get_parameter (waveshaper::emphasis_amount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_emphasis_freq,
  n_stereo_busses,
  param_common (
    "Tone Frq",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::emphasis_freq_tag>()),
  waveshaper::get_parameter (waveshaper::emphasis_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_emphasis_q,
  n_stereo_busses,
  param_common (
    "Tone Q",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::emphasis_q_tag>()),
  waveshaper::get_parameter (waveshaper::emphasis_q_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_attack,
  n_stereo_busses,
  param_common (
    "EF Attack",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_attack_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_attack_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_release,
  n_stereo_busses,
  param_common (
    "EF Rel",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_release_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_release_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_sensitivity,
  n_stereo_busses,
  param_common (
    "EF Sens",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_sensitivity_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_sensitivity_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_mode,
  n_stereo_busses,
  param_common (
    "EF Mode",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_mode_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_to_drive,
  n_stereo_busses,
  param_common (
    "EF:Drive",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_to_drive_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_to_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_to_emphasis_freq,
  n_stereo_busses,
  param_common (
    "EF:Tone F",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_to_emphasis_freq_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_to_emphasis_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_to_emphasis_amount,
  n_stereo_busses,
  param_common (
    "EF:Tone A",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_to_emphasis_amount_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_to_emphasis_amount_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_envfollow_to_dc,
  n_stereo_busses,
  param_common (
    "EF:Asym",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::envfollow_to_dc_tag>()),
  waveshaper::get_parameter (waveshaper::envfollow_to_dc_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_feedback,
  n_stereo_busses,
  param_common (
    "Feedbk",
    declptr<updownsampled<waveshaper>>(),
    declptr<waveshaper::feedback_tag>()),
  waveshaper::get_parameter (waveshaper::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  waveshaper_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<waveshaper>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<waveshaper>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using waveshaper_params = mp_list<
  waveshaper_type,
  waveshaper_trim,
  waveshaper_drive,
  waveshaper_drive_balance,
  waveshaper_feedback,
  waveshaper_emphasis_amount,
  waveshaper_emphasis_freq,
  waveshaper_emphasis_q,

  waveshaper_envfollow_attack,
  waveshaper_envfollow_release,
  waveshaper_envfollow_sensitivity,
  waveshaper_envfollow_mode,
  waveshaper_envfollow_to_drive,
  waveshaper_envfollow_to_dc,
  waveshaper_envfollow_to_emphasis_freq,
  waveshaper_envfollow_to_emphasis_amount,

  waveshaper_lo_mode,
  waveshaper_band_lo,
  waveshaper_hi_mode,
  waveshaper_band_hi,
  waveshaper_band_gain,
  waveshaper_mode,
  waveshaper_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  transient_gate_decay,
  n_stereo_busses,
  param_common (
    "Dec Time",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::decay_tag>()),
  transient_gate_fx::get_parameter (transient_gate_fx::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_decay_shape,
  n_stereo_busses,
  param_common (
    "Dec Shape",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::decay_shape_tag>()),
  transient_gate_fx::get_parameter (transient_gate_fx::decay_shape_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_detect_hipass,
  n_stereo_busses,
  param_common (
    "Det Hipass",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::detector_hipass_tag>()),
  transient_gate_fx::get_parameter (transient_gate_fx::detector_hipass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_detect_recovery,
  n_stereo_busses,
  param_common (
    "Det Recov",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::detector_recovery_tag>()),
  transient_gate_fx::get_parameter (
    transient_gate_fx::detector_recovery_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_detect_channels,
  n_stereo_busses,
  param_common (
    "Det Chnl",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::detector_channels_tag>()),
  transient_gate_fx::get_parameter (
    transient_gate_fx::detector_channels_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_detect_shape,
  n_stereo_busses,
  param_common (
    "Det Shape",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::detector_shape_tag>()),
  transient_gate_fx::get_parameter (transient_gate_fx::detector_shape_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_output,
  n_stereo_busses,
  param_common (
    "Output",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<transient_gate_fx::output_tag>()),
  transient_gate_fx::get_parameter (transient_gate_fx::output_tag {}),
  slider_ext);

parameter_cpp_class_define (
  transient_gate_oversampling,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<transient_gate_fx>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<transient_gate_fx>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using transient_gate_params = mp_list<
  transient_gate_decay,
  transient_gate_decay_shape,
  transient_gate_detect_recovery,
  transient_gate_detect_shape,
  transient_gate_detect_hipass,
  transient_gate_detect_channels,
  transient_gate_output,
  transient_gate_oversampling>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  naive_pitch_amt,
  n_stereo_busses,
  param_common (
    "Amt",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<pitch_shifter::semitones_tag>()),
  pitch_shifter::get_parameter (pitch_shifter::semitones_tag {}),
  slider_ext);

parameter_cpp_class_define (
  naive_pitch_count,
  n_stereo_busses,
  param_common (
    "Count",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<pitch_shifter::count_tag>()),
  pitch_shifter::get_parameter (pitch_shifter::count_tag {}),
  slider_ext);

parameter_cpp_class_define (
  naive_pitch_detune,
  n_stereo_busses,
  param_common (
    "Detune",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<pitch_shifter::detune_tag>()),
  pitch_shifter::get_parameter (pitch_shifter::detune_tag {}),
  slider_ext);

parameter_cpp_class_define (
  naive_pitch_width,
  n_stereo_busses,
  param_common (
    "Width",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<pitch_shifter::width_tag>()),
  pitch_shifter::get_parameter (pitch_shifter::width_tag {}),
  slider_ext);

parameter_cpp_class_define (
  naive_pitch_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<pitch_shifter::mode_tag>()),
  pitch_shifter::get_parameter (pitch_shifter::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  naive_pitch_oversample,
  n_stereo_busses,
  param_common (
    "OverSmpl",
    declptr<updownsampled<pitch_shifter>>(),
    declptr<oversampled_amount_tag>()),
  updownsampled<pitch_shifter>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using naive_pitch_params = mp_list<
  naive_pitch_amt,
  naive_pitch_count,
  naive_pitch_detune,
  naive_pitch_width,
  naive_pitch_mode,
  naive_pitch_oversample>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  rubberband_mode,
  n_stereo_busses,
  param_common ("Mode", declptr<rubberband>(), declptr<rubberband::mode_tag>()),
  rubberband::get_parameter (rubberband::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rubberband_formant_mode,
  n_stereo_busses,
  param_common (
    "Formant",
    declptr<rubberband>(),
    declptr<rubberband::formant_mode_tag>()),
  rubberband::get_parameter (rubberband::formant_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rubberband_semitones,
  n_stereo_busses,
  param_common (
    "Amt",
    declptr<rubberband>(),
    declptr<rubberband::semitones_tag>()),
  rubberband::get_parameter (rubberband::semitones_tag {}),
  slider_ext);

parameter_cpp_class_define (
  rubberband_detune,
  n_stereo_busses,
  param_common (
    "Detune",
    declptr<rubberband>(),
    declptr<rubberband::detune_tag>()),
  rubberband::get_parameter (rubberband::detune_tag {}),
  slider_ext);

using rubberband_params = mp_list<
  rubberband_semitones,
  rubberband_detune,
  rubberband_mode,
  rubberband_formant_mode>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  soundtouch_semitones,
  n_stereo_busses,
  param_common (
    "Amt",
    declptr<soundtouch_fx>(),
    declptr<soundtouch_fx::semitones_tag>()),
  soundtouch_fx::get_parameter (soundtouch_fx::semitones_tag {}),
  slider_ext);

parameter_cpp_class_define (
  soundtouch_detune,
  n_stereo_busses,
  param_common (
    "Detune",
    declptr<soundtouch_fx>(),
    declptr<soundtouch_fx::detune_tag>()),
  soundtouch_fx::get_parameter (soundtouch_fx::detune_tag {}),
  slider_ext);

parameter_cpp_class_define (
  soundtouch_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<soundtouch_fx>(),
    declptr<soundtouch_fx::mode_tag>()),
  soundtouch_fx::get_parameter (soundtouch_fx::mode_tag {}),
  slider_ext);

using soundtouch_params
  = mp_list<soundtouch_semitones, soundtouch_detune, soundtouch_mode>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  reverb_decay,
  n_stereo_busses,
  param_common ("Late Decay", declptr<reverb>(), declptr<reverb::decay_tag>()),
  reverb::get_parameter (reverb::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_size,
  n_stereo_busses,
  param_common ("Late Size", declptr<reverb>(), declptr<reverb::size_tag>()),
  reverb::get_parameter (reverb::size_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_early_size,
  n_stereo_busses,
  param_common ("ER Size", declptr<reverb>(), declptr<reverb::er_size_tag>()),
  reverb::get_parameter (reverb::er_size_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_predelay,
  n_stereo_busses,
  param_common ("Predelay", declptr<reverb>(), declptr<reverb::predelay_tag>()),
  reverb::get_parameter (reverb::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_gap,
  n_stereo_busses,
  param_common ("Gap", declptr<reverb>(), declptr<reverb::gap_tag>()),
  reverb::get_parameter (reverb::gap_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_early_gain,
  n_stereo_busses,
  param_common (
    "ER Gain",
    declptr<reverb>(),
    declptr<reverb::early_gain_tag>()),
  reverb::get_parameter (reverb::early_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_late_gain,
  n_stereo_busses,
  param_common (
    "Late Gain",
    declptr<reverb>(),
    declptr<reverb::late_gain_tag>()),
  reverb::get_parameter (reverb::late_gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_early_2_late_bal,
  n_stereo_busses,
  param_common (
    "ER->Late",
    declptr<reverb>(),
    declptr<reverb::early_2_late_bal_tag>()),
  reverb::get_parameter (reverb::early_2_late_bal_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_mod_freq,
  n_stereo_busses,
  param_common ("Mod Freq", declptr<reverb>(), declptr<reverb::mod_freq_tag>()),
  reverb::get_parameter (reverb::mod_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_mod_depth,
  n_stereo_busses,
  param_common (
    "Mod Depth",
    declptr<reverb>(),
    declptr<reverb::mod_depth_tag>()),
  reverb::get_parameter (reverb::mod_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_mod_mode,
  n_stereo_busses,
  param_common ("Mod Mode", declptr<reverb>(), declptr<reverb::mod_mode_tag>()),
  reverb::get_parameter (reverb::mod_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_mod_spread,
  n_stereo_busses,
  param_common (
    "Mod Spread",
    declptr<reverb>(),
    declptr<reverb::mod_spread_tag>()),
  reverb::get_parameter (reverb::mod_spread_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_diff_in,
  n_stereo_busses,
  param_common (
    "Attack",
    declptr<reverb>(),
    declptr<reverb::in_diffusion_tag>()),
  reverb::get_parameter (reverb::in_diffusion_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_diff_out,
  n_stereo_busses,
  param_common (
    "Diffusion",
    declptr<reverb>(),
    declptr<reverb::out_diffusion_tag>()),
  reverb::get_parameter (reverb::out_diffusion_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_l_sparseness,
  n_stereo_busses,
  param_common (
    "L Sparse",
    declptr<reverb>(),
    declptr<reverb::l_sparseness_tag>()),
  reverb::get_parameter (reverb::l_sparseness_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_r_sparseness,
  n_stereo_busses,
  param_common (
    "R Sparse",
    declptr<reverb>(),
    declptr<reverb::r_sparseness_tag>()),
  reverb::get_parameter (reverb::r_sparseness_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_lr_sparseness,
  n_stereo_busses,
  param_common (
    "LR Sparse",
    declptr<reverb>(),
    declptr<reverb::lr_sparseness_tag>()),
  reverb::get_parameter (reverb::lr_sparseness_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_damp_factor,
  n_stereo_busses,
  param_common (
    "Damp Amt",
    declptr<reverb>(),
    declptr<reverb::damp_factor_tag>()),
  reverb::get_parameter (reverb::damp_factor_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_damp_freq,
  n_stereo_busses,
  param_common (
    "Damp FCut",
    declptr<reverb>(),
    declptr<reverb::damp_freq_tag>()),
  reverb::get_parameter (reverb::damp_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_lf_time_factor,
  n_stereo_busses,
  param_common (
    "Damp Bal",
    declptr<reverb>(),
    declptr<reverb::lf_time_factor_tag>()),
  reverb::get_parameter (reverb::lf_time_factor_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_hp_freq,
  n_stereo_busses,
  param_common ("HP FCut", declptr<reverb>(), declptr<reverb::hp_freq_tag>()),
  reverb::get_parameter (reverb::hp_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_stereo,
  n_stereo_busses,
  param_common ("Stereo", declptr<reverb>(), declptr<reverb::stereo_tag>()),
  reverb::get_parameter (reverb::stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<reverb>(),
    declptr<reverb::ducking_speed_tag>()),
  reverb::get_parameter (reverb::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  reverb_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<reverb>(),
    declptr<reverb::ducking_threshold_tag>()),
  reverb::get_parameter (reverb::ducking_threshold_tag {}),
  slider_ext);

#if 0
parameter_cpp_class_define (
  reverb_test,
  n_stereo_busses,
  param_common ("Test", declptr<reverb>(), declptr<reverb::test_param_tag>()),
  reverb::get_parameter (reverb::test_param_tag {}),
  slider_ext);
#endif

using reverb_params = mp_list<
  reverb_predelay,
  reverb_early_gain,
  reverb_early_size,
  reverb_early_2_late_bal,
  reverb_gap,
  reverb_late_gain,
  reverb_size,
  reverb_decay,

  reverb_damp_freq,
  reverb_damp_factor,
  reverb_lf_time_factor,
  reverb_hp_freq,
  reverb_mod_freq,
  reverb_mod_depth,
  reverb_mod_mode,
  reverb_mod_spread,

  reverb_diff_in,
  reverb_diff_out,
  reverb_l_sparseness,
  reverb_r_sparseness,
  reverb_lr_sparseness,
  reverb_stereo,
  reverb_ducking_threshold,
  reverb_ducking_speed>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  diffuse_delay_sixteenths,
  n_stereo_busses,
  param_common (
    "Sixteenths",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::sixteenths_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::sixteenths_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_gain,
  n_stereo_busses,
  param_common (
    "Out Gain",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::gain_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::gain_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_feedback,
  n_stereo_busses,
  param_common (
    "Feedback",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::feedback_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_damp,
  n_stereo_busses,
  param_common (
    "Damp",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::damp_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::damp_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_freq_spread,
  n_stereo_busses,
  param_common (
    "F Spread",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::freq_spread_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::freq_spread_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_tilt,
  n_stereo_busses,
  param_common (
    "Tilt",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::tilt_db_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::tilt_db_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_mode,
  n_stereo_busses,
  param_common (
    "Mode",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::mode_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_diffusion,
  n_stereo_busses,
  param_common (
    "Diffusion",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::diffusion_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::diffusion_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_ducking_speed,
  n_stereo_busses,
  param_common (
    "Duck Spd",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::ducking_speed_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_ducking_threshold,
  n_stereo_busses,
  param_common (
    "Duck Thres",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::ducking_threshold_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::ducking_threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_mod_freq,
  n_stereo_busses,
  param_common (
    "Mod Freq",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::mod_freq_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::mod_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_mod_depth,
  n_stereo_busses,
  param_common (
    "Mod Amt",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::mod_depth_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::mod_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_mod_mode,
  n_stereo_busses,
  param_common (
    "Mod Mode",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::mod_mode_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::mod_mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_desync,
  n_stereo_busses,
  param_common (
    "Desync",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::desync_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::desync_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_transients,
  n_stereo_busses,
  param_common (
    "Transient",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::transients_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::transients_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_delay_hipass,
  n_stereo_busses,
  param_common (
    "Hipass",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::hipass_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::hipass_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_bp_freq,
  n_stereo_busses,
  param_common (
    "BP Freq",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::bp_freq_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::bp_freq_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_bp_wet_dry,
  n_stereo_busses,
  param_common (
    "BP WetDry",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::bp_wet_dry_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::bp_wet_dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_bp_drive,
  n_stereo_busses,
  param_common (
    "BP Drv",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::bp_drive_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::bp_drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_bp_envfollow,
  n_stereo_busses,
  param_common (
    "BP EnvF",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::bp_envfollow_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::bp_envfollow_tag {}),
  slider_ext);

parameter_cpp_class_define (
  diffuse_bp_reso,
  n_stereo_busses,
  param_common (
    "BP Reso",
    declptr<diffuse_delay>(),
    declptr<diffuse_delay::bp_reso_tag>()),
  diffuse_delay::get_parameter (diffuse_delay::bp_reso_tag {}),
  slider_ext);

using diffuse_delay_params = mp_list<
  diffuse_delay_mode,
  diffuse_delay_sixteenths,
  diffuse_delay_diffusion,
  diffuse_delay_feedback,
  diffuse_delay_transients,
  diffuse_delay_gain,
  diffuse_delay_ducking_threshold,
  diffuse_delay_ducking_speed,

  diffuse_delay_damp,
  diffuse_delay_hipass,
  diffuse_delay_tilt,
  diffuse_delay_freq_spread,
  diffuse_delay_mod_mode,
  diffuse_delay_mod_freq,
  diffuse_delay_desync,
  diffuse_delay_mod_depth,

  diffuse_bp_wet_dry,
  diffuse_bp_freq,
  diffuse_bp_drive,
  diffuse_bp_reso,
  diffuse_bp_envfollow>;
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  mod_stages,
  n_stereo_busses,
  param_common (
    "N Stages",
    declptr<upsampled<mod>>(),
    declptr<mod::stages_tag>()),
  mod::get_parameter (mod::stages_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_lfo_rate,
  n_stereo_busses,
  param_common (
    "LFO Rate",
    declptr<upsampled<mod>>(),
    declptr<mod::lfo_rate_tag>()),
  mod::get_parameter (mod::lfo_rate_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_lfo_time_base,
  n_stereo_busses,
  param_common (
    "LFO Clock",
    declptr<upsampled<mod>>(),
    declptr<mod::lfo_time_base_tag>()),
  mod::get_parameter (mod::lfo_time_base_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_lfo_depth,
  n_stereo_busses,
  param_common (
    "LFO Amt",
    declptr<upsampled<mod>>(),
    declptr<mod::lfo_depth_tag>()),
  mod::get_parameter (mod::lfo_depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_lfo_wave,
  n_stereo_busses,
  param_common (
    "LFO Wave",
    declptr<upsampled<mod>>(),
    declptr<mod::lfo_wave_tag>()),
  mod::get_parameter (mod::lfo_wave_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_lfo_stereo,
  n_stereo_busses,
  param_common (
    "LFO Stereo",
    declptr<upsampled<mod>>(),
    declptr<mod::lfo_stereo_tag>()),
  mod::get_parameter (mod::lfo_stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_center,
  n_stereo_busses,
  param_common (
    "Center",
    declptr<upsampled<mod>>(),
    declptr<mod::center_tag>()),
  mod::get_parameter (mod::center_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_feedback,
  n_stereo_busses,
  param_common (
    "Feedback",
    declptr<upsampled<mod>>(),
    declptr<mod::feedback_tag>()),
  mod::get_parameter (mod::feedback_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_feedback_locut,
  n_stereo_busses,
  param_common (
    "Fb LoCut",
    declptr<upsampled<mod>>(),
    declptr<mod::feedback_locut_tag>()),
  mod::get_parameter (mod::feedback_locut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_feedback_hicut,
  n_stereo_busses,
  param_common (
    "Fb HiCut",
    declptr<upsampled<mod>>(),
    declptr<mod::feedback_hicut_tag>()),
  mod::get_parameter (mod::feedback_hicut_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_drive,
  n_stereo_busses,
  param_common ("Drive", declptr<upsampled<mod>>(), declptr<mod::drive_tag>()),
  mod::get_parameter (mod::drive_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_drive_curve,
  n_stereo_busses,
  param_common (
    "Drv Curve",
    declptr<upsampled<mod>>(),
    declptr<mod::drive_curve_tag>()),
  mod::get_parameter (mod::drive_curve_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_param_a,
  n_stereo_busses,
  param_common ("Param A", declptr<upsampled<mod>>(), declptr<mod::a_tag>()),
  mod::get_parameter (mod::a_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_param_b,
  n_stereo_busses,
  param_common ("Param B", declptr<upsampled<mod>>(), declptr<mod::b_tag>()),
  mod::get_parameter (mod::b_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_spread,
  n_stereo_busses,
  param_common (
    "Spread",
    declptr<upsampled<mod>>(),
    declptr<mod::spread_tag>()),
  mod::get_parameter (mod::spread_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_depth,
  n_stereo_busses,
  param_common ("Depth", declptr<upsampled<mod>>(), declptr<mod::depth_tag>()),
  mod::get_parameter (mod::depth_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_mode,
  n_stereo_busses,
  param_common ("Mode", declptr<upsampled<mod>>(), declptr<mod::mode_tag>()),
  mod::get_parameter (mod::mode_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod_oversampling,
  n_stereo_busses,
  param_common (
    "Upsample",
    declptr<upsampled<mod>>(),
    declptr<oversampled_amount_tag>()),
  upsampled<mod>::get_parameter (oversampled_amount_tag {}),
  slider_ext);

using mod_params = mp_list<
  mod_lfo_rate,
  mod_lfo_depth,
  mod_stages,
  mod_depth,
  mod_center,
  mod_spread,
  mod_param_a,
  mod_param_b,
  mod_mode,
  mod_lfo_stereo,
  mod_lfo_wave,
  mod_feedback,
  mod_drive,
  mod_drive_curve,
  mod_feedback_locut,
  mod_feedback_hicut,
  mod_lfo_time_base,
  mod_oversampling>;
//------------------------------------------------------------------------------

#if 0
parameter_cpp_class_define (
  polyphase_fir_test_gain,
  n_stereo_busses,
  param_common (
    "Gain",
    declptr<polyphase_fir_test>(),
    declptr<polyphase_fir_test::gain_tag>()),
  polyphase_fir_test::get_parameter (polyphase_fir_test::gain_tag {}),
  slider_ext);

using polyphase_fir_test_params = mp_list<polyphase_fir_test_gain>;
#endif
//------------------------------------------------------------------------------
#if 0 // experiments
}
} // namespace artv::parameters

#include "artv-common/dsp/own/fx/experiments.hpp"

namespace artv {
namespace parameters {

parameter_cpp_class_define (
  experiments_a,
  n_stereo_busses,
  param_common (
    "A",
    declptr<experiments>(),
    declptr<experiments::param_a_tag>()),
  experiments::get_parameter (experiments::param_a_tag {}),
  slider_ext);

parameter_cpp_class_define (
  experiments_b,
  n_stereo_busses,
  param_common (
    "B",
    declptr<experiments>(),
    declptr<experiments::param_b_tag>()),
  experiments::get_parameter (experiments::param_b_tag {}),
  slider_ext);

parameter_cpp_class_define (
  experiments_c,
  n_stereo_busses,
  param_common (
    "C",
    declptr<experiments>(),
    declptr<experiments::param_c_tag>()),
  experiments::get_parameter (experiments::param_c_tag {}),
  slider_ext);

parameter_cpp_class_define (
  experiments_d,
  n_stereo_busses,
  param_common (
    "D",
    declptr<experiments>(),
    declptr<experiments::param_d_tag>()),
  experiments::get_parameter (experiments::param_d_tag {}),
  slider_ext);

using experiments_params
  = mp_list<experiments_a, experiments_b, experiments_c, experiments_d>;

#endif // experiments

#define TWEAK_BUILD 1

#if TWEAK_BUILD
using all_fx_typelists = mp_list<
  lr_crossv_params,
  wonky_crossv_params,
  lin_iir_crossv_params,
  mod_params>;

static constexpr auto fx_choices
  = make_cstr_array ("none", "LR", "Wonky", "lin IIR", "1. Mod");
#else
// clang-format off
using all_fx_typelists = mp_list<
  consoles_params,
  bbd_echo_params,
  sound_delay_params,
  echo_cycles_params,
  spring_box_params,
  busscolors4_params,
  nonlinear_params,
  sandwitch_amp_params,
  _1175_params,
  consolidator_params,
  gate_expander_params,
  fairly_childish_params,
  major_tom_params,
  slax_params,
  transience_params,
  event_horizon_2_params,
  _4x4_params,
  bass_professor_params,
  hugebooty_params,
  bbe_params,
  rbj1073_params,
  luftikus_params,
  eq4x_params,
  chow_phaser_params,
  myphaser_params,
  ripple_params,
  zebigchorus3_params,
  zelittlechorus_params,
  atlantis_reverb_params,
  df_er_params,
  df_hall_params,
  df_plate_params,
  df_room_params,
  tal_reverb2_params,
  stereo_bub3_params,
  stereo_tilt_params,
  filter2x_params,
  fdnverb_params,
  track_comp_params,
  signal_crusher_params,
  waveshaper_params,
  lr_crossv_params,
  wonky_crossv_params,
  lin_iir_crossv_params,
  transient_gate_params,
  lin_eq4x_params,
  naive_pitch_params,
  rubberband_params,
  soundtouch_params,
  reverb_params,
  diffuse_delay_params,
  mod_params>;
// clang-format on
//------------------------------------------------------------------------------
// ordering between "fx_choices" and "all_fx_params_typelists" MUST match!
// ordering must not be changed between versions, otherwise it will break old
// versions automation or presets.

static constexpr auto fx_choices = make_cstr_array (
  "-none-",
  ":Compand/Expand Aw Consoles",
  ":Delay BBD",
  ":Delay Sample Delay",
  ":Delay Echo Cycles",
  ":Delay Spring Box",
  ":Distortion Busscolors4",
  ":Distortion Nonlinear",
  ":Distortion Sandwitch Amp",
  ":Dynamics 1175",
  ":Dynamics Consolidator",
  ":Dynamics Gate/Expander",
  ":Dynamics FairlyChildish",
  ":Dynamics Major Tom",
  ":Dynamics Slax (lala)",
  ":Dynamics Transience",
  ":Dynamics Evt Horiz (clip)",
  ":Exciter 4x4",
  ":Exciter Bass Professor",
  ":Exciter Huge Booty",
  ":Exciter Sonic Enhancer",
  ":EQ RBJ 1073",
  ":EQ Luftikus",
  ":EQ 4-band",
  ":Modulation: Chow Phaser",
  ":Modulation: ArtV OldPhaser",
  ":Modulation: Ripple Phaser",
  ":Modulation: Ze Big Chorus",
  ":Modulation: Ze Little Chorus",
  ":Reverb Atlantis",
  ":Reverb Dragonfly Early",
  ":Reverb Dragonfly Hall",
  ":Reverb Dragonfly Plate",
  ":Reverb Dragonfly Room",
  ":Reverb Tal 2",
  ":Stereo Stereo Bub 3",
  ":Stereo Stereo Tilt",
  ":Filter Filters x2",
  ":Reverb FDN Verb Riser",
  ":Dynamics Track Comp",
  ":Distortion Signal Crusher",
  ":Distortion ArtV Shaper",
  ":Filter Crossover LR",
  ":Filter Crossover WTF",
  ":Filter Crossover Lin",
  ":Dynamics Transient Gate",
  ":EQ 4-band LinPh",
  ":Pitch Naive",
  ":Pitch Rubberband",
  ":Pitch Soundtouch",
  ":Reverb Artv Reverb",
  ":Delay Artv DelayVerb",
  ":Modulation: ArtV Mod");

#endif // #if TWEAK_BUILD

static_assert (
  fx_choices.size() - 1 == mp11::mp_size<all_fx_typelists>::value,
  "");

parameter_cpp_class_define (
  fx_type,
  n_stereo_busses,
  param_common ("Effect"),
  choice_param (0, fx_choices, 75, true),
  combobox_ext);
//------------------------------------------------------------------------------
using channel_fx_sliders_typelist = mp11::mp_flatten<all_fx_typelists>;

using main_page_sliders_typelist
  = mp11::mp_list<dry_pan, dry_balance, wet_pan, wet_balance, fx_mix, pan>;

using channel_sliders_typelist
  = mp11::mp_push_front<main_page_sliders_typelist, volume>;

using all_nonfx_sliders_typelist
  = mp11::mp_push_back<channel_sliders_typelist, global_volume>;

using routing_controls_typelist
  = mp11::mp_list<routing, in_selection, out_selection, mixer_sends>;

using non_routing_controls_typelist
  = mp11::mp_list<mute_solo, fx_type, channel_modifs>;

using all_controls_typelist = mp11::mp_flatten<
  mp11::mp_list<routing_controls_typelist, non_routing_controls_typelist>>;

using all_channel_sliders_typelist = mp11::mp_flatten<
  mp11::mp_list<all_nonfx_sliders_typelist, channel_fx_sliders_typelist>>;

using parameters_typelist = mp11::mp_flatten<
  mp11::mp_list<all_controls_typelist, all_channel_sliders_typelist>>;

//------------------------------------------------------------------------------
static constexpr auto channel_text_keys = make_cstr_array (
  "bus1_text",
  "bus2_text",
  "bus3_text",
  "bus4_text",
  "bus5_text",
  "bus6_text",
  "bus7_text",
  "bus8_text");

static_assert (channel_text_keys.size() == n_stereo_busses, "");
//------------------------------------------------------------------------------
// clang-format off
static constexpr char const* about_text =
"v" VERSION_TXT "\n"
"\n"
"CONTROLS\n"
"--------\n"
"\n"
"\"LMB\"= Left Mouse Button. \"RMB\"= Right Mouse Button.\n"
"\n"
"-Double click a knob/slider: Set to default value.\n"
"-RMB on a knob/slider: Keyboard data entry.\n"
"-Ctrl + Dragging a knob/slider: Precise knob/slider adjustment.\n"
"-Mouse Wheel: Same as dragging a knob/slider.\n"
"-RMB drag and drop on FX comboboxes: FX parameters Swap.\n"
"-Ctrl + RMB drag and drop on FX comboboxes: FX parameters Copy.\n"
"-Ctrl + Alt + RMB drag and drop on FX comboboxes: FX + Pan/Bal/MS parameters Copy.\n"
"\n"
"CREDITS\n"
"--------"
"\n"
"Powered by Open Source DSP and libraries. Credits (In no particular order):\n"
"\n"
"From Airwindows (Chris Johnson).\n"
"https://www.patreon.com/airwindows \n"
"\n"
"-BussColors4.\n"
"-Aw Consoles: Console4, 5, 6 and 7.\n"
"-Density2.\n"
"\n"
"From Chokehold.\n"
"http://free.chokehold.net/\n"
"\n"
"-Consolidator.\n"
"-Gate/Expander.\n"
"-Track Comp.\n"
"-Signal Crusher.\n"
"\n"
"From Chow (Jatin Chowdhury).\n"
"https://jatinchowdhury18.medium.com/\n"
"\n"
"-Chow Phaser.\n"
"\n"
"From Geraint Luff:\n"
"Can be reached for commercial work at: https://signalsmith-audio.co.uk/\n"
"https://github.com/geraintluff/jsfx\n"
"\n"
"-Atlantis: Atlantis Reverb.\n"
"-Echo Cycles.\n"
"-Ripple Phaser.\n"
"-Spring Delay.\n"
"-Sandwitch Amp.\n"
"\n"
"From Liteon (Lubomir I. Ivanov):\n"
"\n"
"-Sonic Enhancer.\n"
"-Presence: Presence EQ (Moorer).\n"
"-Nonlinear: Non-linear Processor.\n"
"-Tilt: Tilt (Mono-To-Stereo).\n"
"\n"
"From LKJB :\n"
"https://lkjbdsp.wordpress.com/luftikus/\n"
"https://github.com/lkjbdsp/lkjb-plugins\n"
"\n"
"-Luftikus.\n"
"\n"
"From Michael Willis:\n"
"https://michaelwillis.github.io/dragonfly-reverb/\n"
"\n"
"-Dragonfly Early: Dragonfly Early Reflections.\n"
"-Dragonfly Hall.\n"
"-Dragonfly Plate.\n"
"-Dragonfly Room.\n"
"\n"
"From Saike (Joep Vanlier):\n"
"https://github.com/JoepVanlier\n"
"\n"
"-Korg, Ladder and Steiner filters from Yutani.\n"
"-Stereo Bub 3.\n"
"-Transience.\n"
"\n"
"From Shabronic (S.D.Smith):\n"
"https://github.com/shabtronic\n"
"\n"
"-FDN Reverb Riser.\n"
"\n"
"From Smashed Transistors (Thierry Rochebois):\n"
"https://www.youtube.com/channel/UCAhRo1cCl1r_dFMVYaxTh5Q/videos\n"
"\n"
"-Ze Big Chorus 03.\n"
"-Ze Little Scanner Chorus.\n"
"\n"
"From Sonic Anomaly (Stige T.):\n"
"https://github.com/Sonic-Anomaly/Sonic-Anomaly-JSFX\n"
"\n"
"-Bass Professor: Bass Professor (mkI).\n"
"-Slax (lala): Slax.\n"
"\n"
"From Stillwell (Thomas Scott Stillwell):\n"
"Whose FX included here are primitive JSFX drafts of his excellent commercial\n"
"work (+Fairly priced, +unintrusive copy protection).\n"
"https://www.stillwellaudio.com/\n"
"\n"
"-1175.\n"
"-4x4.\n"
"-Event Horizon 2.\n"
"-FairlyChildish.\n"
"-Major Tim: Major Tom.\n"
"-Huge Booty: Huge Booty Bass Enhancer.\n"
"-RBJ 1073.\n"
"\n"
"From TAL (Patrick Kunz):\n"
"His commercial work found here. He needs no introduction:\n"
"https://tal-software.com/\n"
"\n"
"-TAL Reverb 2.\n"
"\n"
"From Witti:\n"
"https://stash.reaper.fm/v/25168/js_plugins.zip\n"
"\n"
"-BBD Mono: bbd_echo.\n"
"\n"
"Built using the JUCE framework. Thanks for the excellent free support!\n"
"Awesome forums too:\n"
"https://juce.com/\n"
"https://forum.juce.com/\n"
"\n"
"Actually using the LV2 JUCE fork:\n"
"https://github.com/lv2-porting-project/JUCE\n"
"\n"
"Using some pieces of Cockos WDL:\n"
"https://www.cockos.com/wdl/\n"
"\n"
"Using muFFT (for float):\n"
"https://github.com/Themaister/muFFT\n"
"\n"
"Using PFFFT port (for double):\n"
"https://github.com/marton78/pffft\n"
"\n"
"Using XSIMD:\n"
"https://github.com/xtensor-stack/xsimd.git\n"
"\n"
"Using GCEM:\n"
"https://github.com/kthohr/gcem.git\n"
"\n"
"Using FF-meters:\n"
"https://github.com/ffAudio/ff_meters\n"
"\n"
"Using the Rubberband library:\n"
"https://breakfastquay.com/rubberband/\n"
"\n"
"Using the SoundTouch library:\n"
"https://www.surina.net/soundtouch/\n"
"\n"
"Some code is ported from DISHTRO-ports instead of its original repos:\n"
"https://github.com/DISTRHO/DISTRHO-Ports\n"
"\n"
"For the full licenses of each used project see the 3RD-PARTY-LICENSES file:\n"
"https://github.com/RafaGago/artv-audio/blob/master/3RD-PARTY-LICENSES.\n"
"\n"
"The sources of projects have been used as references for doing own\n"
"implementations:\n"
"\n"
"RS-MET library by Robin Schmidt:\n"
"https://github.com/RobinSchmidt/RS-MET\n"
"\n"
"csound:\n"
"https://github.com/csound/csound\n"
"\n"
"ReEq by Justin Johnson:\n"
"https://github.com/Justin-Johnson/ReJJ\n"
"\n"
"ADAA by Jatin Chowdhury (CHOW):\n"
"https://github.com/jatinchowdhury18/ADAA\n"
"\n"
"Liquid DSP:\n"
"https://github.com/jgaeddert/liquid-dsp\n"
"\n"
"Special thanks to the folks on the KVR DSP forum. From the top of my head,\n"
"(not everyone might be mentioned, but I aim to):\n"
"\n"
"-mystran (Teemu Voipio)\n"
"-martinvicanek (Martin Vicanek). The linear phase crossover is impleming\n"
" his \"A New Reverse IIR Filtering Algorithm\" paper."
"-rs-met (Robin Schmidt)\n"
"-Z1202 (Vadim Zavalishin). Author of the also helpful \"The art of va filter\n"
" design\" free book.\n"
"\n"
"A copy of all the licenses for all the third party used by this software can\n"
"be found on:\n"
"https://github.com/RafaGago/artv-audio/blob/master/3RD-PARTY-LICENSES\n"
;
// clang-format on
}} // namespace artv::parameters
