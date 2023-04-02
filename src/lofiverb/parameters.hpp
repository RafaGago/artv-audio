#pragma once

#include "artv-common/dsp/own/classes/mix.hpp"
#include "artv-common/dsp/own/fx/lofiverb.hpp"
#include "artv-common/juce/parameters.hpp"

namespace artv { namespace parameters {
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  algorithm,
  1,
  param_common (
    "Algo",
    declptr<lofiverb>(),
    declptr<lofiverb::algorithm_tag>()),
  lofiverb::get_parameter (lofiverb::algorithm_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  mode,
  1,
  param_common ("Mode", declptr<lofiverb>(), declptr<lofiverb::mode_tag>()),
  lofiverb::get_parameter (lofiverb::mode_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  clock,
  1,
  param_common ("Clock", declptr<lofiverb>(), declptr<lofiverb::clock_tag>()),
  lofiverb::get_parameter (lofiverb::clock_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  character,
  1,
  param_common ("???", declptr<lofiverb>(), declptr<lofiverb::character_tag>()),
  lofiverb::get_parameter (lofiverb::character_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lf_amt,
  1,
  param_common ("Lows", declptr<lofiverb>(), declptr<lofiverb::lf_amt_tag>()),
  lofiverb::get_parameter (lofiverb::lf_amt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  hf_amt,
  1,
  param_common ("Highs", declptr<lofiverb>(), declptr<lofiverb::hf_amt_tag>()),
  lofiverb::get_parameter (lofiverb::hf_amt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  decay,
  1,
  param_common ("Decay", declptr<lofiverb>(), declptr<lofiverb::decay_tag>()),
  lofiverb::get_parameter (lofiverb::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  predelay,
  1,
  param_common (
    "Predelay",
    declptr<lofiverb>(),
    declptr<lofiverb::predelay_tag>()),
  lofiverb::get_parameter (lofiverb::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  operating_range,
  1,
  param_common (
    "Op Range",
    declptr<lofiverb>(),
    declptr<lofiverb::clip_level_tag>()),
  lofiverb::get_parameter (lofiverb::clip_level_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod,
  1,
  param_common ("Mod", declptr<lofiverb>(), declptr<lofiverb::mod_tag>()),
  lofiverb::get_parameter (lofiverb::mod_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo,
  1,
  param_common ("Stereo", declptr<lofiverb>(), declptr<lofiverb::stereo_tag>()),
  lofiverb::get_parameter (lofiverb::stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ducking_speed,
  1,
  param_common (
    "Duck Time",
    declptr<lofiverb>(),
    declptr<lofiverb::ducking_speed_tag>()),
  lofiverb::get_parameter (lofiverb::ducking_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  ducking_threshold,
  1,
  param_common (
    "Duck Thrs",
    declptr<lofiverb>(),
    declptr<lofiverb::ducking_threshold_tag>()),
  lofiverb::get_parameter (lofiverb::ducking_threshold_tag {}),
  slider_ext);

parameter_cpp_class_define (
  dry,
  1,
  param_common (
    "Dry",
    declptr<dry_wet_mixer>(),
    declptr<dry_wet_mixer::dry_tag>()),
  dry_wet_mixer::get_parameter (dry_wet_mixer::dry_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wet,
  1,
  param_common (
    "Wet",
    declptr<dry_wet_mixer>(),
    declptr<dry_wet_mixer::wet_tag>()),
  dry_wet_mixer::get_parameter (dry_wet_mixer::wet_tag {}),
  slider_ext);

parameter_cpp_class_define (
  wet_pan,
  1,
  param_common (
    "Wet Pan",
    declptr<dry_wet_mixer>(),
    declptr<dry_wet_mixer::wet_pan_tag>()),
  dry_wet_mixer::get_parameter (dry_wet_mixer::wet_pan_tag {}),
  slider_ext);

using lofiverb_parameters = mp_list<
  algorithm,
  predelay,
  clock,
  decay,
  lf_amt,
  hf_amt,
  stereo,
  character,
  mod,
  ducking_threshold,
  ducking_speed,
  mode,
  operating_range>;

using mixer_parameters = mp_list<dry, wet, wet_pan>;

using parameters_typelist
  = mp11::mp_flatten<mp11::mp_list<lofiverb_parameters, mixer_parameters>>;

}} // namespace artv::parameters
