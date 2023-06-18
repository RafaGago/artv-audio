#pragma once

#include "artv-common/dsp/own/classes/mix.hpp"
#include "artv-common/dsp/own/fx/turbopaco.hpp"
#include "artv-common/juce/parameters.hpp"

namespace artv { namespace parameters {
//------------------------------------------------------------------------------
parameter_cpp_class_define (
  algorithm,
  1,
  param_common (
    "Algo",
    declptr<turbopaco>(),
    declptr<turbopaco::algorithm_tag>()),
  turbopaco::get_parameter (turbopaco::algorithm_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  mode,
  1,
  param_common ("Mode", declptr<turbopaco>(), declptr<turbopaco::mode_tag>()),
  turbopaco::get_parameter (turbopaco::mode_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  clock,
  1,
  param_common ("Clock", declptr<turbopaco>(), declptr<turbopaco::clock_tag>()),
  turbopaco::get_parameter (turbopaco::clock_tag {}),
  combobox_ext);

parameter_cpp_class_define (
  character,
  1,
  param_common (
    "???",
    declptr<turbopaco>(),
    declptr<turbopaco::character_tag>()),
  turbopaco::get_parameter (turbopaco::character_tag {}),
  slider_ext);

parameter_cpp_class_define (
  lf_amt,
  1,
  param_common ("Lows", declptr<turbopaco>(), declptr<turbopaco::lf_amt_tag>()),
  turbopaco::get_parameter (turbopaco::lf_amt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  hf_amt,
  1,
  param_common (
    "Highs",
    declptr<turbopaco>(),
    declptr<turbopaco::hf_amt_tag>()),
  turbopaco::get_parameter (turbopaco::hf_amt_tag {}),
  slider_ext);

parameter_cpp_class_define (
  decay,
  1,
  param_common ("Decay", declptr<turbopaco>(), declptr<turbopaco::decay_tag>()),
  turbopaco::get_parameter (turbopaco::decay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  predelay,
  1,
  param_common (
    "Predelay",
    declptr<turbopaco>(),
    declptr<turbopaco::predelay_tag>()),
  turbopaco::get_parameter (turbopaco::predelay_tag {}),
  slider_ext);

parameter_cpp_class_define (
  operating_range,
  1,
  param_common (
    "Op Range",
    declptr<turbopaco>(),
    declptr<turbopaco::clip_level_tag>()),
  turbopaco::get_parameter (turbopaco::clip_level_tag {}),
  slider_ext);

parameter_cpp_class_define (
  mod,
  1,
  param_common ("Mod", declptr<turbopaco>(), declptr<turbopaco::mod_tag>()),
  turbopaco::get_parameter (turbopaco::mod_tag {}),
  slider_ext);

parameter_cpp_class_define (
  stereo,
  1,
  param_common (
    "Stereo",
    declptr<turbopaco>(),
    declptr<turbopaco::stereo_tag>()),
  turbopaco::get_parameter (turbopaco::stereo_tag {}),
  slider_ext);

parameter_cpp_class_define (
  dyn_speed,
  1,
  param_common (
    "Dyn Time",
    declptr<turbopaco>(),
    declptr<turbopaco::dyn_speed_tag>()),
  turbopaco::get_parameter (turbopaco::dyn_speed_tag {}),
  slider_ext);

parameter_cpp_class_define (
  dyn_threshold,
  1,
  param_common (
    "Dyn Thres",
    declptr<turbopaco>(),
    declptr<turbopaco::dyn_threshold_tag>()),
  turbopaco::get_parameter (turbopaco::dyn_threshold_tag {}),
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

static constexpr auto wet_param = lambda_forward (
  dry_wet_mixer::get_parameter (dry_wet_mixer::wet_tag {}),
  [] (auto v) {
    v.max      = 12.f;
    v.defaultv = -6.f;
    return v;
  });

parameter_cpp_class_define (
  wet,
  1,
  param_common (
    "Wet",
    declptr<dry_wet_mixer>(),
    declptr<dry_wet_mixer::wet_tag>()),
  wet_param,
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

using turbopaco_parameters = mp_list<
  algorithm,
  predelay,
  clock,
  decay,
  lf_amt,
  hf_amt,
  stereo,
  character,
  mod,
  dyn_threshold,
  dyn_speed,
  mode,
  operating_range>;

using mixer_parameters = mp_list<dry, wet, wet_pan>;

using parameters_typelist
  = mp11::mp_flatten<mp11::mp_list<turbopaco_parameters, mixer_parameters>>;

}} // namespace artv::parameters
