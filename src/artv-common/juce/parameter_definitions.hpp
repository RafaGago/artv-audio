#pragma once

#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/midi.hpp"

namespace artv {

//------------------------------------------------------------------------------
// frequency parameters represented on the Host as MIDI notes but displaying
// both Hz and notes.
static constexpr auto frequency_parameter (
  float min_hz,
  float max_hz,
  float default_hz)
{
  return float_param (
    "",
    constexpr_hz_to_midi_note (min_hz),
    constexpr_hz_to_midi_note (max_hz),
    constexpr_hz_to_midi_note (default_hz),
    0.f,
    1.f,
    false,
    value_string::frequency {});
}
//------------------------------------------------------------------------------
static constexpr auto frequency_parameter_from_zero (
  float max_hz,
  float default_hz = midi_note_to_hz_min)
{
  return frequency_parameter (midi_note_to_hz_min, max_hz, default_hz);
}

}; // namespace artv
