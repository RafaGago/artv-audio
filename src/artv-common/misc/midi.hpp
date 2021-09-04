#pragma once

#include <array>

#include <gcem.hpp>

namespace artv {

// 0Hz is equal to MIDI note -infinity. We cut at "midi_note_to_hz_min" so we
// can reach zero. Note that this affects LFOs.
static constexpr int hz_to_midi_note_min = -80;
// 9.823666718985751e-06 = midi_note_to_hz (hz_to_midi_note_min);
static constexpr double midi_note_to_hz_min = 0.08047547776193117;
//------------------------------------------------------------------------------
static constexpr double constexpr_midi_note_to_hz (double note)
{
  double hz      = gcem::pow (2., (note - 69.) * (1. / 12.)) * 440.;
  bool   is_zero = hz < (midi_note_to_hz_min + (midi_note_to_hz_min * 0.1));
  return is_zero ? 0. : hz;
}
//------------------------------------------------------------------------------
static constexpr double constexpr_hz_to_midi_note (double hz)
{
  return (12. * gcem::log2 (hz * (1. / 440.))) + 69.;
}
//------------------------------------------------------------------------------
static double midi_note_to_hz (double note)
{
  double hz      = exp2 ((note - 69.) * (1. / 12.)) * 440.;
  bool   is_zero = hz < (midi_note_to_hz_min + (midi_note_to_hz_min * 0.1));
  return is_zero ? 0. : hz;
}
//------------------------------------------------------------------------------
static double hz_to_midi_note (double hz)
{
  return (12. * log2 (hz * (1. / 440.))) + 69.;
}
//------------------------------------------------------------------------------
constexpr std::array<char, 4> midi_note_to_str (int midi_note)
{
  constexpr char      note_arr[] = "c-c#d-d#e-f-f#g-g#a-a#b-";
  std::array<char, 4> ret {}; // 0 init

  if (midi_note <= (hz_to_midi_note_min + 1)) {
    // incorrect shortcut. + 1 shouldn't be added
    return ret;
  }

  int octave = midi_note / 12;
  int note   = midi_note % 12;

  note   = note >= 0 ? note : -note;
  ret[0] = note_arr[note * 2];
  ret[1] = note_arr[(note * 2) + 1];

  octave -= 1; // C-0 is note 12.
  ret[2] = (octave >= 0 && octave < 10) ? '0' + octave : 'x';
  return ret;
}
//------------------------------------------------------------------------------
} // namespace artv
