#pragma once

#include <cassert>

#include "artv-common/misc/midi.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

// lambdas can't be default constructed, even without captures, until C++20, so
// I need to wrap them in classes.

namespace value_string {

struct same_text_as_value {
  // AudioParameterFloat/Int take a std::function. nullptr sets them to nothing.
  struct string_from_value {
    static auto get() { return nullptr; }
  };

  struct value_from_string {
    static auto get() { return nullptr; }
  };
};

struct dB {
  struct string_from_value {
    static auto get()
    {
      return [] (float v, int maxlen) { return juce::String {gain_to_db (v)}; };
    }
  };

  struct value_from_string {
    static auto get()
    {
      return
        [] (juce::String const& s) { return db_to_gain (s.getFloatValue()); };
    }
  };
};

// frequency on MIDI note.
struct frequency {
  struct string_from_value {
    static auto get()
    {
      return [] (float v, int maxlen) {
        auto hz = midi_note_to_hz (v);

        juce::String s;
        if (hz < 1.f) {
          s = juce::String {hz, 5};
        }
        else if (hz < 10.f) {
          s = juce::String {hz, 4};
        }
        else if (hz < 100.f) {
          s = juce::String {hz, 3};
        }
        else if (hz < 1000.f) {
          s = juce::String {hz, 2};
        }
        else if (hz < 10000.f) {
          s = juce::String {hz * 0.001, 4};
          s.append ("k", 1);
        }
        else {
          s = juce::String {hz * 0.001, 3};
          s.append ("k", 1);
        }
        s.append ("Hz (", 4);
        s.append (midi_note_to_str ((int) (v + 0.5)).data(), 3);
        s.append (")", 1);
        return s;
      };
    }
  };

  struct value_from_string {
    static auto get()
    {
      return [] (juce::String const& str) {
        auto s   = str.toLowerCase();
        auto end = s.indexOfChar ('h');
        if (end > 0) {
          s = s.substring (0, end);
        }
        float fac = 1.f;
        if (s.endsWith ("k")) {
          fac = 1000.f;
          s   = s.substring (0, s.length() - 1);
        }
        return hz_to_midi_note (s.getFloatValue() * fac);
      };
    }
  };
};

} // namespace value_string

} // namespace artv
