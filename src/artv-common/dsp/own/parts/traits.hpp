#pragma once

namespace artv {

// Tick function related tag
struct single_coeff_set_tag {};

// Filter types
struct lowpass_tag {};
struct highpass_tag {};
struct bandpass_tag {};
struct notch_tag {};
struct peak_tag {};
struct allpass_tag {};
struct bell_tag {}; // as in a parametric EQ
struct lowshelf_tag {}; // as in a parametric EQ
struct highshelf_tag {}; // as in a parametric EQ

// Random tags used elsewhere.
template <uint n>
struct part_tick_tag {};

}; // namespace artv
