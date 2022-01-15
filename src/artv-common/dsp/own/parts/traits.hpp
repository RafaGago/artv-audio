#pragma once

#include "artv-common/misc/util.hpp"

namespace artv {

// can double as compile time and runtime tag by accessing ".value"
template <class type, uint N>
struct tag : public k_uint<N> {};

// Tick function related tag
struct single_coeff_set_tag {};

struct filter_type_tag {};

// Filter types
using lowpass_tag         = tag<filter_type_tag, 0>;
using highpass_tag        = tag<filter_type_tag, 1>;
using bandpass_tag        = tag<filter_type_tag, 2>; // 0dB gain peak
using bandpass_q_gain_tag = tag<filter_type_tag, 3>; // gain depends on Q
using notch_tag           = tag<filter_type_tag, 4>;
using peak_tag            = tag<filter_type_tag, 5>;
using allpass_tag         = tag<filter_type_tag, 6>;
using bell_tag            = tag<filter_type_tag, 7>; // as in a parametric EQ
using bell_bandpass_tag   = tag<filter_type_tag, 8>; // x + in = bell.
using lowshelf_tag        = tag<filter_type_tag, 9>;
using highshelf_tag       = tag<filter_type_tag, 10>;

// Random tags used elsewhere.
template <uint n>
struct part_reset_coeffs_tag {};

template <uint n>
struct part_tick_tag {};

}; // namespace artv
