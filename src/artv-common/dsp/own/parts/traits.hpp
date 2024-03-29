#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

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
using thiran_tag          = tag<filter_type_tag, 11>;
using raw_tag             = tag<filter_type_tag, 12>;
using lowshelf_naive_tag  = tag<filter_type_tag, 13>;
using highshelf_naive_tag = tag<filter_type_tag, 14>;

// Random tags used elsewhere.
template <uint n>
struct part_reset_coeffs_tag {};

template <uint n>
struct part_tick_tag {};

template <uint n>
struct quality_tag {};

struct no_prewarp {};
template <uint n>
struct mode_tag {};

}; // namespace artv
