#pragma once

#include <array>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/juce/value_string_wrappers.hpp"
#include "artv-common/misc/bits.hpp"

namespace artv {

// All these "parameter_xxx" classes have free functions defining them because
// brace initialization doesn't play well with macros, so instead of writing
// constructors we do from functions, which support type deduction on C++14.
// -----------------------------------------------------------------------------
struct float_parameter {
  using value_type = float;

  char const* units;
  float       min, max, defaultv, interval, skew;
  bool        simmetric_skew;
};

template <class ValueStringWrapper = value_string::same_text_as_value>
struct float_parameter_ext : public float_parameter {
  using value_type           = float_parameter::value_type;
  using value_string_wrapper = ValueStringWrapper;
};

template <class ValueStringWrapper = value_string::same_text_as_value>
static constexpr auto float_param (
  char const*        units,
  float              min,
  float              max,
  float              defaultv,
  float              interval       = 0.f,
  float              skew           = 1.0f,
  bool               simmetric_skew = false,
  ValueStringWrapper value_string   = value_string::same_text_as_value {})
{
  return float_parameter_ext<ValueStringWrapper> {
    {units, min, max, defaultv, interval, skew, simmetric_skew}};
}
// -----------------------------------------------------------------------------
struct int_parameter {
  using value_type = int;

  char const* units;
  int         min, max, defaultv;
};

template <class ValueStringWrapper = value_string::same_text_as_value>
struct int_parameter_ext : public int_parameter {
  using value_type           = int_parameter::value_type;
  using value_string_wrapper = ValueStringWrapper;
};

template <class ValueStringWrapper = value_string::same_text_as_value>
static constexpr auto int_param (
  char const*        units,
  int                min,
  int                max,
  int                defaultv,
  ValueStringWrapper value_string = value_string::same_text_as_value {})
{
  return int_parameter_ext<ValueStringWrapper> {units, min, max, defaultv};
}
// -----------------------------------------------------------------------------
// workaround specialized case for use with a combobox, as those require the
// options passed to the widget
template <size_t N>
struct choice_parameter : public int_parameter_ext<> {
  using choices_type = std::array<char const*, N>;

  constexpr choice_parameter (
    int                 defaultv,
    choices_type const& choices_arr,
    uint                n_future_choices,
    bool                sorted)
    : int_parameter_ext {{"", 0, std::max<int> (N, n_future_choices), defaultv}}
    , choices {choices_arr}
    , alphabetically_sorted {sorted}

  {}
  choices_type choices;
  bool         alphabetically_sorted; // only available on for comboboxes
};

namespace detail {
template <class T>
struct is_choice : public std::false_type {};

template <size_t N>
struct is_choice<choice_parameter<N>> : public std::true_type {};

} // namespace detail

template <class T>
static constexpr bool is_choice = detail::is_choice<T>::value;

template <size_t N>
static constexpr choice_parameter<N> choice_param (
  int                               defaultv,
  std::array<char const*, N> const& choices,
  uint                              total_future_choices = 0,
  bool                              sorted               = false) // if 0 == N
{
  return {defaultv, choices, total_future_choices, sorted};
}
// -----------------------------------------------------------------------------
// another specialized case for use with an array of buttons representing just
// a bit on one parameter.
// -----------------------------------------------------------------------------
template <size_t N>
struct toggle_buttons_parameter : public int_parameter_ext<> {
  using texts_type = std::array<char const*, N>;

  // fixed at 23-bit, so they can be expanded without breaking past
  // automation/presets
  constexpr toggle_buttons_parameter (
    int               defaultv,
    texts_type const& texts_arr,
    uint              range_bits)
    // clang-format off
    : int_parameter_ext {{
        "",
        0,
        (int) (range_bits ? lsb_mask<uint> (range_bits) : lsb_mask<uint> (N)),
        defaultv
        }}
    , texts {texts_arr}
  // clang-format on
  {}
  // the parameter has to be mapped in a float parameter, so it can only contain
  // bits on the fraction part of a float.
  static_assert (N <= 23, "Mapping space (float32's fraction bits) exceeded");

  texts_type texts;
};

template <size_t N>
static constexpr toggle_buttons_parameter<N> toggle_buttons_param (
  int                               defaultv,
  std::array<char const*, N> const& texts,
  uint                              parameter_range_bits = 0)
{
  // assert (range_bits <= 23); // can't be constexpr...
  return {defaultv, texts, parameter_range_bits};
}
// -----------------------------------------------------------------------------
} // namespace artv
