#pragma once

#include <array>
#include <cassert>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <type_traits>

#include "artv-common/juce/widgets.hpp"
#include "artv-common/misc/range.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// -----------------------------------------------------------------------------
template <typename... args>
static void set_color (int id, juce::Colour color, args&&... components)
{
  apply (
    [=] (juce::Component& comp) { comp.setColour (id, color); },
    std::forward<args> (components)...);
}
// -----------------------------------------------------------------------------
namespace detail {

template <size_t N>
struct grid_row_idx {
public:
  grid_row_idx<N> operator++()
  {
    ++v;
    v %= N;
    return *this;
  }

  size_t operator()() const { return v; }

private:
  size_t v = 0;
};

// end of recursion
template <class T, size_t N>
static juce::Rectangle<T> grid (
  juce::Rectangle<T> area_remainder,
  juce::Rectangle<T>,
  float,
  std::array<float, N> const&,
  grid_row_idx<N>)
{
  return area_remainder;
}

// variant consuming a grid_pack
template <class T, class C, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  contiguous_range<C>         pack,
  args&&... vargs)
{
  for (int i = 0; i < pack.size(); ++i) {
    if (row() == 0) {
      column = area_remainder.removeFromLeft ((int) width);
    }
    pack[i].setBounds (column.removeFromTop ((int) heights[row()]));
    ++row;
  }
  return grid (
    area_remainder, column, width, heights, row, std::forward<args> (vargs)...);
}

// variant consuming an std::array
template <class T, class C, size_t Cn, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  std::array<C, Cn>&          arr,
  args&&... vargs)
{
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    contiguous_range (arr),
    std::forward<args> (vargs)...);
}

// variant consuming a juce::Component
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  juce::Component&            comp,
  args&&... vargs)
{
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    make_contiguous_range (&comp, 1),
    std::forward<args> (vargs)...);
}

// variant consuming a slider_ex
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  slider_ext&                 s,
  args&&... vargs)
{
  // Assuming slider on top.
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    s.slider,
    s.label,
    std::forward<args> (vargs)...);
}

// variant consuming a button
template <class T, class V, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  button_ext<V>&              b,
  args&&... vargs)
{
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    b.button,
    std::forward<args> (vargs)...);
}

// variant consuming a combobox
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>          area_remainder,
  juce::Rectangle<T>          column,
  float                       width,
  std::array<float, N> const& heights,
  grid_row_idx<N>             row,
  combobox_ext&               c,
  args&&... vargs)
{
  if (!c.prev_next_enabled()) {
    // no button handling
    return grid (
      area_remainder,
      column,
      width,
      heights,
      row,
      c.combo,
      std::forward<args> (vargs)...);
  }
  auto get_next_area = [&]() {
    if (row() == 0) {
      column = area_remainder.removeFromLeft ((int) width);
    }
    auto ret = column.removeFromTop ((int) heights[row()]);
    ++row;
    return ret;
  };

  auto area = get_next_area();
  // if one of these behaviors is undesired use the "juce::Component" overload.
  if (c.prev_next_grid_width_factor() == 0.f) {
    // consume a row for the buttons
    c.combo.setBounds (area);
    area      = get_next_area();
    auto prev = area.removeFromLeft ((int) (width / 2.));
    c.prev.setBounds (prev);
    c.next.setBounds (area);
  }
  else {
    // everything on the same row
    auto next
      = area.removeFromLeft ((int) (width * c.prev_next_grid_width_factor()));
    c.combo.setBounds (next);
    next = area.removeFromLeft ((int) (area.getWidth() / 2.));
    c.prev.setBounds (next);
    c.next.setBounds (area);
  }

  return grid (
    area_remainder, column, width, heights, row, std::forward<args> (vargs)...);
}

} // namespace detail
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// places a list of vargs objects on a grid. Each column is of fixed width. The
// height of the rows are specified on "heights". The elemenst are placed top
// to bottom, left to right.
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>   area,
  float                width,
  std::array<float, N> heights,
  args&&... vargs)
{
  return detail::grid (
    area,
    juce::Rectangle<T>(),
    width,
    heights,
    detail::grid_row_idx<N>(),
    std::forward<args> (vargs)...);
}

// TODO: this belongs somewhere else

// -----------------------------------------------------------------------------
} // namespace artv
