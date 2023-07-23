#pragma once

#include <array>
#include <cassert>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <type_traits>

#include "artv-common/juce/widgets.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

// -----------------------------------------------------------------------------
template <typename... args>
static void set_color (int id, juce::Colour color, args&&... components)
{
  artv::apply (
    [=] (auto& component) {
      if constexpr (std::is_pointer_v<
                      std::remove_reference_t<decltype (component)>>) {
        component->setColour (id, color);
      }
      else {
        component.setColour (id, color);
      }
    },
    std::forward<args> (components)...);
}
// -----------------------------------------------------------------------------
namespace detail {

template <class T>
struct label_on_top {
  T* ref;
};

template <class T>
struct no_buttons {
  T* ref;
};

template <class T>
struct buttons_on_same_row {
  T*    ref;
  float combo_width_factor;
};

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
  xspan<C>                    pack,
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
    xspan (arr),
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
    xspan {&comp, 1},
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

template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>              area_remainder,
  juce::Rectangle<T>              column,
  float                           width,
  std::array<float, N> const&     heights,
  grid_row_idx<N>                 row,
  label_on_top<slider_ext> const& s,
  args&&... vargs)
{
  // Assuming slider on top.
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    s.ref->label,
    s.ref->slider,
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
  auto get_next_area = [&]() {
    if (row() == 0) {
      column = area_remainder.removeFromLeft ((int) width);
    }
    auto ret = column.removeFromTop ((int) heights[row()]);
    ++row;
    return ret;
  };

  auto area = get_next_area();
  // consume a row for the buttons
  c.combo.setBounds (area);
  area      = get_next_area();
  auto prev = area.removeFromLeft ((int) (width / 2.));
  c.prev.setBounds (prev);
  c.next.setBounds (area);

  return grid (
    area_remainder, column, width, heights, row, std::forward<args> (vargs)...);
}

// variant consuming a combobox
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>              area_remainder,
  juce::Rectangle<T>              column,
  float                           width,
  std::array<float, N> const&     heights,
  grid_row_idx<N>                 row,
  no_buttons<combobox_ext> const& c,
  args&&... vargs)
{
  // no button handling
  return grid (
    area_remainder,
    column,
    width,
    heights,
    row,
    c.ref->combo,
    std::forward<args> (vargs)...);
}

// variant consuming a combobox
template <class T, size_t N, class... args>
static juce::Rectangle<T> grid (
  juce::Rectangle<T>                       area_remainder,
  juce::Rectangle<T>                       column,
  float                                    width,
  std::array<float, N> const&              heights,
  grid_row_idx<N>                          row,
  buttons_on_same_row<combobox_ext> const& c,
  args&&... vargs)
{
  float wf        = c.combo_width_factor;
  int   prev_conn = juce::TextButton::ConnectedOnRight;
  if (wf <= 0.f || wf >= 1.f) {
    wf = 0.f;
  }
  else {
    prev_conn |= juce::TextButton::ConnectedOnLeft;
  }
  c.ref->prev.setConnectedEdges (prev_conn);

  auto get_next_area = [&]() {
    if (row() == 0) {
      column = area_remainder.removeFromLeft ((int) width);
    }
    auto ret = column.removeFromTop ((int) heights[row()]);
    ++row;
    return ret;
  };

  auto area = get_next_area();
  // combo can't set connected edges...

  // everything on the same row
  auto next = area.removeFromLeft ((int) (width * wf));
  c.ref->combo.setBounds (next);
  next = area.removeFromLeft ((int) (area.getWidth() / 2.));
  c.ref->prev.setBounds (next);
  c.ref->next.setBounds (area);

  return grid (
    area_remainder, column, width, heights, row, std::forward<args> (vargs)...);
}

} // namespace detail
// -----------------------------------------------------------------------------
inline detail::label_on_top<slider_ext> grid_label_on_top (slider_ext& v)
{
  return {&v};
}
// -----------------------------------------------------------------------------
inline detail::no_buttons<combobox_ext> grid_no_buttons (combobox_ext& v)
{
  return {&v};
}
// -----------------------------------------------------------------------------
inline detail::buttons_on_same_row<combobox_ext> grid_buttons_on_same_row (
  combobox_ext& v,
  float         combo_width_ratio)
{
  assert (combo_width_ratio >= 0 && combo_width_ratio <= 1.f);
  return {&v, combo_width_ratio};
}
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
