#pragma once

#include <array>
#include <cassert>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

static constexpr char slider_track_bg_from_center[]
  = "slider_track_bg_from_center";
struct ellipse {
  float rx, ry, rw;
};

static ellipse get_ellipse (float cx, float cy, float rad)
{
  return {cx - rad, cy - rad, rad * 2.0f};
}

static void draw_rotary_1 (
  juce::Graphics& g,
  int             x,
  int             y,
  int             width,
  int             height,
  float           pos,
  float const     start_angle,
  float const     end_angle,
  juce::Slider&   s)
{

  auto radius   = (float) juce::jmin (width / 2, height / 2) - 4.0f;
  auto centre_x = (float) x + (float) width * 0.5f;
  auto centre_y = (float) y + (float) height * 0.5f;

  ellipse out = get_ellipse (centre_x, centre_y, radius);

  // external
  auto bg = s.findColour (juce::Slider::backgroundColourId);
  g.setGradientFill (juce::ColourGradient::vertical (
    bg.brighter (0.20), out.ry, bg.darker (0.2), out.ry + out.rw));
  g.fillEllipse (out.rx, out.ry, out.rw, out.rw);

  // outline
  auto            thumb          = s.findColour (juce::Slider::thumbColourId);
  constexpr float outline_factor = 0.73f;
  ellipse outline = get_ellipse (centre_x, centre_y, radius * outline_factor);

  g.setGradientFill (juce::ColourGradient::vertical (
    thumb.brighter (0.40),
    outline.ry,
    thumb.darker (0.2),
    outline.ry + outline.rw));
  g.fillEllipse (outline.rx, outline.ry, outline.rw, outline.rw);

  // inside
  constexpr float in_factor = 0.63f;

  ellipse in = get_ellipse (centre_x, centre_y, radius * in_factor);
#if 1 // flashier: more cpu
  g.setGradientFill (juce::ColourGradient::vertical (
    thumb, in.ry, thumb.brighter (0.65), in.ry + in.rw));
#else
  g.setColour (thumb);
#endif
  g.fillEllipse (in.rx, in.ry, in.rw, in.rw);

  // pointer
  juce::Path p;
  auto       angle     = start_angle + pos * (end_angle - start_angle);
  auto       in_radius = in.rw * 0.5f;
  auto       p_length  = in_radius * 0.33f;
  auto       p_width   = p_length * 0.33f;
  p.addRectangle (-p_width * 0.5f, -in_radius, p_width, p_length);
  p.applyTransform (
    juce::AffineTransform::rotation (angle).translated (centre_x, centre_y));
  g.setColour (thumb.brighter (1.80));
  g.fillPath (p);

  bool center_track = s.getProperties().contains (slider_track_bg_from_center);

  // arc
  p.clear();
  constexpr auto arc_factor = 0.80f;
  p.addCentredArc (
    centre_x,
    centre_y,
    radius * arc_factor,
    radius * arc_factor,
    0.f,
    center_track ? start_angle + 0.5f * (end_angle - start_angle) : start_angle,
    angle,
    true);
  g.setColour (s.findColour (juce::Slider::trackColourId));
  g.strokePath (p, juce::PathStrokeType (2.f));
}

// Get rid of inheritance legacy
class saner_look_and_feel : public juce::LookAndFeel_V4 {
public:
  // This class is big, so add methods on an as-needed basis.
  // ---------------------------------------------------------------------------
  void drawRotarySlider (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           pos,
    float const     start_angle,
    float const     end_angle,
    juce::Slider&   s) override
  {
    if (on_rotary_draw) {
      on_rotary_draw (g, x, y, width, height, pos, start_angle, end_angle, s);
    }
    else {
      juce::LookAndFeel_V4::drawRotarySlider (
        g, x, y, width, height, pos, start_angle, end_angle, s);
    }
  }
  std::function<void (
    juce::Graphics&,
    int,
    int,
    int,
    int,
    float,
    float const,
    float const,
    juce::Slider&)>
    on_rotary_draw;
  // ---------------------------------------------------------------------------
  void drawLinearSlider (
    juce::Graphics&                 g,
    int                             x,
    int                             y,
    int                             width,
    int                             height,
    float                           pos,
    float                           min_pos,
    float                           max_pos,
    const juce::Slider::SliderStyle style,
    juce::Slider&                   s) override
  {
    if (on_linear_slider_draw) {
      on_linear_slider_draw (
        g, x, y, width, height, pos, min_pos, max_pos, style, s);
    }
    else {
      LookAndFeel_V4::drawLinearSlider (
        g, x, y, width, height, pos, min_pos, max_pos, style, s);
    }
  }
  std::function<void (
    juce::Graphics&,
    int,
    int,
    int,
    int,
    float,
    float,
    float,
    const juce::Slider::SliderStyle,
    juce::Slider&)>
    on_linear_slider_draw;
  // ---------------------------------------------------------------------------
  juce::Font getLabelFont (juce::Label& obj) override
  {
    if (on_get_label_font) {
      // avoid recursing, the original V4 function can be called on the callback
      auto fn           = std::move (on_get_label_font);
      auto ret          = fn (obj);
      on_get_label_font = std::move (fn);
      return ret;
    }
    return juce::LookAndFeel_V4::getLabelFont (obj);
  }
  std::function<juce::Font (juce::Label&)> on_get_label_font;
  // ---------------------------------------------------------------------------
  juce::Font getComboBoxFont (juce::ComboBox& obj) override
  {
    if (on_get_combobox_font) {
      // avoid recursing, the original V4 function can be called on the callback
      auto fn              = std::move (on_get_combobox_font);
      auto ret             = fn (obj);
      on_get_combobox_font = std::move (fn);
      return ret;
    }
    return juce::LookAndFeel_V4::getComboBoxFont (obj);
  }
  std::function<juce::Font (juce::ComboBox&)> on_get_combobox_font;
  // ---------------------------------------------------------------------------
  juce::Font getPopupMenuFont() override
  {
    if (on_get_popup_menu_font) {
      // avoid recursing, the original V4 function can be called on the callback
      auto fn                = std::move (on_get_popup_menu_font);
      auto ret               = fn();
      on_get_popup_menu_font = std::move (fn);
      return ret;
    }
    return juce::LookAndFeel_V4::getPopupMenuFont();
  }
  std::function<juce::Font()> on_get_popup_menu_font;
  // ---------------------------------------------------------------------------
  void fillTextEditorBackground (
    juce::Graphics&   g,
    int               width,
    int               height,
    juce::TextEditor& te) override
  {
    if (on_fill_text_editor_background) {
      on_fill_text_editor_background (g, width, height, te);
    }
    else {
      return juce::LookAndFeel_V4::fillTextEditorBackground (
        g, width, height, te);
    }
  }
  std::function<void (juce::Graphics&, int, int, juce::TextEditor&)>
    on_fill_text_editor_background;
  // ---------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------
} // namespace artv
