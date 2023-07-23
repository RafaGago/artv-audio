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
private:
  template <class F, class... Ts>
  auto invoke (F& func, Ts&&... args)
  {
    // allow calling the original v4 functions from the callback
    auto f = std::move (func);
    if constexpr (std::is_same_v<typename F::result_type, void>) {
      f (std::forward<Ts> (args)...);
      func = std::move (f);
    }
    else {
      auto ret = f (std::forward<Ts> (args)...);
      func     = std::move (f);
      return ret;
    }
  }

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
      invoke (
        on_rotary_draw, g, x, y, width, height, pos, start_angle, end_angle, s);
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
      invoke (
        on_linear_slider_draw,
        g,
        x,
        y,
        width,
        height,
        pos,
        min_pos,
        max_pos,
        style,
        s);
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
      return invoke (on_get_label_font, obj);
    }
    return juce::LookAndFeel_V4::getLabelFont (obj);
  }
  std::function<juce::Font (juce::Label&)> on_get_label_font;
  // ---------------------------------------------------------------------------
  juce::Font getComboBoxFont (juce::ComboBox& obj) override
  {
    if (on_get_combobox_font) {
      return invoke (on_get_combobox_font, obj);
    }
    return juce::LookAndFeel_V4::getComboBoxFont (obj);
  }
  std::function<juce::Font (juce::ComboBox&)> on_get_combobox_font;
  // ---------------------------------------------------------------------------
  juce::Font getPopupMenuFont() override
  {
    if (on_get_popup_menu_font) {
      return invoke (on_get_popup_menu_font);
    }
    return juce::LookAndFeel_V4::getPopupMenuFont();
  }
  std::function<juce::Font()> on_get_popup_menu_font;
  // ---------------------------------------------------------------------------
  juce::Font getSidePanelTitleFont (juce::SidePanel& obj) override
  {
    if (on_get_side_panel_title_font) {
      return invoke (on_get_side_panel_title_font, obj);
    }
    return juce::LookAndFeel_V4::getSidePanelTitleFont (obj);
  }
  std::function<juce::Font (juce::SidePanel& obj)> on_get_side_panel_title_font;
  // ---------------------------------------------------------------------------
  void fillTextEditorBackground (
    juce::Graphics&   g,
    int               width,
    int               height,
    juce::TextEditor& te) override
  {
    if (on_fill_text_editor_background) {
      invoke (on_fill_text_editor_background, g, width, height, te);
    }
    else {
      return juce::LookAndFeel_V4::fillTextEditorBackground (
        g, width, height, te);
    }
  }
  std::function<void (juce::Graphics&, int, int, juce::TextEditor&)>
    on_fill_text_editor_background;
  // ---------------------------------------------------------------------------
  juce::PopupMenu::Options getOptionsForComboBoxPopupMenu (
    juce::ComboBox& cb,
    juce::Label&    l) override
  {
    if (on_get_options_for_combobox_popup_menu) {
      return invoke (on_get_options_for_combobox_popup_menu, cb, l);
    }
    else {
      return juce::LookAndFeel_V4::getOptionsForComboBoxPopupMenu (cb, l);
    }
  }
  std::function<juce::PopupMenu::Options (juce::ComboBox&, juce::Label&)>
    on_get_options_for_combobox_popup_menu;
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // override functions that are as of now not very interesting to expose
  // ---------------------------------------------------------------------------
  void drawPopupMenuItem (
    juce::Graphics&             g,
    juce::Rectangle<int> const& area,
    bool const                  isSeparator,
    bool const                  isActive,
    bool const                  isHighlighted,
    bool const                  isTicked,
    bool const                  hasSubMenu,
    juce::String const&         text,
    juce::String const&         shortcutKeyText,
    juce::Drawable const*       icon,
    juce::Colour const* const   textColourToUse) override
  {
    using namespace juce;
    if (isSeparator) {
      auto r = area.reduced (5, 0);
      r.removeFromTop (roundToInt (((float) r.getHeight() * 0.5f) - 0.5f));

      g.setColour (findColour (PopupMenu::textColourId).withAlpha (0.3f));
      g.fillRect (r.removeFromTop (1));
    }
    else {
      auto textColour
        = (textColourToUse == nullptr ? findColour (PopupMenu::textColourId) : *textColourToUse);

      auto r = area.reduced (1);

      if (isHighlighted && isActive) {
        g.setColour (findColour (PopupMenu::highlightedBackgroundColourId));
        g.fillRect (r);

        g.setColour (findColour (PopupMenu::highlightedTextColourId));
      }
      else {
        g.setColour (textColour.withMultipliedAlpha (isActive ? 1.0f : 0.5f));
      }

      r.reduce (jmin (5, area.getWidth() / 20), 0);

      auto font = getPopupMenuFont();
#if LF_ORIGINAL_CODE
      auto maxFontHeight = (float) r.getHeight() / 1.3f;

      if (font.getHeight() > maxFontHeight)
        font.setHeight (maxFontHeight);
#else
      auto maxFontHeight = (float) r.getHeight();
#endif

      g.setFont (font);

      auto iconArea = r.removeFromLeft (roundToInt (maxFontHeight)).toFloat();

      if (icon != nullptr) {
        icon->drawWithin (
          g,
          iconArea,
          RectanglePlacement::centred | RectanglePlacement::onlyReduceInSize,
          1.0f);
        r.removeFromLeft (roundToInt (maxFontHeight * 0.5f));
      }
      else if (isTicked) {
        auto tick = getTickShape (1.0f);
        g.fillPath (
          tick,
          tick.getTransformToScaleToFit (
            iconArea.reduced (iconArea.getWidth() / 5, 0).toFloat(), true));
      }

      if (hasSubMenu) {
        auto arrowH = 0.6f * getPopupMenuFont().getAscent();

        auto x = static_cast<float> (r.removeFromRight ((int) arrowH).getX());
        auto halfH = static_cast<float> (r.getCentreY());

        Path path;
        path.startNewSubPath (x, halfH - arrowH * 0.5f);
        path.lineTo (x + arrowH * 0.6f, halfH);
        path.lineTo (x, halfH + arrowH * 0.5f);

        g.strokePath (path, PathStrokeType (2.0f));
      }

      r.removeFromRight (3);
      g.drawFittedText (text, r, Justification::centredLeft, 1);

      if (shortcutKeyText.isNotEmpty()) {
        auto f2 = font;
        f2.setHeight (f2.getHeight() * 0.75f);
        f2.setHorizontalScale (0.95f);
        g.setFont (f2);

        g.drawText (shortcutKeyText, r, Justification::centredRight, true);
      }
    }
  }
  // ---------------------------------------------------------------------------
  // a copy from lookandfeel v4
  void getIdealPopupMenuItemSize (
    juce::String const& text,
    bool const          isSeparator,
    int                 standardMenuItemHeight,
    int&                idealWidth,
    int&                idealHeight) override
  {
    using namespace juce;
    if (isSeparator) {
      idealWidth = 50;
      idealHeight
        = standardMenuItemHeight > 0 ? standardMenuItemHeight / 10 : 10;
    }
    else {
      auto font = getPopupMenuFont();
#if LF_ORIGINAL_CODE
      if (
        standardMenuItemHeight > 0
        && font.getHeight() > (float) standardMenuItemHeight / 1.3f)
        font.setHeight ((float) standardMenuItemHeight / 1.3f);

      idealHeight = standardMenuItemHeight > 0
        ? standardMenuItemHeight
        : roundToInt (font.getHeight() * 1.3f);
#else
      idealHeight        = font.getHeight() * 1.1f;
#endif
      idealWidth = font.getStringWidth (text) + idealHeight * 2;
    }
  }
};
// -----------------------------------------------------------------------------
} // namespace artv
