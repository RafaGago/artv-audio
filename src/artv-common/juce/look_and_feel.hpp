#pragma once

#include <array>
#include <cassert>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/range.hpp"

namespace artv {

//------------------------------------------------------------------------------
class look_and_feel : public juce::LookAndFeel_V4 {
public:
  static constexpr char slider_track_bg_from_center[]
    = "slider_track_bg_from_center";

  look_and_feel()
  {
    using namespace juce;

    setColourScheme (getMidnightColourScheme());
    auto bg_color = findColour (ResizableWindow::backgroundColourId);
    // setColourScheme (getGreyColourScheme());
    setColour (TextButton::textColourOnId, bg_color);
    setColour (
      TextButton::buttonOnColourId, juce::Colours::teal.brighter (0.6));
    setColour (TextButton::textColourOffId, juce::Colours::white);
    setColour (TextButton::buttonColourId, bg_color);
    //
    //        setColour (ResizableWindow::backgroundColourId, Colour
    //        (0xfffcfcf7)); setColour (Slider::backgroundColourId, Colour
    //        (0xfffcfdf7)); setColour(
    //            Slider::textBoxOutlineColourId,
    //            findColour (ResizableWindow::backgroundColourId)
    //            );
    // notice:
    // juce::Array<juce::Font> fonts;
    // juce::Font::findFonts (fonts);
    // fonts[x].font.getTypefaceName()

    // setDefaultSansSerifTypefaceName ("TODO");
    setDefaultLookAndFeel (this);
  }

  ~look_and_feel() { setDefaultLookAndFeel (nullptr); }

  juce::Font getLabelFont (juce::Label& obj) override
  {
    auto f = LookAndFeel_V4::getLabelFont (obj);
    f.setHeight (obj.getHeight() - 4);
    return f;
  }

#if 0
    juce::Font getComboBoxFont (juce::ComboBox& obj) override
    {
        auto f = LookAndFeel_V4::getComboBoxFont (obj);
        f.setHeight (obj.getHeight() / 2);
        return f;
    }


    juce::Font getPopupMenuFont () override
    {
        auto f = LookAndFeel_V4::getPopupMenuFont();
        f.setHeight (f.getHeight() / 2);
        return f;
    }
#endif
  juce::Font getSidePanelTitleFont (juce::SidePanel& obj) override
  {
    auto f = LookAndFeel_V4::getSidePanelTitleFont (obj);
    // f.setHeight (obj.getHeight() / 2);
    f.setHeight (5);
    return f;
  }

  struct ellipse {
    float rx, ry, rw;
  };

  ellipse get_ellipse (float cx, float cy, float rad)
  {
    return {cx - rad, cy - rad, rad * 2.0f};
  }

  void drawRotarySlider (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           sliderPos,
    const float     rotaryStartAngle,
    const float     rotaryEndAngle,
    juce::Slider&   s) override
  {
    drawDefaultRotary (
      g, x, y, width, height, sliderPos, rotaryStartAngle, rotaryEndAngle, s);
  }

  void drawLinearSlider (
    juce::Graphics&                 g,
    int                             x,
    int                             y,
    int                             width,
    int                             height,
    float                           sliderPos,
    float                           minSliderPos,
    float                           maxSliderPos,
    const juce::Slider::SliderStyle style,
    juce::Slider&                   slider) override
  {
    using namespace juce;

    // TODO: implement for real.
    juce::Colour track_color;
    bool         center_track
      = slider.getProperties().contains (slider_track_bg_from_center);

    if (center_track) {
      track_color = slider.findColour (Slider::trackColourId);
      slider.setColour (
        Slider::trackColourId, slider.findColour (Slider::backgroundColourId));
    }

    LookAndFeel_V4::drawLinearSlider (
      g,
      x,
      y,
      width,
      height,
      sliderPos,
      minSliderPos,
      maxSliderPos,
      style,
      slider);

    if (center_track) {
      slider.setColour (Slider::trackColourId, track_color);
    }
  }

private:
  void drawDefaultRotary (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           sliderPos,
    const float     rotaryStartAngle,
    const float     rotaryEndAngle,
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
    auto       angle
      = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
    auto in_radius = in.rw * 0.5f;
    auto p_length  = in_radius * 0.33f;
    auto p_width   = p_length * 0.33f;
    p.addRectangle (-p_width * 0.5f, -in_radius, p_width, p_length);
    p.applyTransform (
      juce::AffineTransform::rotation (angle).translated (centre_x, centre_y));
    g.setColour (thumb.brighter (1.80));
    g.fillPath (p);

    // arc
    constexpr auto arc_factor = 0.80f;
    p.addCentredArc (
      centre_x,
      centre_y,
      radius * arc_factor,
      radius * arc_factor,
      0.f,
      rotaryStartAngle,
      angle,
      true);
    g.strokePath (p, juce::PathStrokeType (2.f));
  }

  // needs more work...
  void drawGradientBarRotary (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           sliderPos,
    const float     rotaryStartAngle,
    const float     rotaryEndAngle,
    juce::Slider&   s)
  {

    // auto radius   = (float) juce::jmin (width / 2, height / 2) - 4.0f;
    auto centre_x = (float) x + (float) width * 0.5f;
    auto centre_y = (float) y + (float) height * 0.5f;

    // const auto pi = juce::MathConstants<float>::pi;

    juce::ColourGradient gradient (
      juce::Colours::purple, x, y, juce::Colours::yellow, x + width, y, false);
    gradient.addColour (0.125, juce::Colours::purple);
    gradient.addColour (0.875, juce::Colours::yellow);
    g.setGradientFill (gradient);

    auto angle
      = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
    juce::Path p;
    p.addCentredArc (
      centre_x, centre_y, 35.f, 35.f, 0.f, rotaryStartAngle, angle, true);
    g.strokePath (p, juce::PathStrokeType (2.f));
  }

  void drawRotaryVariation1 (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           sliderPos,
    const float     rotaryStartAngle,
    const float     rotaryEndAngle,
    juce::Slider&   s)
  {
    auto radius   = (float) juce::jmin (width / 2, height / 2) - 4.0f;
    auto centre_x = (float) x + (float) width * 0.5f;
    auto centre_y = (float) y + (float) height * 0.5f;

    ellipse out = get_ellipse (centre_x, centre_y, radius);
    // external
    auto bg  = s.findColour (juce::Slider::backgroundColourId);
    auto cbg = bg.brighter (0.05);
    g.setGradientFill (
      juce::ColourGradient::vertical (cbg, out.ry, bg, out.ry + out.rw));
    g.fillEllipse (out.rx, out.ry, out.rw, out.rw);

    // outline
    auto thumb     = s.findColour (juce::Slider::thumbColourId);
    auto outline_h = thumb.brighter (0.31);
    auto outline_l = outline_h.brighter (0.41);

    constexpr float outline_factor = 0.70f;
    ellipse outline = get_ellipse (centre_x, centre_y, radius * outline_factor);

    g.setGradientFill (juce::ColourGradient::vertical (
      outline_h, outline.ry, outline_l, outline.ry + outline.rw));
    g.fillEllipse (outline.rx, outline.ry, outline.rw, outline.rw);

    // inside
    auto cthumb = thumb.brighter (0.37);

    constexpr float in_factor = 0.61f;
    ellipse         in = get_ellipse (centre_x, centre_y, radius * in_factor);
    g.setGradientFill (
      juce::ColourGradient::vertical (thumb, in.ry, cthumb, in.ry + in.rw));
    g.fillEllipse (in.rx, in.ry, in.rw, in.rw);

    // pointer
    juce::Path p;
    auto       angle
      = rotaryStartAngle + sliderPos * (rotaryEndAngle - rotaryStartAngle);
    auto in_radius = in.rw * 0.5f;
    auto p_length  = in_radius * 0.33f;
    auto p_width   = p_length * 0.33f;
    p.addRectangle (-p_width * 0.5f, -in_radius, p_width, p_length);
    p.applyTransform (
      juce::AffineTransform::rotation (angle).translated (centre_x, centre_y));
    g.setColour (thumb.brighter (1.03));
    // TODO: a lo mejor g.setFillType()
    g.fillPath (p);
  }
};

// -----------------------------------------------------------------------------
} // namespace artv
