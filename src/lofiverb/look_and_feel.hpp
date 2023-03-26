#pragma once

#include "artv-common/juce/look_and_feel.hpp"

#include <array>
#include <cassert>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
//------------------------------------------------------------------------------
class look_and_feel : public juce::LookAndFeel_V4 {
public:
  look_and_feel()
  {
    using namespace juce;
    setColourScheme (getGreyColourScheme());
    setDefaultLookAndFeel (this);
  }

  ~look_and_feel() { setDefaultLookAndFeel (nullptr); }

  juce::Font getLabelFont (juce::Label& obj) override
  {
    auto f = LookAndFeel_V4::getLabelFont (obj);
    f.setHeight (obj.getHeight() - 4);
    return f;
  }

  juce::Font getSidePanelTitleFont (juce::SidePanel& obj) override
  {
    auto f = LookAndFeel_V4::getSidePanelTitleFont (obj);
    // f.setHeight (obj.getHeight() / 2);
    f.setHeight (5);
    return f;
  }

  void drawRotarySlider (
    juce::Graphics& g,
    int             x,
    int             y,
    int             width,
    int             height,
    float           sliderPos,
    float const     rotaryStartAngle,
    float const     rotaryEndAngle,
    juce::Slider&   s) override
  {
    draw_rotary_1 (
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
};

} // namespace artv
