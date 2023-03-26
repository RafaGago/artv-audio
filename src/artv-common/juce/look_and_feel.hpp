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
  float           sliderPos,
  float const     rotaryStartAngle,
  float const     rotaryEndAngle,
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
    center_track ? rotaryStartAngle + 0.5f * (rotaryEndAngle - rotaryStartAngle)
                 : rotaryStartAngle,
    angle,
    true);
  g.setColour (s.findColour (juce::Slider::trackColourId));
  g.strokePath (p, juce::PathStrokeType (2.f));
}

// -----------------------------------------------------------------------------
} // namespace artv
