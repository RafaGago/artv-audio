#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <stdint.h>
#include <utility>

#include <boost/hana.hpp>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/gui_util.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/parameters.hpp"
#include "artv-common/juce/value_tree_attachments.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "lofiverb/look_and_feel.hpp"
#include "lofiverb/parameters.hpp"

#define VERSION_INT VERSION_GET (VERSION_MAJOR, VERSION_MINOR, VERSION_REV)

namespace artv {

// -----------------------------------------------------------------------------
struct vertical_line : public juce::Component {
  //----------------------------------------------------------------------------
  bool left = true;
  //----------------------------------------------------------------------------
  void paint (juce::Graphics& g) override
  {
    auto c = findColour (juce::Slider::backgroundColourId).brighter (0.5);
    g.setColour (c);
    auto lb = getLocalBounds();

    float x = lb.getX();
    if (!left) {
      x += lb.getWidth();
    }
    juce::Line<float> line {
      juce::Point<float> (x, lb.getY() + 4.f),
      juce::Point<float> (x, lb.getY() + lb.getHeight() - 4.f)};

    g.drawLine (line, 2.f);
  }
  //----------------------------------------------------------------------------
};
// -----------------------------------------------------------------------------
struct panel : public juce::Component {
  // reuse, so we don't clash...
  enum { backgroundColourId = juce::ResizableWindow::backgroundColourId };
  //----------------------------------------------------------------------------
  bool upper_line = true;
  //----------------------------------------------------------------------------
  void paint (juce::Graphics& g) override
  {
    auto bg = findColour (backgroundColourId);
    g.setColour (bg);
    auto lb = getLocalBounds();
    g.fillRect (lb);

    if (!upper_line) {
      return;
    }
    auto line = lb.removeFromTop (2);
    g.setColour (bg.brighter (0.3));
    g.fillRect (line);
  }
  //----------------------------------------------------------------------------
};
// -----------------------------------------------------------------------------
class editor : public juce::AudioProcessorEditor,
               public juce::DragAndDropContainer {
public:
  //----------------------------------------------------------------------------
  explicit editor (
    juce::AudioProcessor&               p,
    juce::AudioProcessorValueTreeState& params,
    juce::ValueTree&)
    : AudioProcessorEditor (p), _processor (p)
  {
    //        id, juce::Colours::green.brighter (0.2),
    //        juce::Colours::orange.brighter (0.07),
    //        juce::Colours::teal,
    //        juce::Colours::darkgoldenrod,
    //        juce::Colours::red.brighter (0.4),
    //        juce::Colours::sienna,

    static constexpr uint up    = juce::TextButton::ConnectedOnTop;
    static constexpr uint down  = juce::TextButton::ConnectedOnBottom;
    static constexpr uint left  = juce::TextButton::ConnectedOnLeft;
    static constexpr uint right = juce::TextButton::ConnectedOnRight;

    //_header.setLookAndFeel (&_lookfeel);
    //_display_frame.setLookAndFeel (&_lookfeel);
    //_display_value.setLookAndFeel (&_lookfeel);

    // setLookAndFeel (&_lookfeel);
    _params.init_widgets (*this, params);

    register_mouse_events();

    // set slider styles
    auto& dry          = *_params.p_get (parameters::dry {})[0];
    auto& wet          = *_params.p_get (parameters::wet {})[0];
    auto& wet_pan      = *_params.p_get (parameters::wet_pan {})[0];
    auto& algorithm    = *_params.p_get (parameters::algorithm {})[0];
    auto& predelay     = *_params.p_get (parameters::predelay {})[0];
    auto& decay        = *_params.p_get (parameters::decay {})[0];
    auto& damp         = *_params.p_get (parameters::damp {})[0];
    auto& stereo       = *_params.p_get (parameters::stereo {})[0];
    auto& character    = *_params.p_get (parameters::character {})[0];
    auto& freq_balance = *_params.p_get (parameters::freq_balance {})[0];
    auto& ducking_threshold
      = *_params.p_get (parameters::ducking_threshold {})[0];
    auto& ducking_speed   = *_params.p_get (parameters::ducking_speed {})[0];
    auto& mode            = *_params.p_get (parameters::mode {})[0];
    auto& mod             = *_params.p_get (parameters::mod {})[0];
    auto& operating_range = *_params.p_get (parameters::operating_range {})[0];

    dry.slider.setSliderStyle (juce::Slider::LinearVertical);
    wet.slider.setSliderStyle (juce::Slider::LinearVertical);

    addAndMakeVisible (_header);
    addAndMakeVisible (_display_frame);
    addAndMakeVisible (_display_value);

    _display_value.setJustificationType (juce::Justification::left);
    //_display_value.setFont (juce::Font {
    //  juce::Font::getDefaultMonospacedFontName(), 12, juce::Font::bold});

    auto opaque      = juce::Colour {0xff000000};
    auto transparent = juce::Colour {0};
    auto dark_grey   = juce::Colour::fromString {0xff000000 | 0x121f1f};
    auto light_grey  = juce::Colour::fromString {0xff000000 | 0x1c2b2b};
    auto panel       = juce::Colour::fromString {0xff000000 | 0x487576};

    set_color (juce::TextButton::buttonColourId, panel, _display_frame);
    set_color (juce::Label::textColourId, panel.darker (0.7f), _display_value);
    set_color (juce::Label::backgroundColourId, transparent, _display_value);
    set_color (panel::backgroundColourId, dark_grey, _header);

    // colors
    auto color = make_array<juce::Colour> (
      juce::Colours::teal,
      juce::Colours::orange.brighter (0.07),
      juce::Colours::red.brighter (0.4),
      juce::Colour (0xff00cc99),
      juce::Colour (0xff0099cc), // 5
      juce::Colour (0xffff794d),
      juce::Colour (0xffd22d6f),
      juce::Colour (0xff9cb946));

    dry.slider.setColour (juce::Slider::thumbColourId, color[0]);
    wet.slider.setColour (juce::Slider::thumbColourId, color[0]);

    predelay.slider.setColour (juce::Slider::thumbColourId, color[0]);
    operating_range.slider.setColour (juce::Slider::thumbColourId, color[0]);

    decay.slider.setColour (juce::Slider::thumbColourId, color[1]);

    character.slider.setColour (juce::Slider::thumbColourId, color[0]);
    mod.slider.setColour (juce::Slider::thumbColourId, color[0]);

    freq_balance.slider.setColour (juce::Slider::thumbColourId, color[0]);
    damp.slider.setColour (juce::Slider::thumbColourId, color[0]);

    wet_pan.slider.setColour (juce::Slider::thumbColourId, color[0]);
    stereo.slider.setColour (juce::Slider::thumbColourId, color[0]);

    ducking_threshold.slider.setColour (juce::Slider::thumbColourId, color[0]);
    ducking_speed.slider.setColour (juce::Slider::thumbColourId, color[0]);

    // pan and stereo from center
    wet_pan.slider.getProperties().set (slider_track_bg_from_center, true);
    stereo.slider.getProperties().set (slider_track_bg_from_center, true);

    // size
    constexpr float ratio  = sizes::total_w_divs / sizes::total_h_divs;
    constexpr float factor = 0.25f;
    setResizable (true, true);
    getConstrainer()->setFixedAspectRatio (ratio);
    auto area = juce::Desktop::getInstance()
                  .getDisplays()
                  .getDisplayForPoint ({0, 0})
                  ->userArea;
    float height = factor * (float) area.getHeight();
    float width  = height * ratio;
    if (width > area.getWidth()) {
      width  = factor * (float) area.getWidth();
      height = width / ratio;
    }
    setSize (width, height);
  }
  //----------------------------------------------------------------------------
  ~editor() override
  {
    // this lookanfeel thing could really be a bit smarter...
    _params.clear_widgets();
    _display_value.setLookAndFeel (nullptr);
    _display_frame.setLookAndFeel (nullptr);
    setLookAndFeel (nullptr);
  }
  //----------------------------------------------------------------------------
  struct sizes {

    static constexpr float main_knob_wh_divs = 4.f;
    static constexpr float main_label_w_divs = main_knob_wh_divs;
    static constexpr float main_label_h_divs = main_knob_wh_divs / 4.f;

    static constexpr float w_separation_divs = 1.f;
    static constexpr float h_separation_divs = 1.f;

    static constexpr float header_h_divs     = 0.75f * main_knob_wh_divs;
    static constexpr float header_sep_h_divs = 1.f;

    static constexpr float header_margin_w_divs  = 4.f;
    static constexpr float header_display_w_divs = 10.f;

    static constexpr float header_w_divs
      = header_margin_w_divs + header_display_w_divs + header_display_w_divs;

    static constexpr float main_w_separation_divs = 1.f;
    static constexpr float main_w_divs            = main_knob_wh_divs;
    static constexpr float main_h_divs
      = 2 * main_knob_wh_divs + 2 * main_label_h_divs;

    static constexpr float drywet_w_divs = main_w_divs * 1.5f;
    static constexpr float drywet_h_divs = main_h_divs;

    static constexpr float decay_area_w_divs     = main_w_divs * 1.75f;
    static constexpr float decay_area_h_divs     = main_h_divs;
    static constexpr float decay_combobox_h_divs = main_label_h_divs;
    static constexpr float decay_label_h_divs    = main_label_h_divs;
    static constexpr float decay_knob_h_divs
      = main_h_divs - 2 * decay_combobox_h_divs - decay_label_h_divs;

    static constexpr float total_w_divs = //
      w_separation_divs + //
      drywet_w_divs + //
      w_separation_divs + //
      main_w_divs + // predelay + oprange
      w_separation_divs + //
      decay_area_w_divs + //
      w_separation_divs + //
      main_w_divs + // charac + mod
      w_separation_divs + //
      main_w_divs + // freq
      w_separation_divs + //
      main_w_divs + // stereo
      w_separation_divs + //
      main_w_divs + // ducking
      w_separation_divs //
      ;
    static constexpr float total_h_divs = //
      h_separation_divs + //
      header_h_divs + //
      h_separation_divs + //
      main_h_divs + //
      h_separation_divs;
    ;
  };
  //----------------------------------------------------------------------------
  void resized() override
  {
    auto area = getLocalBounds();

    auto h_units = (float) area.getHeight();
    auto w_units = (float) area.getWidth();

    auto h = (float) h_units / sizes::total_h_divs;
    auto w = (float) w_units / sizes::total_w_divs;

    auto sep_w = sizes::w_separation_divs * w;
    auto sep_h = sizes::h_separation_divs * h;

    auto label_h = sizes::main_label_h_divs * h;
    auto knob_h  = sizes::main_knob_wh_divs * h;

    // header
    float h_h             = sizes::header_h_divs * h;
    auto  h_w             = (float) w_units / sizes::header_w_divs;
    auto  header          = area.removeFromTop (sep_h + h_h + sep_h);
    auto  header_margin_w = h_w * sizes::header_margin_w_divs;
    _header.setBounds (header);

    header.removeFromTop (sep_h); //  upper margin
    header.removeFromBottom (sep_h); //  lower margin

    header.removeFromLeft (header_margin_w);
    header.removeFromRight (header_margin_w);

    auto display = header;

    _display_frame.setBounds (display);
    auto pvfbounds = _display_frame.getBounds().reduced (h_h / 10.f);
    // make the font smaller to avoid resizings
    pvfbounds.reduce (0, pvfbounds.getHeight() / 9);
    _display_value.setBounds (pvfbounds);

    area.removeFromTop (sep_h); //  header margin
    auto main = area.removeFromTop (h * (float) sizes::main_h_divs);

    main.removeFromLeft (sep_w);
    auto drywet = main.removeFromLeft (w * (float) sizes::drywet_w_divs);

    grid (
      drywet,
      ((float) drywet.getWidth()) / 2.f,
      make_array (label_h, label_h + 2 * knob_h),
      grid_label_on_top (*_params.p_get (parameters::dry {})[0]),
      grid_label_on_top (*_params.p_get (parameters::wet {})[0]));

    main.removeFromLeft (sep_w);
    auto predel = main.removeFromLeft (w * (float) sizes::main_w_divs);

    grid (
      predel,
      predel.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::predelay {})[0]),
      grid_label_on_top (*_params.p_get (parameters::operating_range {})[0]));

    main.removeFromLeft (sep_w);
    auto decay = main.removeFromLeft (w * (float) sizes::decay_area_w_divs);

    auto decay_combobox_h = sizes::decay_combobox_h_divs * h;
    auto decay_knob_h     = sizes::decay_knob_h_divs * h;
    auto decay_label_h    = sizes::decay_label_h_divs * h;

    grid (
      decay,
      decay.getWidth(),
      make_array (
        decay_label_h, decay_knob_h, decay_combobox_h, decay_combobox_h),
      grid_label_on_top (*_params.p_get (parameters::decay {})[0]),
      grid_buttons_on_same_row (
        *_params.p_get (parameters::algorithm {})[0], 0.8f),
      grid_buttons_on_same_row (*_params.p_get (parameters::mode {})[0], 0.8f));

    main.removeFromLeft (sep_w);
    auto mod = main.removeFromLeft (w * (float) sizes::main_w_divs);

    grid (
      mod,
      mod.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::mod {})[0]),
      grid_label_on_top (*_params.p_get (parameters::character {})[0]));

    main.removeFromLeft (sep_w);
    auto freq = main.removeFromLeft (w * (float) sizes::main_w_divs);

    grid (
      freq,
      freq.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::freq_balance {})[0]),
      grid_label_on_top (*_params.p_get (parameters::damp {})[0]));

    main.removeFromLeft (sep_w);
    auto stereo = main.removeFromLeft (w * (float) sizes::main_w_divs);

    grid (
      stereo,
      stereo.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::stereo {})[0]),
      grid_label_on_top (*_params.p_get (parameters::wet_pan {})[0]));

    main.removeFromLeft (sep_w);
    auto duck = main.removeFromLeft (w * (float) sizes::main_w_divs);

    grid (
      duck,
      duck.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::ducking_threshold {})[0]),
      grid_label_on_top (*_params.p_get (parameters::ducking_speed {})[0]));

    // area.removeFromTop (sep_h); //  bottom margin
  }
//----------------------------------------------------------------------------
#if 0
  void paint (juce::Graphics& g) override
  {
    g.fillAll (
      getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
  }
#endif
  //----------------------------------------------------------------------------
  void mouse_event (
    event_common::source    src,
    mouse_event::type       event_type,
    juce::MouseEvent const& event,
    mouse_event::data       data)
  {
    juce::String           text;
    juce::Component const* component;

    if (std::holds_alternative<juce::Slider const*> (src)) {
      auto& slider = *std::get<juce::Slider const*> (src);
      component    = &slider;

      juce::String value = const_cast<juce::Slider&> (slider).getTextFromValue (
        slider.getValue());

      text = slider.getName();

      auto suffix = slider.getTextValueSuffix();
      if (suffix.length() > 0) {
        text += " (" + suffix + ")";
      }
      text += ": " + value.removeCharacters (suffix);
      text = text.trim();
    }
    else if (std::holds_alternative<juce::Button const*> (src)) {
      auto& button = *std::get<juce::Button const*> (src);
      component    = &button;

      text = button.getName();
      if (!text.endsWith (" Next") && !text.endsWith (" Prev")) {
        text += ": ";
        text += button.getToggleState() ? "On" : "Off";
      }
    }
    else if (std::holds_alternative<juce::ComboBox const*> (src)) {
      auto& combobox = *std::get<juce::ComboBox const*> (src);
      component      = &combobox;

      text = combobox.getName();
      text += ": ";
      text += combobox.getText();
    }
    else {
      // what?
      return;
    }
    if (event_type == mouse_event::exit) {
      text.clear();
    }
    if (!component->isEnabled()) {
      text = "Disabled";
    }
    _display_value.setText (text, juce::NotificationType::dontSendNotification);
  }
  //----------------------------------------------------------------------------
  void register_mouse_events()
  {
    using namespace std::placeholders;

    _params.pforeach ([this] (auto key, auto& warray) {
      using type = typename std::remove_reference<decltype (*warray[0])>::type;

      for (auto& w : warray) {
        if constexpr (is_slider_ext<type>::value) {
          w->slider.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
        }
        if constexpr (is_button_ext<type>::value) {
          w->button.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
        }
        if constexpr (is_toggle_buttons<type>::value) {
          for (auto& btn : w->buttons) {
            btn.on_mouse_event
              = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
          }
        }
        if constexpr (is_combobox_ext<type>::value) {
          w->combo.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
          w->prev.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
          w->next.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
        }
      }
    });
  }
  //----------------------------------------------------------------------------
private:
  juce::AudioProcessor& _processor; // unused
  // look_and_feel         _lookfeel;

  editor_apvts_widgets<parameters::parameters_typelist> _params;

  panel                _header;
  std::array<panel, 8> _main;
  juce::Label          _display_value;
  juce::TextButton     _display_frame;
  //  std::array<vertical_line, 2>                   _side_lines;
  //  std::array<vertical_line, n_stereo_busses * 2> _fx_lines;
}; // namespace artv
//------------------------------------------------------------------------------
juce::AudioProcessorEditor* new_editor (
  juce::AudioProcessor&               p,
  juce::AudioProcessorValueTreeState& params,
  juce::ValueTree&                    gui_params)
{
  return new artv::editor (p, params, gui_params);
}
// -----------------------------------------------------------------------------
} // namespace artv
// -----------------------------------------------------------------------------
