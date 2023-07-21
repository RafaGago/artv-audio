#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <stdint.h>
#include <type_traits>
#include <utility>

#include <boost/hana.hpp>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "artv-common/juce/look_and_feel.hpp"

#include "fonts/OpenSans-Bold.ttf.hpp"
#include "fonts/PassionsConflict-Regular.ttf.hpp"

#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/gui_util.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/parameters.hpp"
#include "artv-common/juce/value_tree_attachments.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "turbopaco/logo.hpp"
#include "turbopaco/parameters.hpp"

#define VERSION_INT VERSION_GET (VERSION_MAJOR, VERSION_MINOR, VERSION_REV)

namespace artv {

juce::Typeface::Ptr const default_typeface {
  juce::Typeface::createSystemTypefaceFor (
    font::OpenSans_Bold_ttf,
    sizeof font::OpenSans_Bold_ttf)};

juce::Typeface::Ptr const lcd_typeface = default_typeface;

juce::Typeface::Ptr const title_typeface {
  juce::Typeface::createSystemTypefaceFor (
    font::PassionsConflict_Regular_ttf,
    sizeof font::PassionsConflict_Regular_ttf)};
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
  bool upper_line = false;
  //----------------------------------------------------------------------------
  void paint (juce::Graphics& g) override
  {
    auto lb = getLocalBounds();
    auto bg = findColour (backgroundColourId);
    g.setColour (bg.darker (0.3f));
    g.fillRect (lb);

    auto inner = lb.reduced (1);
    g.setColour (bg);
    g.fillRect (inner);

    if (!upper_line) {
      return;
    }
    auto line = lb.removeFromTop (1);
    line      = line.reduced (1, 0);
    g.setColour (bg.brighter (0.2));
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
    // look and feel stuff
    _lf.setDefaultSansSerifTypeface (default_typeface);
    _lf.setDefaultLookAndFeel (&_lf);

    _lf.on_get_label_font = [this] (juce::Label& obj) {
      auto  f = obj.getFont();
      float h = obj.getTopLevelComponent()->getHeight();

      if (&obj == &_title) {
        f.setHeight (h * 0.14f);
      }
      else if (&obj == &_display_value) {
        f.setHeight (h * 0.12f);
      }
      else {
        f.setHeight (h * 0.055f);
        // there are many labels, so using "container_of" might be bug prone,
        // brute force instead
        slider_ext* sl = nullptr;
        _params.pforeach ([&, this] (auto type, auto& warray) {
          using T = std::decay_t<decltype (*warray[0])>;
          if constexpr (std::is_same_v<T, slider_ext>) {
            if (&obj == &warray[0]->label) {
              sl = &*warray[0];
            }
          }
        });
        if (sl) {
          auto& ed = sl->slider.edit;
          ed.setFont (f);
        }
      }
      return f;
    };
    _lf.on_rotary_draw = draw_rotary_1;

    _lf.on_get_combobox_font = [this] (juce::ComboBox& obj) {
      auto  f = _lf.getComboBoxFont (obj); // call LF V4
      float h = obj.getTopLevelComponent()->getHeight();
      f.setHeight (h * 0.2f);
      return f;
    };

    _lf.on_get_popup_menu_font = [this]() {
      auto f = _lf.getPopupMenuFont(); // call LF V4
      f.setHeight (this->getHeight() * 0.05f);
      return f;
    };

    // init
    _display_value.setFont (juce::Font {lcd_typeface});
    _title.setFont (juce::Font {title_typeface});
    _title.setText ("TurboPaco", juce::dontSendNotification);
    _logo = juce::Drawable::createFromImageData (logo_svg, sizeof logo_svg);

    register_mouse_events();
    // notice: the order in which they are added is the rendering order. This
    // is important to keep components on top without having to make complex
    // hierarchies, which suits simple UIs like this one.
    addAndMakeVisible (_header);
    addAndMakeVisible (_display_frame);
    addAndMakeVisible (_display_value);
    addAndMakeVisible (_title);
    addAndMakeVisible (*_logo);
    for (auto& v : _mainpanels) {
      addAndMakeVisible (v);
    }
    addAndMakeVisible (_footer);
    _params.init_widgets (*this, params); // knobs/sliders on top
    _display_value.setJustificationType (juce::Justification::centred);
    _title.setJustificationType (juce::Justification::centred);
    //_display_value.setFont (juce::Font {
    //  juce::Font::getDefaultMonospacedFontName(), 12, juce::Font::bold});

    auto opaque      = juce::Colour {0xff000000};
    auto transparent = juce::Colour {0};
    auto dark_grey   = juce::Colour {0xff000000 | 0x1A1C1E}.brighter (0.03f);
    auto light_grey  = juce::Colour {0xff000000 | 0x222526}.brighter (0.03f);
    auto display     = juce::Colour {0xff000000 | 0xB5BEAC};
    auto knob        = light_grey.brighter (0.2);
    auto knob_bg     = light_grey.darker (0.2);
    auto decay_color = juce::Colours::orange.brighter (0.07);
    auto track       = decay_color.brighter (0.8f);
    auto label_txt   = light_grey.brighter (1.8);
    auto outline     = dark_grey; // light_grey.brighter (0.2f);

    _logo->replaceColour (juce::Colours::black, label_txt);

    _lf.setColour (juce::PopupMenu::backgroundColourId, light_grey);
    _lf.setColour (juce::PopupMenu::highlightedBackgroundColourId, dark_grey);
    _lf.setColour (juce::Label::textColourId, label_txt);
    _lf.setColour (juce::Slider::trackColourId, track);
    _lf.setColour (juce::Slider::thumbColourId, knob);
    _lf.setColour (juce::Slider::backgroundColourId, knob_bg);

    _lf.setColour (juce::ComboBox::backgroundColourId, light_grey);
    _lf.setColour (juce::ComboBox::outlineColourId, outline);
    _lf.setColour (juce::ComboBox::arrowColourId, light_grey);
    _lf.setColour (juce::ComboBox::textColourId, label_txt);
    _lf.setColour (juce::PopupMenu::textColourId, label_txt);
    _lf.setColour (juce::PopupMenu::highlightedTextColourId, label_txt);
    _lf.setColour (juce::PopupMenu::backgroundColourId, light_grey);
    _lf.setColour (juce::PopupMenu::highlightedBackgroundColourId, dark_grey);

    _lf.setColour (juce::TextButton::buttonColourId, light_grey);
    _lf.setColour (juce::TextButton::textColourOffId, label_txt);

    _lf.setColour (juce::TextEditor::backgroundColourId, light_grey);
    _lf.setColour (juce::TextEditor::textColourId, label_txt);

    // different colors and styles
    set_color (
      juce::TextButton::buttonColourId, display.darker (0.1f), _display_frame);
    set_color (
      juce::Label::textColourId, display.darker (1.3f), _display_value);
    set_color (juce::Label::backgroundColourId, transparent, _display_value);
    set_color (juce::Label::textColourId, label_txt, _title);
    set_color (juce::Label::backgroundColourId, transparent, _title);
    set_color (panel::backgroundColourId, dark_grey, _header, _footer);
    for (auto& v : _mainpanels) {
      set_color (panel::backgroundColourId, light_grey, v);
    }

    using sliders = mp_list<parameters::dry, parameters::wet>;
    _params.pforeach (sliders {}, [=] (auto key, auto& warray) {
      warray[0]->slider.setSliderStyle (juce::Slider::LinearVertical);
    });

    using pans = mp_list<parameters::wet_pan, parameters::stereo>;
    _params.pforeach (pans {}, [=] (auto key, auto& warray) {
      warray[0]->slider.getProperties().set (slider_track_bg_from_center, true);
    });

    using combos
      = mp_list<parameters::algorithm, parameters::mode, parameters::clock>;
    _params.pforeach (combos {}, [=] (auto key, auto& warray) {
      warray[0]->combo.setJustificationType (juce::Justification::centred);
    });

    auto& decay = *_params.p_get (parameters::decay {})[0];
    decay.slider.setColour (juce::Slider::thumbColourId, decay_color);

    // changing text on chorus
    _params.p_get (parameters::algorithm {})[0]->on_change
      = [this] (combobox_ext& cb) {
          auto& decay   = *_params.p_get (parameters::decay {})[0];
          auto& mod     = *_params.p_get (parameters::mod {})[0];
          auto  nonotif = juce::NotificationType::dontSendNotification;
          bool  is_cho  = cb.combo.getText().containsIgnoreCase ("chorus");
          decay.label.setText (is_cho ? "Amt" : "Decay", nonotif);
          mod.label.setText (is_cho ? "Freq" : "Mod", nonotif);
        };

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

    mouse_event (
      event_common::source {}, mouse_event::type::exit, mouse_event::data {});
  }
  //----------------------------------------------------------------------------
  struct sizes {

    static constexpr float main_knob_wh_divs = 16.f;
    static constexpr float main_label_w_divs = main_knob_wh_divs;
    static constexpr float main_label_h_divs = main_knob_wh_divs / 3.f;

    static constexpr float w_separation_divs = 2.f;
    static constexpr float h_separation_divs = 2.f;

    static constexpr float header_h_divs = 0.75f * main_knob_wh_divs;
    static constexpr float footer_h_divs = 1.5f * main_label_h_divs;

    static constexpr float header_margin_w_divs  = 6.f;
    static constexpr float header_display_w_divs = 10.f;

    static constexpr float header_w_divs
      = header_margin_w_divs + header_display_w_divs + header_display_w_divs;

    static constexpr float footer_separation_w_divs = w_separation_divs;
    static constexpr float footer_section_w_divs    = 2 * main_knob_wh_divs;

    static constexpr float footer_w_divs = //
      footer_separation_w_divs + //
      footer_separation_w_divs + //
      footer_section_w_divs + //
      footer_separation_w_divs + //
      footer_separation_w_divs + //
      footer_section_w_divs + //
      footer_separation_w_divs + //
      footer_separation_w_divs + //
      footer_section_w_divs + //
      footer_separation_w_divs + //
      footer_separation_w_divs;

    static constexpr float main_w_separation_divs = 1.f;
    static constexpr float main_w_divs            = main_knob_wh_divs;
    static constexpr float main_h_divs
      = 2 * main_knob_wh_divs + 2 * main_label_h_divs;

    static constexpr float drywet_w_divs = main_w_divs * 1.5f;
    static constexpr float drywet_h_divs = main_h_divs;

    static constexpr float decay_area_w_divs  = main_w_divs * 1.75f;
    static constexpr float decay_area_h_divs  = main_h_divs;
    static constexpr float decay_label_h_divs = main_label_h_divs;
    static constexpr float decay_knob_h_divs
      = decay_area_h_divs - decay_label_h_divs;

    static constexpr float total_w_divs = //
      w_separation_divs + //
      w_separation_divs + //
      drywet_w_divs + //
      w_separation_divs + //
      w_separation_divs + //
      main_w_divs + // predelay + oprange
      w_separation_divs + //
      w_separation_divs + //
      decay_area_w_divs + //
      w_separation_divs + //
      w_separation_divs + //
      main_w_divs + // charac + mod
      w_separation_divs + //
      w_separation_divs + //
      main_w_divs + // freq
      w_separation_divs + //
      w_separation_divs + //
      main_w_divs + // stereo
      w_separation_divs + //
      w_separation_divs + //
      main_w_divs + // ducking
      w_separation_divs + //
      w_separation_divs //
      ;
    static constexpr float total_h_divs = //
      h_separation_divs + //
      header_h_divs + //
      h_separation_divs + //
      h_separation_divs + //
      main_h_divs + //
      h_separation_divs + //
      h_separation_divs + //
      footer_h_divs + //
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

    auto title = header.removeFromLeft (header_margin_w);
    _title.setBounds (title);

    //    auto tf = _title.getFont();
    //    tf.setHeight ((float) title.getHeight() * 113.25f);
    //    _title.setFont (tf);
    auto logo = header.removeFromRight (header_margin_w).toFloat();
    logo.removeFromRight (logo.getWidth() * 0.05f);

    _logo->setTransformToFit (
      logo.reduced (logo.getHeight() * 0.13f),
      juce::RectanglePlacement {juce::RectanglePlacement::Flags::centred});

    auto display = header;
    _display_frame.setBounds (display);
    auto pvfbounds = _display_frame.getBounds().reduced (h_h / 10.f);

    pvfbounds.reduce (0, pvfbounds.getHeight() / 9);
    _display_value.setBounds (pvfbounds);

    //    auto f = _display_value.getFont();
    //    f.setHeight ((float) pvfbounds.getHeight() * 3.25f);
    //    _display_value.setFont (f);

    auto main
      = area.removeFromTop (sep_h + (h * (float) sizes::main_h_divs) + sep_h);

    auto drywet = main.removeFromLeft (
      2.f * sep_w + (w * (float) sizes::drywet_w_divs) + sep_w);
    _mainpanels[0].setBounds (drywet);
    drywet.removeFromTop (sep_h);
    drywet.removeFromBottom (sep_h);
    drywet.removeFromLeft (2.f * sep_w);
    drywet.removeFromRight (sep_w);
    grid (
      drywet,
      ((float) drywet.getWidth()) / 2.f,
      make_array (label_h, label_h + 2 * knob_h),
      grid_label_on_top (*_params.p_get (parameters::dry {})[0]),
      grid_label_on_top (*_params.p_get (parameters::wet {})[0]));

    auto predel
      = main.removeFromLeft (sep_w + (w * (float) sizes::main_w_divs) + sep_w);
    _mainpanels[1].setBounds (predel);
    predel.removeFromTop (sep_h);
    predel.removeFromBottom (sep_h);
    predel.removeFromLeft (sep_w);
    predel.removeFromRight (sep_w);
    grid (
      predel,
      predel.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::predelay {})[0]),
      grid_label_on_top (*_params.p_get (parameters::operating_range {})[0]));

    auto decay = main.removeFromLeft (
      sep_w + (w * (float) sizes::decay_area_w_divs) + sep_w);
    _mainpanels[2].setBounds (decay);
    decay.removeFromTop (sep_h);
    decay.removeFromBottom (sep_h);
    decay.removeFromLeft (sep_w);
    decay.removeFromRight (sep_w);
    auto decay_knob_h  = sizes::decay_knob_h_divs * h;
    auto decay_label_h = sizes::decay_label_h_divs * h;
    grid (
      decay,
      decay.getWidth(),
      make_array (decay_label_h, decay_knob_h),
      grid_label_on_top (*_params.p_get (parameters::decay {})[0]));

    auto mod
      = main.removeFromLeft (sep_w + (w * (float) sizes::main_w_divs) + sep_w);
    _mainpanels[3].setBounds (mod);
    mod.removeFromTop (sep_h);
    mod.removeFromBottom (sep_h);
    mod.removeFromLeft (sep_w);
    mod.removeFromRight (sep_w);
    grid (
      mod,
      mod.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::mod {})[0]),
      grid_label_on_top (*_params.p_get (parameters::character {})[0]));

    auto freq
      = main.removeFromLeft (sep_w + (w * (float) sizes::main_w_divs) + sep_w);
    _mainpanels[4].setBounds (freq);
    freq.removeFromTop (sep_h);
    freq.removeFromBottom (sep_h);
    freq.removeFromLeft (sep_w);
    freq.removeFromRight (sep_w);
    grid (
      freq,
      freq.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::hf_amt {})[0]),
      grid_label_on_top (*_params.p_get (parameters::lf_amt {})[0]));

    auto stereo
      = main.removeFromLeft (sep_w + (w * (float) sizes::main_w_divs) + sep_w);
    _mainpanels[5].setBounds (stereo);
    stereo.removeFromTop (sep_h);
    stereo.removeFromBottom (sep_h);
    stereo.removeFromLeft (sep_w);
    stereo.removeFromRight (sep_w);
    grid (
      stereo,
      stereo.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::stereo {})[0]),
      grid_label_on_top (*_params.p_get (parameters::wet_pan {})[0]));

    auto duck = main.removeFromLeft (main.getWidth());
    _mainpanels[6].setBounds (duck);
    duck.removeFromTop (sep_h);
    duck.removeFromBottom (sep_h);
    duck.removeFromLeft (sep_w);
    duck.removeFromRight (duck.getWidth() - (w * (float) sizes::main_w_divs));
    grid (
      duck,
      duck.getWidth(),
      make_array (label_h, knob_h, label_h, knob_h),
      grid_label_on_top (*_params.p_get (parameters::dyn_threshold {})[0]),
      grid_label_on_top (*_params.p_get (parameters::dyn_speed {})[0]));

    auto footer
      = area.removeFromTop (sep_h + (sizes::footer_h_divs * h) + sep_h);
    _footer.setBounds (footer);
    footer.removeFromTop (sep_h);
    footer.removeFromBottom (sep_h);

    auto f_w    = (float) w_units / sizes::footer_w_divs;
    auto f_sep  = f_w * sizes::footer_separation_w_divs;
    auto f_sect = f_w * sizes::footer_section_w_divs;

    auto s1 = footer.removeFromLeft (2 * f_sep + f_sect + f_sep);
    s1.removeFromLeft (2 * f_sep);
    auto clock = s1.removeFromLeft (f_sect);
    grid (
      clock,
      (float) clock.getWidth(),
      make_array ((float) clock.getHeight()),
      grid_buttons_on_same_row (
        *_params.p_get (parameters::clock {})[0], 0.8f));

    auto s2 = footer.removeFromLeft (f_sep + f_sect + f_sep);
    s2.removeFromLeft (f_sep);
    auto algo = s2.removeFromLeft (f_sect);
    grid (
      algo,
      (float) algo.getWidth(),
      make_array ((float) algo.getHeight()),
      grid_buttons_on_same_row (
        *_params.p_get (parameters::algorithm {})[0], 0.8f));

    auto s3 = footer.removeFromLeft (f_sep + f_sect + f_sep * 2);
    s3.removeFromLeft (f_sep);
    auto mode = s3.removeFromLeft (f_sect);
    grid (
      mode,
      (float) mode.getWidth(),
      make_array ((float) mode.getHeight()),
      grid_buttons_on_same_row (*_params.p_get (parameters::mode {})[0], 0.8f));
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
    event_common::source src,
    mouse_event::type    event_type,
    mouse_event::data    data)
  {
    juce::String     text;
    juce::Component* component {};

    auto slider_2ptr   = std::get_if<juce::Slider*> (&src);
    auto button_2ptr   = std::get_if<juce::Button*> (&src);
    auto combobox_2ptr = std::get_if<juce::ComboBox*> (&src);
    auto label_2ptr    = std::get_if<juce::Label*> (&src);

    if (slider_2ptr) {
      auto& slider       = **slider_2ptr;
      component          = &slider;
      juce::String value = slider.getTextFromValue (slider.getValue());
      text               = slider.getName();
      auto suffix        = slider.getTextValueSuffix();
      if (suffix.length() > 0) {
        text += " (" + suffix + ")";
      }
      text += ": " + value.removeCharacters (suffix);
      text = text.trim();
    }
    else if (button_2ptr) {
      auto& button = **button_2ptr;
      component    = &button;
      text         = button.getName();
      if (!text.endsWith (" Next") && !text.endsWith (" Prev")) {
        text += ": ";
        text += button.getToggleState() ? "On" : "Off";
      }
    }
    else if (combobox_2ptr) {
      auto& combobox = **combobox_2ptr;
      component      = &combobox;
      text           = combobox.getName();
      text += ": ";
      text += combobox.getText();
    }
    else if (label_2ptr) {
      auto& label = **label_2ptr;
      if (&label == &_title) {
        text = "| Turbo Paco v" VERSION_TXT " |";
      }
    }
    else {}
    if (event_type == mouse_event::exit) {
      text = ":)";
    }
    if (component && !component->isEnabled()) {
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
            = std::bind (&editor::mouse_event, this, _1, _2, _4);
        }
        if constexpr (is_button_ext<type>::value) {
          w->button.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _4);
        }
        if constexpr (is_toggle_buttons<type>::value) {
          for (auto& btn : w->buttons) {
            btn.on_mouse_event
              = std::bind (&editor::mouse_event, this, _1, _2, _4);
          }
        }
        if constexpr (is_combobox_ext<type>::value) {
          w->combo.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _4);
          w->prev.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _4);
          w->next.on_mouse_event
            = std::bind (&editor::mouse_event, this, _1, _2, _4);
        }
      }
    });
    _title.on_mouse_event = std::bind (&editor::mouse_event, this, _1, _2, _4);
  }
  //----------------------------------------------------------------------------
private:
  juce::AudioProcessor& _processor; // unused
  saner_look_and_feel   _lf;

  editor_apvts_widgets<parameters::parameters_typelist> _params;

  panel                                 _header;
  std::array<panel, 7>                  _mainpanels;
  panel                                 _footer;
  juce::Label                           _display_value;
  add_juce_mouse_callbacks<juce::Label> _title;
  std::unique_ptr<juce::Drawable>       _logo;
  juce::TextButton                      _display_frame;
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
