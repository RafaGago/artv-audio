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

#include <ff_meters.h>

#include "artv-common/juce/effect_base.hpp"
#include "artv-common/juce/gui_util.hpp"
#include "artv-common/juce/look_and_feel.hpp"
#include "artv-common/juce/math.hpp"
#include "artv-common/juce/parameters.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

#include "mix-maxtrix/parameters.hpp"

#define VERSION_INT VERSION_GET (VERSION_MAJOR, VERSION_MINOR, VERSION_REV)

namespace artv {

// -----------------------------------------------------------------------------
struct VerticalLine : public juce::Component {
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
struct radio_id {
  enum {
    solo_mute_beg  = (0 * parameters::n_stereo_busses) + 1,
    stereo_cfg_beg = (1 * parameters::n_stereo_busses) + 1,
    page_beg       = (2 * parameters::n_stereo_busses) + 1,
  };
};
// -----------------------------------------------------------------------------
class editor : public juce::AudioProcessorEditor,
               public juce::DragAndDropContainer,
               private has_editor_params<parameters::parameters_typelist> {
public:
  //----------------------------------------------------------------------------
  explicit editor (
    juce::AudioProcessor&               p,
    juce::AudioProcessorValueTreeState& params,
    crange<foleys::LevelMeterSource>    meter_srcs)
    : AudioProcessorEditor (p), _processor (p)
  {
    //        id, juce::Colours::green.brighter (0.2),
    //        juce::Colours::orange.brighter (0.07),
    //        juce::Colours::teal,
    //        juce::Colours::darkgoldenrod,
    //        juce::Colours::red.brighter (0.4),
    //        juce::Colours::sienna,

    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      _meters[i].setLookAndFeel (&_meters_look_feel[i]);
      _meters[i].setMeterSource (&meter_srcs[i]);
      meter_srcs[i].setSuspended (false);
      _meters[i].setMeterFlags (foleys::LevelMeter::MeterFlags::Minimal);
      addAndMakeVisible (_meters[i]);
    }
    _meter_srcs = meter_srcs;

    setLookAndFeel (&_lookfeel);

    // setting up the dummy lines and group elements
    for (auto& ln : _side_lines) {
      ln.setLookAndFeel (&_lookfeel);
      addAndMakeVisible (ln);
    }
    _side_lines[0].left = false;
    _side_lines[1].left = true;

    // lines for the fx modifiers
    for (uint i = 0; i < _fx_lines.size(); ++i) {
      _fx_lines[i].left = !(i % 2);
      _fx_lines[i].setLookAndFeel (&_lookfeel);
      addAndMakeVisible (_fx_lines[i]);
    }

    init_widgets (*this, params);

    register_mouse_events();

    // set slider styles
    p_get (parameters::global_volume {})[0]->slider.setSliderStyle (
      juce::Slider::LinearHorizontal);

    for (auto& s : p_get (parameters::volume {})) {
      s->slider.setSliderStyle (juce::Slider::LinearVertical);
    }

    // setting descriptive names to the IO buttons
    auto& ins  = p_get (parameters::in_selection {});
    auto& outs = p_get (parameters::out_selection {});
    for (uint src = 0; src < ins.size(); ++src) {
      for (uint dst = 0; dst < ins.size(); ++dst) {
        // TODO: get {fmt}
        char str[64];
        str[0] = 0;
        snprintf (str, sizeof str, "In%u to Mix%u", dst + 1, src + 1);
        ins[src]->buttons[dst].setName (str);

        str[0] = 0;
        snprintf (str, sizeof str, "Mix%u to Out%u", src + 1, dst + 1);
        outs[src]->buttons[dst].setName (str);
      }
    }

    // setting descriptive names to the channel modifiers
    auto& mods = p_get (parameters::channel_modifs {});
    for (auto& widg_ptr : mods) {
      widg_ptr->buttons[0].setName ("Phase inv");
      widg_ptr->buttons[1].setName ("Mono L");
      widg_ptr->buttons[2].setName ("Mono R");
      widg_ptr->buttons[3].setName ("L/R swap");
    }

    // setting descriptive names to the mixer to mixer send buttons
    for (uint i = 1; i < parameters::n_stereo_busses; ++i) {
      // TODO: get {fmt}
      char str[64];
      str[0] = 0;
      snprintf (str, sizeof str, "Mix%u to Mix%u", i, i + 1);
      p_get (parameters::mixer_sends {})[0]->buttons[i - 1].setName (str);
    }

    // setting a better description to the Dry/Wet parameters modifiers
    for (auto& s : p_get (parameters::dry_balance {})) {
      s->slider.setName ("Dry M/S");
    }
    for (auto& s : p_get (parameters::wet_balance {})) {
      s->slider.setName ("Wet M/S");
    }
    for (auto& s : p_get (parameters::dry_pan {})) {
      s->slider.setName ("Dry Pan");
    }
    for (auto& s : p_get (parameters::wet_pan {})) {
      s->slider.setName ("Wet Pan");
    }
    for (auto& s : p_get (parameters::fx_mix {})) {
      s->slider.setName ("Fx Mix");
    }
    for (auto& s : p_get (parameters::pan {})) {
      s->slider.setName ("Main Pan");
    }

    // setup disabled sliders not bound to parameters.
    for (auto& arr : _fx_off_sliders) {
      for (slider_ext& s : arr) {
        s.init ("", "", *this);
        s.slider.setRange (0.0, 1., 1.);
        s.slider.setName ("Disabled");
        s.slider.setEnabled (false);
        using namespace std::placeholders;
        s.slider.on_mouse_event
          = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
      }
    }

    // setup page buttons.
    for (uint chnl = 0; chnl < _page_buttons.size(); ++chnl) {
      uint pb_size = _page_buttons[0].size();
      for (uint b = 0; b < pb_size; ++b) {
        // TODO: get {fmt}
        auto& btn = _page_buttons[chnl][b];
        char  str[64];
        str[0] = 0;
        snprintf (str, sizeof str, "Page %u", b + 1);
        btn.setButtonText (str);
        snprintf (str, sizeof str, "Show fx page %u", b + 1);
        btn.setName (str);
        btn.setClickingTogglesState (true);
        btn.setToggleState (
          b == 0, juce::NotificationType::dontSendNotification);
        btn.onStateChange = [=]() { on_fx_type_or_page_change (chnl, false); };
        uint right = b < (pb_size - 1) ? juce::TextButton::ConnectedOnRight : 0;
        uint left  = (b != 0) ? juce::TextButton::ConnectedOnLeft : 0;
        btn.setConnectedEdges (left | right);
        using namespace std::placeholders;
        btn.setLookAndFeel (&_lookfeel);
        btn.on_mouse_event
          = std::bind (&editor::mouse_event, this, _1, _2, _3, _4);
        addAndMakeVisible (btn);
      }
    }

    // initialize header elements label
    _about.setLookAndFeel (&_lookfeel);
    _parameter_value_frame.setLookAndFeel (&_lookfeel);
    _parameter_value.setLookAndFeel (&_lookfeel);

    _about.setButtonText ("Help/Credits");
    _about.onClick = [=] { this->show_about_popup(); };

    addAndMakeVisible (_about);
    addAndMakeVisible (_parameter_value_frame);
    addAndMakeVisible (_parameter_value);

    _parameter_value.setJustificationType (juce::Justification::left);
    //_parameter_value.setFont (juce::Font {
    //  juce::Font::getDefaultMonospacedFontName(), 12, juce::Font::bold});

    auto display_color = juce::Colours::teal;

    set_color (
      juce::TextButton::buttonColourId,
      display_color.brighter (1.),
      _parameter_value_frame);

    set_color (
      juce::Label::textColourId, display_color.darker (0.7), _parameter_value);

    set_color (
      juce::Label::backgroundColourId, juce::Colour (0), _parameter_value);

#if 0
    set_color (
      juce::Label::outlineColourId,
      juce::Colours::teal.darker (2.),
      _parameter_name,
      _parameter_value_frame);
#endif
    // colors
    auto ch_color = make_array<juce::Colour> (
      juce::Colours::teal,
      juce::Colours::orange.brighter (0.07),
      juce::Colours::red.brighter (0.4),
      juce::Colour (0xff00cc99),
      juce::Colour (0xff0099cc), // 5
      juce::Colour (0xffff794d),
      juce::Colour (0xffd22d6f),
      juce::Colour (0xff9cb946));

    pforeach (
      parameters::all_channel_sliders_typelist {},
      [=] (auto type, auto& warray) {
        for (uint i = 0; i < warray.size(); ++i) {
          warray[i]->slider.setColour (
            juce::Slider::thumbColourId, ch_color[i]);
        }
      });

    // make the FX off sliders lighter
    for (auto& arr : _fx_off_sliders) {
      for (slider_ext& s : arr) {
        s.slider.setColour (
          juce::Slider::thumbColourId,
          // s.slider.findColour (juce::Slider::thumbColourId).brighter
          // (0.3)); juce::Colours::grey);
          //_lookfeel.findColour (juce::ResizableWindow::backgroundColourId)
          _lookfeel.findColour (juce::Slider::backgroundColourId));
      }
    }

    // meter colors
    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmMeterGradientLowColour,
        ch_color[i].brighter (1.1));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmMeterGradientMidColour,
        ch_color[i].brighter (0.8));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmMeterGradientMaxColour,
        ch_color[i].brighter (0.8));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmMeterBackgroundColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmOutlineColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmBackgroundClipColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmTicksColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmBackgroundColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmBackgroundColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
      _meters_look_feel[i].setColour (
        foleys::LevelMeter::lmMeterOutlineColour,
        _lookfeel.findColour (juce::ResizableWindow::backgroundColourId));
    }

    auto& mutesolo = p_get (parameters::mute_solo {});

    // button colors
    constexpr uint nbuses = parameters::n_stereo_busses;
    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      set_color (
        juce::TextButton::buttonOnColourId,
        ch_color[i].brighter (0.5),
        ins[0]->buttons[i],
        ins[1]->buttons[i],
        ins[2]->buttons[i],
        ins[3]->buttons[i],
        outs[0]->buttons[i],
        outs[1]->buttons[i],
        outs[2]->buttons[i],
        outs[3]->buttons[i],
#if MIXMAXTRIX_CHANNELS == 8
        ins[4]->buttons[i],
        ins[5]->buttons[i],
        ins[6]->buttons[i],
        ins[7]->buttons[i],
        outs[4]->buttons[i],
        outs[5]->buttons[i],
        outs[6]->buttons[i],
        outs[7]->buttons[i],
#endif
        mods[i]->buttons,
        mutesolo[i]->buttons,
        make_contiguous_range (_page_buttons[i]));
    }
    set_color (
      juce::TextButton::buttonOnColourId,
      _lookfeel.findColour (juce::Slider::backgroundColourId).brighter (0.5),
      p_get (parameters::mixer_sends {})[0]->buttons);

    // radio grouping
    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      mutesolo[i]->set_radio_group (radio_id::solo_mute_beg + i, true);
      mods[i]->set_radio_group (radio_id::stereo_cfg_beg + i, 1, 3, true);
      for (auto& btn : _page_buttons[i]) {
        btn.setRadioGroupId (
          radio_id::page_beg + i, juce::NotificationType::dontSendNotification);
      }
    }

    // pan slider bg's.
    mp11::mp_for_each<mp_list<parameters::dry_balance, parameters::pan>> (
      [this] (auto w) {
        for (auto& s : p_get (w)) {
          s->slider.getProperties().set (
            look_and_feel::slider_track_bg_from_center, true);
        }
      });

    // fx type events
    for (uint i = 0; i < parameters::n_stereo_busses; ++i) {
      auto& combo = p_get (parameters::fx_type {})[i]->combo;
      // fx change
      combo.onChange = [=] { on_fx_type_or_page_change (i, true); };
      on_fx_type_or_page_change (i, true);
      // drag drop
      using namespace std::placeholders;
      combo.is_drag_source = true;
      combo.on_drop
        = std::bind (&editor::fx_drag_and_drop_event, this, i, _1, _2, _3);
    }

    // size
    constexpr float ratio  = sizes::w_divs / sizes::h_divs;
    constexpr float factor = 0.9f;

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
      height = width / factor;
    }
    setSize (width, height);
  }
  //----------------------------------------------------------------------------
  ~editor() override
  {
    // this lookanfeel thing could really be a bit smarter...
    clear_widgets();
    _parameter_value.setLookAndFeel (nullptr);
    _parameter_value_frame.setLookAndFeel (nullptr);
    for (auto& a1 : _fx_off_sliders) {
      for (auto& s : a1) {
        s.clear();
      }
    }
    for (auto& btnarr : _page_buttons) {
      for (auto& btn : btnarr) {
        btn.setLookAndFeel (nullptr);
      }
    }
    for (auto& ln : _side_lines) {
      ln.setLookAndFeel (nullptr);
    }
    for (auto& m : _meters) {
      m.setLookAndFeel (nullptr);
    }
    for (auto& src : _meter_srcs) {
      src.setSuspended (true);
    }
    setLookAndFeel (nullptr);
  }
  //----------------------------------------------------------------------------
  struct sizes {
    static constexpr float columns = parameters::n_stereo_busses;

    static constexpr float sqr_btn_divs         = 4.f;
    static constexpr float fx_combobox_row_divs = sqr_btn_divs * 0.66;
    static constexpr float fx_combobox_prevnext_divs
      = fx_combobox_row_divs * 0.95;

    static constexpr float col_separation_divs = 2.5f;
    static constexpr float row_separation_divs = 1.f;

    static constexpr float header_row_divs = 1 * sqr_btn_divs;

    static constexpr float header_sep_row_divs      = 1.;
    static constexpr float about_row_divs           = 9.;
    static constexpr float parameter_value_row_divs = 16.;
    static constexpr float routing_row_divs         = 4.;
    static constexpr float global_volume_row_divs   = 4.;

    // clang-format off
    static constexpr float header_w_divs =
      about_row_divs +
      header_sep_row_divs +
      parameter_value_row_divs +
      header_sep_row_divs +
      routing_row_divs +
      header_sep_row_divs +
      global_volume_row_divs;
    // clang-format on

    static constexpr float btns_per_column = 4.;

    // referencing to the width of the 4 buttons
    static constexpr float col_divs = sqr_btn_divs * btns_per_column;

    static constexpr float dry_wet_knobs_slider_row_divs = col_divs / 3.;
    static constexpr float dry_wet_knobs_label_row_divs  = sqr_btn_divs / 2.;
    static constexpr float dry_wet_knobs_row_divs
      = dry_wet_knobs_slider_row_divs + dry_wet_knobs_label_row_divs;

    static constexpr float fx_param_page_row_divs = sqr_btn_divs;

    static constexpr float fx_param_row_divs       = col_divs / 2.;
    static constexpr float fx_param_label_row_divs = sqr_btn_divs / 2.;
    // TODO: this 4 here is half the fx params
    static constexpr float fx_params_row_divs
      = 4 * (fx_param_row_divs + fx_param_label_row_divs);

    static constexpr float fader_row_divs       = 1.5 * col_divs;
    static constexpr float fader_label_row_divs = fx_param_label_row_divs;
    static constexpr float fader_sect_row_divs
      = fader_row_divs + fader_label_row_divs;

    static constexpr float pan_width_slider_row_divs = sqr_btn_divs;
    static constexpr float pan_width_label_row_divs  = fader_label_row_divs;
    static constexpr float pan_width_row_divs
      = (pan_width_slider_row_divs + pan_width_label_row_divs) * 2.;

    static constexpr float mute_solo_button_row_divs = sqr_btn_divs;
    static constexpr float mute_solo_button_col_divs = sqr_btn_divs * 2.;

    static constexpr float io_btns_row_divs
      = sqr_btn_divs * (parameters::n_stereo_busses / 4);

    static constexpr float w_divs
      = (columns * col_divs) + ((columns + 1.f) * col_separation_divs);

    static constexpr float h_divs = row_separation_divs + // upper margin
      header_row_divs + // upper info bar
      row_separation_divs + // margin
      io_btns_row_divs + // input selection buttons
      row_separation_divs + // margin
      fx_combobox_row_divs + // Fx
      fx_combobox_prevnext_divs + // Fx prev/next
      dry_wet_knobs_row_divs + // Fx mod rotaries
      row_separation_divs + // margin
      fx_param_page_row_divs + //
      row_separation_divs + // margin
      fx_params_row_divs + // params
      row_separation_divs + // margin
      fader_sect_row_divs + // fader
      row_separation_divs + // margin
      dry_wet_knobs_row_divs + // Fx mod rotaries
      row_separation_divs + // margin
      mute_solo_button_row_divs + // mute/solo
      row_separation_divs + // margin
      sqr_btn_divs + // stereo section
      row_separation_divs + // margin
      io_btns_row_divs + // output selection buttons
      row_separation_divs; // lower margin
  };
  //----------------------------------------------------------------------------
  void resized() override
  {
    auto area = getLocalBounds();

    auto h_units = (float) area.getHeight();
    auto w_units = (float) area.getWidth();

    auto h = (float) h_units / (float) sizes::h_divs;
    auto w = (float) w_units / (float) sizes::w_divs;

    auto sqr_btn_w = sizes::sqr_btn_divs * w;
    auto sqr_btn_h = sizes::sqr_btn_divs * h;

    auto sep_w = sizes::col_separation_divs * w;
    auto sep_h = sizes::row_separation_divs * h;

    auto col_w = sizes::col_divs * w;

    auto get_columns = [=] (auto& area) {
      std::array<juce::Rectangle<int>, parameters::n_stereo_busses> areas;
      area.removeFromLeft (sep_w); // left margin
      for (size_t i = 0; i < areas.size(); ++i) {
        areas[i] = area.removeFromLeft (col_w);
        area.removeFromLeft (sep_w);
      }
      return areas;
    };

    area.removeFromTop (sep_h); //  upper margin

    // upper bar
    float header_h = sizes::header_row_divs * h;
    auto  h_w = ((float) w_units - 2 * sep_w) / (float) sizes::header_w_divs;

    // auto  header_w = hw sizes::header_column_divs;
    auto header = area.removeFromTop (header_h);
    header.removeFromLeft (sep_w);
    header.removeFromRight (sep_w);

    _about.setBounds (header.removeFromLeft (h_w * sizes::about_row_divs));

    header.removeFromLeft (h_w * sizes::header_sep_row_divs);

    _parameter_value_frame.setBounds (
      header.removeFromLeft (h_w * sizes::parameter_value_row_divs));

    header.removeFromLeft (h_w * sizes::header_sep_row_divs);

    auto pvfbounds = _parameter_value_frame.getBounds().reduced (header_h / 10);
    // make the font smaller to avoid resizings
    pvfbounds.reduce (0, pvfbounds.getHeight() / 9);
    _parameter_value.setBounds (pvfbounds);

    auto routing = header.removeFromLeft (h_w * sizes::routing_row_divs);

    grid (
      routing,
      (float) routing.getWidth(),
      make_array (header_h / 2.f, header_h / 2.f),
      *p_get (parameters::routing {})[0]);

    header.removeFromLeft (h_w * sizes::header_sep_row_divs);

    auto main_gain
      = header.removeFromLeft (h_w * sizes::global_volume_row_divs);

    grid (
      main_gain,
      (float) main_gain.getWidth(),
      make_array (header_h / 2.f, header_h / 2.f),
      *p_get (parameters::global_volume {})[0]);

    area.removeFromTop (sep_h); // separator

    // input routings
    auto io_h        = sizes::io_btns_row_divs * h;
    auto in_routings = area.removeFromTop (io_h);
    auto columns     = get_columns (in_routings);

    for (uint i = 0; i < columns.size(); ++i) {
      auto& in_btns = p_get (parameters::in_selection {})[i]->buttons;
      grid (
        columns[i],
        sqr_btn_w,
#if MIXMAXTRIX_CHANNELS == 4
        make_array (sqr_btn_h),
        in_btns);
#elif MIXMAXTRIX_CHANNELS == 8
        make_array (sqr_btn_h, sqr_btn_h),
        in_btns[0],
        in_btns[4],
        in_btns[1],
        in_btns[5],
        in_btns[2],
        in_btns[6],
        in_btns[3],
        in_btns[7]);
#endif
    }

    area.removeFromTop (sep_h); // separator

    // fx
    auto fx_type_h      = sizes::fx_combobox_row_divs * h;
    auto fx_prev_next_h = sizes::fx_combobox_prevnext_divs * h;
    auto fx             = area.removeFromTop (fx_type_h + fx_prev_next_h);
    columns             = get_columns (fx);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        (float) columns[0].getWidth(),
        make_array (fx_type_h, fx_prev_next_h),
        *p_get (parameters::fx_type {})[i]);
    }

    area.removeFromTop (sep_h); // separator

    // fx modifiers
    auto dry_wet_mods_slider_h = sizes::dry_wet_knobs_slider_row_divs * h;
    auto dry_wet_mods_label_h  = sizes::dry_wet_knobs_label_row_divs * h;
    auto fx_dry_wet_mods
      = area.removeFromTop (dry_wet_mods_slider_h + dry_wet_mods_label_h);
    columns = get_columns (fx_dry_wet_mods);

    for (uint i = 0; i < columns.size(); ++i) {
      _fx_lines[i * 2].setBounds (columns[i]);
      _fx_lines[(i * 2) + 1].setBounds (columns[i]);
      grid (
        columns[i],
        (float) columns[0].getWidth() / 3.,
        make_array (dry_wet_mods_slider_h, dry_wet_mods_label_h),
        *p_get (parameters::wet_balance {})[i],
        *p_get (parameters::wet_pan {})[i],
        *p_get (parameters::fx_mix {})[i]);
    }

    area.removeFromTop (sep_h); // separator

    // fx page
    auto fx_page_h = sizes::fx_param_page_row_divs * h;
    auto fx_page   = area.removeFromTop (fx_page_h);
    columns        = get_columns (fx_page);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        (float) columns[i].getWidth() / 2.,
        make_array (fx_page_h),
        make_contiguous_range (_page_buttons[i]));
    }

    area.removeFromTop (sep_h); // separator

    auto mix_snd_btn_rm_header
      = getLocalBounds().getHeight() - area.getHeight();

    // fx params
    auto fxpars_h         = sizes::fx_params_row_divs * h;
    auto fx_param_h       = sizes::fx_param_row_divs * h;
    auto fx_param_label_h = sizes::fx_param_label_row_divs * h;
    auto fxpars           = area.removeFromTop (fxpars_h);
    columns               = get_columns (fxpars);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        (float) columns[i].getWidth() / 2,
        make_array (
          fx_param_h,
          fx_param_label_h,
          fx_param_h,
          fx_param_label_h,
          fx_param_h,
          fx_param_label_h,
          fx_param_h,
          fx_param_label_h),
        _fx_off_sliders[i][0],
        _fx_off_sliders[i][2],
        _fx_off_sliders[i][4],
        _fx_off_sliders[i][6],
        _fx_off_sliders[i][1],
        _fx_off_sliders[i][3],
        _fx_off_sliders[i][5],
        _fx_off_sliders[i][7]);
    }
    // setting all the other dynamically shown/hiddnen FX sliders onto the
    // same positions than the fx-off sliders.
    auto cp_slider_ex_bounds
      = [&] (uint row, uint col, auto t) { using dst_type = decltype (t); };

    // positioning  all fx parameter sliders on top of each respective
    // _fx_off_slider object
    mp11::mp_for_each<parameters::all_fx_typelists> ([&] (auto fx_tlist) {
      uint param_idx = 0;
      mp11::mp_for_each<decltype (fx_tlist)> ([&] (auto param) {
        auto& param_arr = p_get (param);
        for (uint chnl = 0; chnl < param_arr.size(); ++chnl) {

          param_arr[chnl]->slider.setBounds (
            _fx_off_sliders[chnl][param_idx % num_fx_page_params]
              .slider.getBounds());

          param_arr[chnl]->label.setBounds (
            _fx_off_sliders[chnl][param_idx % num_fx_page_params]
              .label.getBounds());
        }
        ++param_idx;
      });
    });

    // Fader section
    area.removeFromTop (sep_h); // separator

    auto fader_h            = sizes::fader_row_divs * h;
    auto fader_label_h      = sizes::fader_label_row_divs * h;
    auto fader_label_sect_h = sizes::fader_sect_row_divs * h;
    auto fader              = area.removeFromTop (fader_label_sect_h);
    columns                 = get_columns (fader);

    for (uint i = 0; i < columns.size(); ++i) {
      auto fader_area = columns[i];
      auto meter_area = columns[i];

      meter_area.removeFromBottom (fader_label_h);
      meter_area.removeFromLeft ((float) fader_area.getWidth() / 2.8);
      meter_area.removeFromRight ((float) fader_area.getWidth() / 2.8);
      meter_area = meter_area.reduced (0, (float) meter_area.getHeight() / 25.);
      _meters[i].setBounds (meter_area);

      grid (
        fader_area,
        (float) fader_area.getWidth(),
        make_array (fader_h, fader_label_h),
        *p_get (parameters::volume {})[i]);
    }

    area.removeFromTop (sep_h); // separator
    // pan + width
    auto dry_mods
      = area.removeFromTop (dry_wet_mods_slider_h + dry_wet_mods_label_h);
    columns = get_columns (dry_mods);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        (float) columns[0].getWidth() / 3.,
        make_array (dry_wet_mods_slider_h, dry_wet_mods_label_h),
        *p_get (parameters::dry_balance {})[i],
        *p_get (parameters::dry_pan {})[i],
        *p_get (parameters::pan {})[i]);
    }

    auto mix_snd_btn_rm_footer = area.getHeight();

    area.removeFromTop (sep_h); // separator

    // Mute / solo
    auto mute_solo_button_w = sizes::mute_solo_button_col_divs * w;
    auto mute_solo_h        = sizes::mute_solo_button_row_divs * h;
    auto mute_solo_area     = area.removeFromTop (mute_solo_h);
    columns                 = get_columns (mute_solo_area);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        mute_solo_button_w,
        make_array (mute_solo_h),
        p_get (parameters::mute_solo {})[i]->buttons);
    }

    area.removeFromTop (sep_h); // separator

    // Channel modifiers
    auto stereo_buttons = area.removeFromTop (sqr_btn_h);
    columns             = get_columns (stereo_buttons);

    for (uint i = 0; i < columns.size(); ++i) {
      grid (
        columns[i],
        sqr_btn_w,
        make_array (sqr_btn_h),
        p_get (parameters::channel_modifs {})[i]->buttons);
    }

    area.removeFromTop (sep_h); // separator

    // Output destination selection section
    auto outs = area.removeFromTop (io_h);
    columns   = get_columns (outs);

    for (uint i = 0; i < columns.size(); ++i) {
      auto& out_btns = p_get (parameters::out_selection {})[i]->buttons;
      grid (
        columns[i],
        sqr_btn_w,
#if MIXMAXTRIX_CHANNELS == 4
        make_array (sqr_btn_h),
        out_btns);
#elif MIXMAXTRIX_CHANNELS == 8
        make_array (sqr_btn_h, sqr_btn_h),
        out_btns[0],
        out_btns[4],
        out_btns[1],
        out_btns[5],
        out_btns[2],
        out_btns[6],
        out_btns[3],
        out_btns[7]);
#endif
    }

    // place the mixer send buttons
    area = getLocalBounds();
    area.removeFromTop (mix_snd_btn_rm_header);
    area.removeFromBottom (mix_snd_btn_rm_footer);

    std::array<juce::Rectangle<int>, parameters::n_stereo_busses + 1> btn_areas;
    for (uint i = 0; i < btn_areas.size(); ++i) {
      btn_areas[i] = area.removeFromLeft (sep_w);
      area.removeFromLeft (col_w);
    }

    auto& mixsend_buttons = p_get (parameters::mixer_sends {})[0]->buttons;
    _side_lines[0].setBounds (btn_areas[0]);
    for (uint i = 1; i < btn_areas.size() - 1; ++i) {
      mixsend_buttons[i - 1].setBounds (btn_areas[i]);
    }
    _side_lines[1].setBounds (btn_areas.back());
  }
  //----------------------------------------------------------------------------
  void paint (juce::Graphics& g) override
  {
    g.fillAll (
      getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
  }
  //----------------------------------------------------------------------------
  void on_fx_type_or_page_change (uint chnl, bool fx_change)
  {
    auto set_visible = [] (juce::Component& c) { c.setVisible (true); };
    auto set_hidden  = [] (juce::Component& c) { c.setVisible (false); };

    for (uint i = 0; i < num_fx_page_params; ++i) {
      _fx_off_sliders[chnl][i].foreach_component (set_hidden);
    }

    mp11::mp_for_each<parameters::all_fx_typelists> ([=] (auto fxtl) {
      mp11::mp_for_each<decltype (fxtl)> ([=] (auto param) {
        p_get (param)[chnl]->foreach_component (set_hidden);
      });
    });

    auto&          combo  = p_get (parameters::fx_type {})[chnl]->combo;
    auto           fx_id  = combo.getSelectedId();
    constexpr auto fx_val = parameters::fx_type {};
    if (fx_id <= 1 || fx_id > fx_val.type.choices.size()) {
      // No fx value: dummy sliders
      combo.setSelectedId (1, juce::NotificationType::dontSendNotification);
      for (uint i = 0; i < num_fx_page_params; ++i) {
        _fx_off_sliders[chnl][i].foreach_component (set_visible);
      }

      for (auto& b : _page_buttons[chnl]) {
        b.setToggleState (false, juce::NotificationType::dontSendNotification);
        b.setEnabled (false);
      }
      return;
    }

    static_assert (num_fx_pages <= 2, "");
    uint page_idx      = _page_buttons[chnl][1].getToggleState();
    uint param_idx_beg = page_idx * num_fx_page_params;
    uint param_idx_end = param_idx_beg + num_fx_page_params;

    uint fx_idx = 0;
    mp11::mp_for_each<parameters::all_fx_typelists> ([=, &fx_idx] (auto fxtl) {
      // 2 = combobox unset (juce) + no fx value.
      if (fx_idx++ != fx_id - 2) {
        return; // next
      }

      uint p_idx  = 0;
      uint p_last = 0;
      mp11::mp_for_each<decltype (fxtl)> ([=, &p_idx, &p_last] (auto param) {
        if (p_idx >= param_idx_beg && p_idx < param_idx_end) {
          p_get (param)[chnl]->foreach_component (set_visible);
          ++p_last;
        }
        ++p_idx;
      });
      for (; (param_idx_beg + p_last) < param_idx_end; ++p_last) {
        _fx_off_sliders[chnl][p_last % num_fx_page_params].foreach_component (
          set_visible);
      }
      static_assert (num_fx_pages <= 2, "");
      _page_buttons[chnl][0].setEnabled (true);
      _page_buttons[chnl][1].setEnabled (true);

      if (fx_change) {
        _page_buttons[chnl][0].setToggleState (
          true, juce::NotificationType::sendNotification);
      }
      if (p_idx <= num_fx_page_params) {
        _page_buttons[chnl][0].setToggleState (
          true, juce::NotificationType::sendNotification);
        _page_buttons[chnl][1].setEnabled (false);
        _page_buttons[chnl][1].setToggleState (
          false, juce::NotificationType::dontSendNotification);
      }
    });
  }
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
    _parameter_value.setText (
      text, juce::NotificationType::dontSendNotification);
  }
  //----------------------------------------------------------------------------
  void register_mouse_events()
  {
    using namespace std::placeholders;

    pforeach ([this] (auto key, auto& warray) {
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
  bool fx_drag_and_drop_event (
    uint dst_chnl_idx,
    event_common::source,
    juce::DragAndDropTarget::SourceDetails const& drag_src,
    drop_event::type                              t)
  {
    using namespace drop_event;
    if (t != type::request_interest && t != type::dropped) {
      return false;
    }
    juce::Component* dd_src = drag_src.sourceComponent.get();
    if (!dd_src) {
      return false;
    }

    auto& fx_type_arr  = p_get (parameters::fx_type {});
    uint  src_chnl_idx = 0;
    for (; src_chnl_idx < fx_type_arr.size(); ++src_chnl_idx) {
      juce::Component* c = &fx_type_arr[src_chnl_idx]->combo;
      if (dd_src == c) {
        break;
      }
    }
    if (src_chnl_idx == fx_type_arr.size() || src_chnl_idx == dst_chnl_idx) {
      return false;
    }
    if (t == type::request_interest) {
      return true;
    }
    if (t == type::dropped) {
      bool ctrl = juce::ModifierKeys::getCurrentModifiers().isCtrlDown();
      bool alt  = juce::ModifierKeys::getCurrentModifiers().isAltDown();
      if (ctrl && alt) {
        slider_list_cp_or_swap (
          src_chnl_idx,
          dst_chnl_idx,
          true,
          parameters::channel_sliders_typelist {});
        ctrl = false; // swap
      }
      channel_fx_cp_or_swap (src_chnl_idx, dst_chnl_idx, !ctrl);
      return false;
    }
    return false;
  };
  //----------------------------------------------------------------------------
  void channel_fx_cp_or_swap (uint src_chnl, uint dst_chnl, bool swap)
  {
    slider_list_cp_or_swap (
      src_chnl, dst_chnl, swap, parameters::channel_fx_sliders_typelist {});

    constexpr auto notif = juce::NotificationType::sendNotificationAsync;

    auto&           fx_type = p_get (parameters::fx_type {});
    juce::ComboBox& src     = fx_type[src_chnl]->combo;
    juce::ComboBox& dst     = fx_type[dst_chnl]->combo;
    int             src_val = src.getSelectedId();
    if (swap) {
      src.setSelectedId (dst.getSelectedId(), notif);
    }
    dst.setSelectedId (src_val, notif);
  }
  //----------------------------------------------------------------------------
  template <class slider_typelist>
  void slider_list_cp_or_swap (
    uint src_chnl,
    uint dst_chnl,
    bool swap,
    slider_typelist)
  {
    constexpr auto notif = juce::NotificationType::sendNotificationAsync;

    mp11::mp_for_each<slider_typelist> ([=] (auto paramtype) {
      auto&         param = p_get (paramtype);
      juce::Slider& src   = param[src_chnl]->slider;
      juce::Slider& dst   = param[dst_chnl]->slider;

      double dst_val = dst.getValue();
      dst.setValue (src.getValue(), notif);
      if (swap) {
        src.setValue (dst_val, notif);
      }
    });
  }
  //----------------------------------------------------------------------------
  void show_about_popup()
  {
    juce::DialogWindow::LaunchOptions lo;

    lo.dialogTitle                  = "About";
    lo.escapeKeyTriggersCloseButton = true;
    lo.resizable                    = false;
    lo.useNativeTitleBar            = false;
    lo.useBottomRightCornerResizer  = false;
    lo.componentToCentreAround      = this;

    lo.content.setOwned (new juce::TextEditor {});
    auto& te = static_cast<juce::TextEditor&> (*lo.content);
    te.setBounds (this->getBounds().reduced (40));
    te.setMultiLine (true);
    te.setReadOnly (true);
    te.setScrollbarsShown (true);
    te.setText (parameters::about_text, juce::dontSendNotification);
    te.setLookAndFeel (&_lookfeel);
    te.getFont().getStyleFlags();
    // TODO: font size not working...
    auto cf = te.getFont();
    te.setFont (juce::Font (
      cf.getTypefaceName(), cf.getTypefaceStyle(), cf.getHeight() * 3));

    if (auto wdw = lo.launchAsync()) {
      wdw->setLookAndFeel (&_lookfeel);
    }
  }
  //----------------------------------------------------------------------------
private:
  static constexpr uint num_fx_params      = 16;
  static constexpr uint num_fx_pages       = 2;
  static constexpr uint num_fx_page_params = num_fx_params / num_fx_pages;

  juce::AudioProcessor& _processor; // unused
  look_and_feel         _lookfeel;

  using chnl_slider_array = std::array<slider_ext, num_fx_page_params>;
  using chnl_fx_page_btns
    = std::array<add_juce_callbacks<juce::TextButton>, num_fx_pages>;

  std::array<chnl_slider_array, parameters::n_stereo_busses> _fx_off_sliders;
  std::array<chnl_fx_page_btns, parameters::n_stereo_busses> _page_buttons;

  juce::Label                 _parameter_value;
  juce::TextButton            _parameter_value_frame;
  juce::TextButton            _about;
  std::array<VerticalLine, 2> _side_lines;
  std::array<VerticalLine, parameters::n_stereo_busses * 2> _fx_lines;

  std::array<foleys::LevelMeterLookAndFeel, parameters::n_stereo_busses>
                                                              _meters_look_feel;
  std::array<foleys::LevelMeter, parameters::n_stereo_busses> _meters;
  crange<foleys::LevelMeterSource>                            _meter_srcs;

}; // namespace artv
//------------------------------------------------------------------------------
juce::AudioProcessorEditor* new_editor (
  juce::AudioProcessor&               p,
  juce::AudioProcessorValueTreeState& params,
  crange<foleys::LevelMeterSource>    meter_srcs)
{
  return new artv::editor (p, params, meter_srcs);
}
// -----------------------------------------------------------------------------
} // namespace artv
// -----------------------------------------------------------------------------
