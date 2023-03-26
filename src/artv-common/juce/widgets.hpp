#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <juce_audio_processors/juce_audio_processors.h>
#include <optional>
#include <set>
#include <type_traits>
#include <variant>
#include <vector>
// #include <stdio.h>

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
// -----------------------------------------------------------------------------
namespace event_common {
// this one uses pointers because "reference_wrappers" via cref would inhibit
// implicit conversion to base classes.
using source = std::variant<
  std::monostate,
  juce::Slider const*,
  juce::Button const*,
  juce::ComboBox const*>;
} // namespace event_common
// -----------------------------------------------------------------------------
namespace mouse_event {
// -----------------------------------------------------------------------------
enum type {
  move,
  enter,
  exit,
  down,
  drag,
  up,
  double_click,
  wheel_move,
  magnify
};

using data = std::variant<
  std::monostate,
  std::reference_wrapper<const juce::MouseWheelDetails>,
  float>; // maybe this needs a tagged class if more floats are added.
// -----------------------------------------------------------------------------
} // namespace mouse_event
// -----------------------------------------------------------------------------
namespace drop_event {
// -----------------------------------------------------------------------------
enum type {
  request_interest,
  enter,
  move,
  exit,
  dropped,
};
// -----------------------------------------------------------------------------
} // namespace drop_event
// Move from an inheritance-based way of working to a callback based one.
template <class T>
struct add_juce_callbacks : public T, public juce::DragAndDropTarget {

  static_assert (std::is_base_of<juce::Component, T>::value, "");

  using T::T;
  // mouse wrapping ------------------------------------------------------------
  std::function<void (
    event_common::source,
    mouse_event::type,
    juce::MouseEvent const&,
    mouse_event::data)>
    on_mouse_event;

  void mouseMove (juce::MouseEvent const& event) override
  {
    T::mouseMove (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::move, event, mouse_event::data {});
    }
  }

  void mouseEnter (juce::MouseEvent const& event) override
  {
    T::mouseEnter (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::enter, event, mouse_event::data {});
    }
  }

  void mouseExit (juce::MouseEvent const& event) override
  {
    T::mouseExit (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::exit, event, mouse_event::data {});
    }
  }

  void mouseDown (juce::MouseEvent const& event) override
  {
    T::mouseDown (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::down, event, mouse_event::data {});
    }
  }

  void mouseDrag (juce::MouseEvent const& event) override
  {
    T::mouseDrag (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::drag, event, mouse_event::data {});
    }
    if (!is_drag_source) {
      return;
    }
    juce::Component* self = static_cast<T*> (this);

    juce::DragAndDropContainer* parent
      = juce::DragAndDropContainer::findParentDragContainerFor (self);

    if (parent && !parent->isDragAndDropActive()) {
      parent->startDragging (self->getName(), self);
    }
  }

  void mouseUp (juce::MouseEvent const& event) override
  {
    T::mouseUp (event);
    if (on_mouse_event) {
      on_mouse_event (this, mouse_event::type::up, event, mouse_event::data {});
    }
  }

  void mouseDoubleClick (juce::MouseEvent const& event) override
  {
    T::mouseDoubleClick (event);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::double_click, event, mouse_event::data {});
    }
  }

  void mouseWheelMove (
    juce::MouseEvent const&        event,
    juce::MouseWheelDetails const& wheel) override
  {
    T::mouseWheelMove (event, wheel);
    if (on_mouse_event) {
      on_mouse_event (
        this, mouse_event::type::wheel_move, event, std::cref (wheel));
    }
  }

  void mouseMagnify (juce::MouseEvent const& event, float scaleFactor) override
  {
    T::mouseMagnify (event, scaleFactor);
    if (on_mouse_event) {
      on_mouse_event (this, mouse_event::type::magnify, event, scaleFactor);
    }
  }
  // drag and drop -------------------------------------------------------------
  bool is_drag_source = false; // enable drag

  std::function<bool (
    event_common::source, // component on which the event is happening
    juce::DragAndDropTarget::SourceDetails const&,
    drop_event::type)>
    on_drop {};

  bool isInterestedInDragSource (
    juce::DragAndDropTarget::SourceDetails const& src) final
  {
    // what to do?
    bool decide = false;
    bool ret    = isInterestedInDragSource (src, decide);
    if (!decide) {
      return ret;
    }
    return on_drop ? on_drop (this, src, drop_event::type::request_interest)
                   : false;
  }

  virtual bool isInterestedInDragSource (
    juce::DragAndDropTarget::SourceDetails const& src,
    bool&                                         let_function_handler_decide)
  {
    let_function_handler_decide = true;
    return false;
  }

  void itemDragEnter (
    juce::DragAndDropTarget::SourceDetails const& src) override
  {
    if (on_drop) {
      on_drop (this, src, drop_event::type::enter);
    }
  }

  void itemDragMove (juce::DragAndDropTarget::SourceDetails const& src) override
  {
    if (on_drop) {
      on_drop (this, src, drop_event::type::move);
    }
  }

  void itemDragExit (juce::DragAndDropTarget::SourceDetails const& src) override
  {
    if (on_drop) {
      on_drop (this, src, drop_event::type::exit);
    }
  }

  void itemDropped (juce::DragAndDropTarget::SourceDetails const& src) override
  {
    if (on_drop) {
      on_drop (this, src, drop_event::type::dropped);
    }
  }
};
// -----------------------------------------------------------------------------
namespace detail {
// -----------------------------------------------------------------------------
template <class T>
struct has_attachment {
  using attachment_type = T;

  template <class U>
  void init_attachment (
    char const*                         param_id,
    juce::AudioProcessorValueTreeState& params,
    U&                                  component)
  {
    assert (params.getParameter (param_id) && "Bug: unknown id!");
    attachment.emplace (params, param_id, component);
  }

  void clear_attachment() { attachment.reset(); }

  std::optional<T> attachment;
};
// -----------------------------------------------------------------------------
using has_slider_attachment
  = has_attachment<juce::AudioProcessorValueTreeState::SliderAttachment>;

using has_button_attachment
  = has_attachment<juce::AudioProcessorValueTreeState::ButtonAttachment>;

using has_combobox_attachment
  = has_attachment<juce::AudioProcessorValueTreeState::ComboBoxAttachment>;
// -----------------------------------------------------------------------------
template <class Derived>
struct widget_base {
  template <class Functor>
  void foreach_component (Functor f)
  {
    for (auto v : static_cast<Derived*> (this)->get_components()) {
      f (*static_cast<juce::Component*> (v));
    }
  }
};
// -----------------------------------------------------------------------------
// having to inherit everything feels so 90's...
// A window that when modal gets destroyed by clicking outside
class resizable_win_destroyed_clicking : public juce::ResizableWindow {
public:
  using ResizableWindow::ResizableWindow;

  virtual void inputAttemptWhenModal() override
  {
    setVisible (false);
    exitModalState (0);
  }
};
// A slider that allows data entry on the center of the slider and doesn't have
// the width limited to the actual slider size.
class slider_w_data_entry : public juce::Slider {
public:
  template <class... Ts>
  slider_w_data_entry (Ts&&... args)
  {
    _edit.setScrollbarsShown (false);
    _edit.setReadOnly (false);
    _edit.onEscapeKey = [this]() {
      _edit_win->setVisible (false);
      _edit_win->exitModalState (0);
    };
    _edit.onReturnKey = [this]() {
      setValue (getValueFromText (_edit.getText()));
      _edit.onEscapeKey();
    };
    _edit.onFocusLost = [this]() { _edit.onEscapeKey(); };
  }

  ~slider_w_data_entry()
  {
    if (_edit_win) {
      _edit_win->setLookAndFeel (nullptr);
    }
    _edit.setLookAndFeel (nullptr);
  }

  void lookAndFeelChanged() override
  {
    Slider::lookAndFeelChanged();
    _edit.setLookAndFeel (&getLookAndFeel());
    if (_edit_win) {
      _edit_win->setLookAndFeel (&getLookAndFeel());
    }
  }

  void parentHierarchyChanged() override
  {
    // make the window be children of the main window
    Slider::parentHierarchyChanged();
    auto top = getTopLevelComponent();
    if (top && _edit_win) {
      top->addChildComponent (*_edit_win);
    }
  }

  void mouseDown (juce::MouseEvent const& e) override
  {
    juce::ModifierKeys mods = juce::ModifierKeys::getCurrentModifiersRealtime();
    if (mods.isRightButtonDown() && isEnabled()) {
      // lazy creation, Windows are expensive, at least on Linux
      if (!_edit_win) {
        _edit_win.emplace ("", true);
        _edit_win->setLookAndFeel (&getLookAndFeel());
        if (auto top = getTopLevelComponent()) {
          top->addChildComponent (*_edit_win);
        }
      }
      adjust_positions();
      _edit_win->setContentNonOwned (&_edit, false);
      _edit_win->setVisible (true);
      _edit_win->enterModalState (true, nullptr, false);
      _edit.grabKeyboardFocus();
    }
    else {
      Slider::mouseDown (e);
    }
  }

private:
  void adjust_positions()
  {
    auto b = getBounds();
    _edit.setText (
      filter_number (getTextFromValue (getValue())),
      juce::dontSendNotification);
    b = b.withSizeKeepingCentre (
      std::max (
        b.getWidth(),
        _edit.getTextWidth() + _edit.getBorder().getLeftAndRight()),
      _edit.getTextHeight());
    if (b.getX() < 0) {
      b = b.withX (0);
    }
    auto top = getTopLevelComponent();
    if (top) {
      auto rdiff = top->getLocalBounds().getRight() - b.getRight();
      if (rdiff < 0) {
        b = b.withX (b.getX() + rdiff);
      }
    }
    _edit_win->setBounds (b);
    _edit.setBounds (b);
  }

  static juce::String filter_number (juce::String const& s)
  {
    bool has_dot = false;
    int  last    = 0;
    for (; last < s.length(); ++last) {
      if (s[last] == '.') {
        has_dot = true;
        continue;
      }
      if (!std::isdigit (s[last])) {
        if (s[last] == '-' && last == 0) {
          continue;
        }
        break;
      }
    }
    if (!has_dot) {
      // not a number.
      return s;
    }
    auto num  = s.substring (0, last);
    auto tail = s.substring (last, s.length());
    if (tail.startsWith ("kHz")) {
      // hack, only for Hz
      num = juce::String {num.getDoubleValue() * 1000.};
    }
    return num;
  }

  std::optional<resizable_win_destroyed_clicking> _edit_win;
  juce::TextEditor                                _edit;
};

// -----------------------------------------------------------------------------
} // namespace detail
// -----------------------------------------------------------------------------
struct slider_ext
  : public detail::has_slider_attachment,
    private juce::Slider::Listener,
    public detail::widget_base<slider_ext> {

  ~slider_ext()
  {
    // has to be destroyed before the gui elems it points to.
    clear_attachment();
  }

  void init (
    char const*                         param_id,
    char const*                         text,
    char const*                         suffix,
    double                              max_value,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    common_init (text, suffix, parent);
    _paramslider.addListener (this);
    init_attachment (param_id, params, _paramslider);
    slider.addListener (this);
    copy_paramslider_settings (max_value);
    sliderValueChanged (&_paramslider);
  }

  // detached, used for disabled sliders, while still being able to access
  // everything provided by slider_ex
  void init (char const* text, char const* suffix, juce::Component& parent)
  {
    common_init (text, suffix, parent);
  }

  void clear()
  {
    slider.setLookAndFeel (nullptr);
    label.setLookAndFeel (nullptr);
    _paramslider.setLookAndFeel (nullptr);
  }

  std::array<juce::Component*, 2> get_components() { return {&slider, &label}; }

  add_juce_callbacks<detail::slider_w_data_entry> slider;
  juce::Label                                     label;

private:
  void sliderValueChanged (juce::Slider* ptr) final
  {
    if (_feedback) {
      // we can't use direct value comparison, as both components may have
      // different ranges;
      return;
    }
    _feedback        = true;
    bool host_change = (ptr == &_paramslider);
    if (host_change) {
      double v = _paramslider.getValue();
      if (v > slider.getMaximum()) {
        v = slider.getDoubleClickReturnValue();
      }
      slider.setValue (v, juce::sendNotificationSync);
    }
    else {
      _paramslider.setValue (slider.getValue(), juce::sendNotificationSync);
    }
    _feedback = false;
  }

  void common_init (
    char const*      text,
    char const*      suffix,
    juce::Component& parent)
  {
    slider.setName (text);
    slider.setTextValueSuffix (suffix);
    slider.setSliderStyle (juce::Slider::RotaryVerticalDrag);
    slider.setTextBoxStyle (juce::Slider::NoTextBox, false, 0, 0);

    // slider.setPopupDisplayEnabled (true, true, &parent);
    slider.setLookAndFeel (&parent.getLookAndFeel());

    parent.addAndMakeVisible (slider);

    label.setText (text, juce::dontSendNotification);
    label.setJustificationType (juce::Justification::centred);
    label.setLookAndFeel (&parent.getLookAndFeel());
    parent.addAndMakeVisible (label);
  }

  void copy_paramslider_settings (double max_value)
  {
    auto& ps = _paramslider;
    slider.setRange (ps.getMinimum(), max_value, ps.getInterval());
    slider.setDoubleClickReturnValue (true, ps.getDoubleClickReturnValue());
    slider.setSkewFactor (ps.getSkewFactor(), ps.isSymmetricSkew());
    slider.addListener (this);
    slider.valueFromTextFunction = [=] (juce::String const& v) {
      return _paramslider.getValueFromText (v);
    };
    slider.textFromValueFunction
      = [=] (double v) { return _paramslider.getTextFromValue (v); };
  }

  // this allows us to have a different range on the real parameter than on the
  // GUI. Just done for choice parameters.
  juce::Slider _paramslider;
  bool         _feedback = false;
};

template <class T>
using is_slider_ext = std::is_same<T, slider_ext>;
// -----------------------------------------------------------------------------
template <class T, bool Toggle = true>
struct button_ext : public detail::has_button_attachment,
                    public detail::widget_base<button_ext<T, Toggle>> {

  ~button_ext()
  {
    // has to be destroyed before the gui elems it points to.
    clear_attachment();
  }

  void init (
    char const*                         param_id,
    char const*                         text,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    button.setName (text);
    button.setButtonText (text);
    button.setLookAndFeel (&parent.getLookAndFeel());
    button.setClickingTogglesState (Toggle);
    parent.addAndMakeVisible (button);
    init_attachment (param_id, params, button);
  }

  void clear() { button.setLookAndFeel (nullptr); }

  std::array<juce::Component*, 1> get_components() { return {&button}; }

  add_juce_callbacks<T> button;
};

template <class T>
struct is_button_ext : public std::false_type {};

template <class T, bool Toggle>
struct is_button_ext<button_ext<T, Toggle>> : public std::true_type {};
// -----------------------------------------------------------------------------
struct combobox_ext
  : public detail::has_combobox_attachment,
    public juce::ComboBox::Listener,
    public detail::widget_base<combobox_ext> {
  combobox_ext() = default;
  // there is a lambda capturing "this": no copy/move
  combobox_ext (combobox_ext&&)                 = delete;
  combobox_ext (combobox_ext const&)            = delete;
  combobox_ext& operator= (combobox_ext&&)      = delete;
  combobox_ext& operator= (combobox_ext const&) = delete;

  ~combobox_ext()
  {
    // has to be destroyed before the gui elems it points to.
    clear_attachment();
  }

  void init (
    char const*                         param_id,
    char const*                         text,
    char const* const*                  choices,
    size_t                              n_choices,
    size_t                              n_future_choices,
    uint                                default_index,
    bool                                item_sort,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    assert (default_index < (n_choices + n_future_choices));
    _default_id = default_index;

    combo.setLookAndFeel (&parent.getLookAndFeel());
    init_combobox_items (choices, n_choices, n_future_choices, item_sort);

    combo.setScrollWheelEnabled (true);
    combo.setJustificationType (juce::Justification::Flags::centred);
    parent.addAndMakeVisible (combo);
    init_attachment (param_id, params, _paramcombo);

    _paramcombo.addListener (this);
    combo.addListener (this);

    combo.setName (text);
    comboBoxChanged (&_paramcombo);

    prev.setLookAndFeel (&parent.getLookAndFeel());
    prev.setButtonText ("<");
    prev.setName (juce::String {text} + " Prev");
    prev.onClick = [this] {
      int idx = combo.getSelectedItemIndex() - 1;
      do {
        if (idx < 0) {
          idx = combo.getNumItems() - 1;
        }
        if (!combo.isItemEnabled (combo.getItemId (idx))) {
          --idx;
          continue;
        }
        combo.setSelectedItemIndex (idx);
        break;
      } while (combo.getSelectedItemIndex() != idx);
    };
    parent.addAndMakeVisible (prev);

    next.setLookAndFeel (&parent.getLookAndFeel());
    next.setButtonText (">");
    next.setName (juce::String {text} + " Next");
    next.onClick = [this] {
      int idx = combo.getSelectedItemIndex() + 1;
      do {
        if (idx == combo.getNumItems()) {
          idx = 0;
        }
        if (!combo.isItemEnabled (combo.getItemId (idx))) {
          ++idx;
          continue;
        }
        combo.setSelectedItemIndex (idx);
        break;
      } while (combo.getSelectedItemIndex() != idx);
    };
    next.setConnectedEdges (juce::TextButton::ConnectedOnLeft);
    parent.addAndMakeVisible (next);

    int id = combo.getSelectedItemIndex();
    if (id == -1) {
      // only overwrite if the attachment didn't set a value
      combo.setSelectedItemIndex (_default_id);
    }
  }

  void clear()
  {
    combo.setLookAndFeel (nullptr);
    prev.setLookAndFeel (nullptr);
    next.setLookAndFeel (nullptr);
  }

  std::array<juce::Component*, 3> get_components()
  {
    return {&combo, &prev, &next};
  }

  // if this is 0, the "grid" function will consume one row for the buttons,
  // otherwise it will place them at the right according to this value. Max =
  // 1f. E.g. 0.75f makes the combobox consume 75% of the width.

  add_juce_callbacks<juce::ComboBox> combo;
  // There are more fancy buttons but they aren't default constructible, as
  // these members are public I don't want to have optionals on some components
  // sometimes...
  add_juce_callbacks<juce::TextButton> prev;
  add_juce_callbacks<juce::TextButton> next;

private:
  virtual void comboBoxChanged (juce::ComboBox* ptr) final
  {
    if (_feedback) {
      // we can't use direct value comparison, as both components may have
      // different ranges;
      return;
    }
    _feedback        = true;
    bool host_change = ptr == &_paramcombo;
    if (host_change) {
      auto id = _paramcombo.getSelectedId();
      if (id > combo.getNumItems()) { // there is an unaccounted id 0
        id = _default_id;
      }
      combo.setSelectedId (id, juce::sendNotificationSync);
    }
    else {
      _paramcombo.setSelectedId (
        combo.getSelectedId(), juce::sendNotificationSync);
    }
    _feedback = false;
  }

  struct comboitem {
    char const* str;
    uint        id;

    bool operator<(comboitem const& other)
    {
      return strcmp (str, other.str) < 0;
    }
  };

  void init_combobox_items (
    char const* const* choices,
    size_t             n_choices,
    size_t             n_future_choices,
    bool               sort)
  {
    std::set<std::string>  sections {};
    std::vector<comboitem> values;

    for (uint i = 0; i < n_choices; ++i) {
      values.emplace_back (comboitem {choices[i], i + 1});
    }
    if (sort) {
      std::sort (values.begin(), values.end());
    }

    int i = 0;
    for (; i < n_choices; ++i) {
      std::string item {values[i].str};

      if (item[0] == ':') {
        auto        sect_start = 1;
        auto        sect_end   = item.find (" ");
        std::string sect = item.substr (sect_start, sect_end - sect_start);
        item             = item.substr (sect_end + 1);
        if (sections.find (sect) == sections.end()) {
          combo.addSectionHeading (sect);
          sections.emplace (std::move (sect));
        }
      }

      combo.addItem (item, values[i].id);
      _paramcombo.addItem (choices[i], i + 1);
    }
    for (; i < (n_choices + n_future_choices); ++i) {
      _paramcombo.addItem ("reserved", i + 1);
    }
  }

  // "_paramcombo" is there only to have an easy implementation of range
  // reservation for future expansion (to don't break user automation), so
  // "_paramcombo" is linked to the APVTS and has reserved values, while "combo"
  // is linked via callbacks to "_paramcombo" but shown to the user without the
  // reserved values.
  juce::ComboBox _paramcombo;
  bool           _prev_next_enabled           = true;
  float          _prev_next_grid_width_factor = 0.f;
  uint           _default_id;
  bool           _feedback = false;
};

template <class T>
using is_combobox_ext = std::is_same<T, combobox_ext>;
//------------------------------------------------------------------------------
namespace detail {
// -----------------------------------------------------------------------------
// multiple buttons as a bitfield via using a dummy slider (for easier impl).
// Lazy but effective. For simplifying plugins with loads of plugins.
//
// Placement of the buttons has to be done externally as usual.
// -----------------------------------------------------------------------------
class toggle_buttons_impl
  : private juce::Button::Listener,
    private juce::Slider::Listener,
    public detail::has_slider_attachment,
    public detail::widget_base<toggle_buttons_impl> {
private:
  using bitfield_type = uint32_t;

public:
  toggle_buttons_impl (xspan<add_juce_callbacks<juce::TextButton>> txt_buttons)
    : buttons (txt_buttons)
  {}

  ~toggle_buttons_impl()
  {
    // has to be destroyed before the gui elems it points to.
    clear_attachment();
  }

  void init (
    char const*                         param_id,
    char const* const*                  buttons_text,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    for (size_t i = 0; i < buttons.size(); ++i) {
      auto& button = buttons[i];
      button.setName (buttons_text[i]);
      button.addListener (this);
      button.setButtonText (buttons_text[i]);
      button.setLookAndFeel (&parent.getLookAndFeel());
      button.setClickingTogglesState (true);
      parent.addAndMakeVisible (button);
    }
    _hidden_slider.addListener (this);
    init_attachment (param_id, params, _hidden_slider);
  }

  void clear()
  {
    for (auto& button : buttons) {
      button.setLookAndFeel (nullptr);
    }
  }

  void buttonClicked (juce::Button* candidate_button) final
  {
    // This could be done in O(1), with just offsetoff and pointer arith,
    // just going for the simple impl...
    bitfield_type i;
    for (i = 0; i < buttons.size(); ++i) {
      if (candidate_button == static_cast<juce::Button*> (&buttons[i])) {
        break;
      }
    }
    if (i == buttons.size()) {
      return;
    }
    auto& button   = *candidate_button;
    bool  on       = button.getToggleState();
    auto  idx_mask = bit<bitfield_type> (i);

    if (idx_mask & _on_radio) {
      if (!on) {
        // updates on "_value" are not atomic. "_value" will be/ unset
        // when this radio does the set operation on the button that
        // triggered this. Assuming only one mouse pointer, so no
        // keeping track of many values changing at once, which would
        // require a change list.
        _radio_swap_source |= idx_mask;
        return;
      }
      else {
        if (!(idx_mask & _value)) {
          // clearing the callback that was ignored on the branch
          // above, so the processor can't see an intermediate state.
          _value &= ~_radio_swap_source;
          _radio_swap_source = 0;
        }
        else if (idx_mask & _on_radio_allowing_all_unset) {
          // user clicked but the radio group prevented disabling, as
          // the button previous state hasn't changed.
          button.setToggleState (false, juce::dontSendNotification);
          on = false;
        }
      }
    }
    set_bit (_value, i, on);
    _hidden_slider.setValue ((double) _value);
  }

  void sliderValueChanged (juce::Slider* ptr) final
  {
    if (ptr != &_hidden_slider) {
      return;
    }

    bitfield_type new_value = (bitfield_type) _hidden_slider.getValue();
    // the slider range is fixed at 23 bits for allowing expanding the
    // parameters. Mask it.
    new_value &= lsb_mask<uint> (buttons.size());
    if (new_value == _value) {
      return;
    }
    bitfield_type change_list = new_value ^ _value;

    iterate_set_bits (change_list, [&] (unsigned i) {
      bitfield_type mask = bit<bitfield_type> (i);
      bool          on   = !!(new_value & mask);
      buttons[i].setToggleState (on, juce::dontSendNotification);
      set_bit (_value, i, on);
    });
  }

  void set_radio_group (
    int           juce_radio_id,
    bitfield_type first,
    bitfield_type last,
    bool          allow_all_unset = false)
  {
    assert (first <= last);
    assert (last < buttons.size());

    bitfield_type members = lsb_mask<bitfield_type> (last - first + 1);
    members <<= first;
    if (juce_radio_id != 0) {
      // assert if some members already on a radio. This implementation
      // wants to be trivial and won't track which button is on which
      // group id with allowing all unset or not.
      assert (!(members & _on_radio));
      _on_radio |= members;
      if (allow_all_unset) {
        _on_radio_allowing_all_unset |= members;
      }
    }
    else {
      _on_radio_allowing_all_unset &= ~members;
      _on_radio &= ~members;
    }
    iterate_set_bits (members, [&] (unsigned i) {
      buttons[i].setRadioGroupId (juce_radio_id);
    });
  }

  void set_radio_group (int juce_radio_id, bool allow_all_unset = false)
  {
    set_radio_group (juce_radio_id, 0, buttons.size() - 1, allow_all_unset);
  }

  // has to return ptrs...
  auto get_components() { return buttons; }

  xspan<add_juce_callbacks<juce::TextButton>> buttons;

private:
  juce::Slider  _hidden_slider;
  bitfield_type _value                       = 0;
  bitfield_type _on_radio                    = 0;
  bitfield_type _on_radio_allowing_all_unset = 0;
  bitfield_type _radio_swap_source           = 0;
};
//------------------------------------------------------------------------------
} // namespace detail

//------------------------------------------------------------------------------
template <size_t N>
class toggle_buttons : public detail::toggle_buttons_impl {
public:
  static_assert (N <= 23 && "JUCE uses float for storing(?).");
  toggle_buttons() : toggle_buttons_impl {xspan {_button_storage}} {}

  // in the current state moving or copying invalidates the base, easy to fix
  // if required
  toggle_buttons (toggle_buttons const&)            = delete;
  toggle_buttons& operator= (toggle_buttons const&) = delete;
  toggle_buttons (toggle_buttons&&)                 = delete;
  toggle_buttons& operator= (toggle_buttons&&)      = delete;

  auto get_components()
  {
    std::array<juce::Component*, N> ret;
    for (uint i = 0; i < N; ++i) {
      ret[i] = &_button_storage[i];
    }
    return ret;
  }

private:
  std::array<add_juce_callbacks<juce::TextButton>, N> _button_storage;
};

// TODO: is there a better way? Maybe all these should just take traits...
template <class T>
struct is_toggle_buttons : public std::false_type {};

template <size_t N>
struct is_toggle_buttons<toggle_buttons<N>> : public std::true_type {};
// -----------------------------------------------------------------------------
} // namespace artv
