#pragma once

// Based on
//
// https://github.com/ffAudio/ffGuiAttachments
//
// See:
// https://forum.juce.com/t/best-practice-for-storing-text-on-the-apvts/49291

#include <string>

#include <juce_core/juce_core.h>
#include <juce_data_structures/juce_data_structures.h>
#include <juce_gui_basics/juce_gui_basics.h>

namespace artv {
//------------------------------------------------------------------------------
class value_tree_label_attachment
  : public juce::Label::Listener,
    public juce::ValueTree::Listener {
public:
  //----------------------------------------------------------------------------
  void reset (
    juce::ValueTree&   tree,
    juce::Label&       label,
    juce::Identifier   property_key,
    juce::UndoManager* undo      = nullptr,
    std::string        empty_val = "")
  {
    jassert (tree.isValid()); // Don't attach an invalid valuetree!

    _tree      = tree;
    _label     = &label;
    _property  = property_key;
    _undo      = undo;
    _updating  = false;
    _empty_val = empty_val;

    if (tree.hasProperty (_property)) {
      set_label_text_from_property();
    }
    else {
      tree.setProperty (_property, _empty_val.c_str(), _undo);
      _label->setText (_empty_val.c_str(), juce::dontSendNotification);
    }
    _tree.removeListener (this);
    _tree.addListener (this);
    _label->removeListener (this);
    _label->addListener (this);
  }
  //----------------------------------------------------------------------------
  ~value_tree_label_attachment()
  {
    _tree.removeListener (this);
    if (_label) {
      _label->removeListener (this);
    }
  }
  //----------------------------------------------------------------------------
  void labelTextChanged (juce::Label* src) override
  {
    if (_updating) {
      return;
    }
    if (src == _label) {
      _updating = true;
      _tree.setProperty (_property, _label->getText(), _undo);
      _updating = false;
    }
  }
  //----------------------------------------------------------------------------
  void valueTreePropertyChanged (
    juce::ValueTree&        src,
    juce::Identifier const& property) override
  {
    if (_updating) {
      return;
    }
    if (src == _tree && _label && _property == property) {
      set_label_text_from_property();
    }
  }
  //----------------------------------------------------------------------------
  void valueTreeChildAdded (juce::ValueTree& parent, juce::ValueTree& child)
    override
  {}
  //----------------------------------------------------------------------------
  void valueTreeChildRemoved (
    juce::ValueTree& parent,
    juce::ValueTree& child,
    int              child_removal_idx) override
  {}
  //----------------------------------------------------------------------------
  void valueTreeChildOrderChanged (
    juce::ValueTree& parent,
    int              old_idx,
    int              new_idx) override
  {}
  //----------------------------------------------------------------------------
  void valueTreeParentChanged (juce::ValueTree& tree) override {}
  //----------------------------------------------------------------------------
  void valueTreeRedirected (juce::ValueTree& tree) override
  {
    // never tested...
    reset (tree, *_label, _property, _undo, _empty_val);
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void set_label_text_from_property()
  {
    auto str  = _tree.getProperty (_property).toString();
    str       = str.isEmpty() ? juce::String {_empty_val.c_str()} : str;
    _updating = true;
    _label->setText (str, juce::dontSendNotification);
    _updating = false;
  }
  //----------------------------------------------------------------------------
  juce::ValueTree                           _tree;
  juce::Component::SafePointer<juce::Label> _label;
  juce::Identifier                          _property;
  juce::UndoManager*                        _undo = nullptr;
  std::string                               _empty_val;
  bool                                      _updating = false;
};
//------------------------------------------------------------------------------
} // namespace artv
