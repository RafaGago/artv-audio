#pragma once

#include <atomic>
#include <cstddef>
#include <tuple>
#include <type_traits>

#include <juce_audio_processors/juce_audio_processors.h>

#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/juce/widgets.hpp"
#include "artv-common/misc/hana.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/xspan.hpp"

/*
How does this work:

You define FX parameters with the "parameter_cpp_class_define" macro. The macro
defines a parameter type representing a parameter value. Notice that you have
many parameter values of the same type, see the second parameter of
"parameter_cpp_class_define".

These parameter types contain parameter traits (juce id string(s), ranges, text
strings, etc). They double as parameter definitions an tag classes.

The editor has to inherit "editor_apvts_widgets" and the processor has to
inherit "has_processor_params". Both are passed the typelist with all the
parameters/tag_classes named above.

On both the processor and editor a member "boost::hana::map" is created when
inheriting "has_content_map". It has as key the tag types defined by the
"parameter_cpp_class_define" macro and contains etiher widget arrays or builtin
type arrays in the case of the editor and processor respectively.

"has_processor_params" inherits "has_content_map". It expands it with
methods to fetch parameters from an externally owned AudioProcessorTree.

Then "editor_apvts_widgets" inherits "has_content_map" too. It expands it with
methods to initialize the widgets from an AudioProcessorTree.

These two classes remove a lot of the AudioProcessorTree boilerplate and Juce
doing-too-many-things-on-construction uglyness at the expense of requiring
metaprogramming, a capable compiler, compile times and probably code size to
some extent.

Notice that on FX's with hundreds of parameters doing a lot of searches on the
"boost::hana::map" makes the compile times to grow exponentially (around 300
parameters compiling in 2 minutes). Using as small parameter maps as possible is
encouraged, so e.g. if you have a synth with three oscillators doing a frequency
parameter like below can save you boilearplate and map space:

parameter_cpp_class_define (oscillator_freq, 3, param_common(..),
param_float(...));

"boost::hana" makes metraprogramming really accessible but the compile times
might suffer if not used carefully.
*/

namespace artv {

// All these "parameter_xxx" classes have free functions defining them because
// brace initialization doesn't play well with macros, so instead of writing
// constructors we do from functions, which support type deduction on C++14.
// -----------------------------------------------------------------------------
// Dsp_type: the DSP class for what this parameter is for
// Dsp_param: the tag (dummy) type to access the parameter on the DSP type.
template <class Dsp_type, class Dsp_param>
struct parameter_common {
  char const* const gui_text;
  using dsp_type  = Dsp_type;
  using dsp_param = Dsp_param;
};

template <class Dsp_type = void, class Dsp_param = void>
constexpr parameter_common<Dsp_type, Dsp_param> param_common (
  char const* gui_text = "",
  Dsp_type*   dt       = nullptr,
  Dsp_param*  dp       = nullptr)
{
  // passing pointers is the best I could come up with to be able to pass
  // classes with non constexpr constructors as parameters by deducing the type
  // on a constexpr context.
  return {gui_text};
}
// -----------------------------------------------------------------------------
namespace detail {
// -----------------------------------------------------------------------------
template <class derived>
struct float_parameter_base {
  static auto make_parameter()
  {
    derived d {};

    static_assert (
      is_same_template_v<parameter_common, decltype (d.common)>, "");

    juce::NormalisableRange<float> range {
      d.type.min,
      d.type.max,
      d.type.interval,
      d.type.skew,
      d.type.simmetric_skew};

    using vsw = typename decltype (derived::type)::value_string_wrapper;

    std::array<std::unique_ptr<juce::AudioParameterFloat>, d.juce_ids.size()>
      ret;
    for (uint i = 0; i < d.juce_ids.size(); ++i) {
      ret[i] = std::make_unique<juce::AudioParameterFloat> (
        d.juce_ids[i],
        d.juce_ids[i],
        range,
        d.type.defaultv,
        juce::String(),
        juce::AudioProcessorParameter::genericParameter,
        vsw::string_from_value::get(),
        vsw::value_from_string::get());
    }
    return ret;
  }
};
// -----------------------------------------------------------------------------
template <class derived>
struct int_parameter_base {
  static auto make_parameter()
  {
    derived d {};
    static_assert (
      is_same_template_v<parameter_common, decltype (d.common)>, "");

    using vsw = typename decltype (derived::type)::value_string_wrapper;

    std::array<std::unique_ptr<juce::AudioParameterInt>, d.juce_ids.size()> ret;
    for (uint i = 0; i < d.juce_ids.size(); ++i) {
      ret[i] = std::make_unique<juce::AudioParameterInt> (
        d.juce_ids[i],
        d.juce_ids[i],
        d.type.min,
        d.type.max,
        d.type.defaultv,
        juce::String(),
        vsw::string_from_value::get(),
        vsw::value_from_string::get());
    }
    return ret;
  }
};
// -----------------------------------------------------------------------------
template <class derived>
struct choice_parameter_base {
  static auto make_parameter()
  {
    derived d {};
    static_assert (
      is_same_template_v<parameter_common, decltype (d.common)>, "");

    juce::StringArray choices;
    uint              i;
    for (i = 0; i < d.type.choices.size(); ++i) {
      choices.add (derived::type.choices[i]);
    }
    for (; i < d.type.max; ++i) {
      choices.add ("future choice");
    }

    std::array<std::unique_ptr<juce::AudioParameterChoice>, d.juce_ids.size()>
      ret;
    for (uint i = 0; i < d.juce_ids.size(); ++i) {
      ret[i] = std::make_unique<juce::AudioParameterChoice> (
        d.juce_ids[i], d.juce_ids[i], choices, d.type.defaultv);
    }
    return ret;
  }
};
//------------------------------------------------------------------------------
template <class key_, class value_>
struct pair_with_key_as_type {
  using key = key_;
  value_ value;
};
// -----------------------------------------------------------------------------
template <class key>
struct has_comptime_key {
  template <class T>
  struct trait : public std::is_same<typename T::key, key> {};
};
// -----------------------------------------------------------------------------
struct take_widget_from_user {};

template <class T>
struct get_parameter_base_class;

template <class ValueStringWrapper>
struct get_parameter_base_class<float_parameter_ext<ValueStringWrapper>> {
  template <class T>
  using base         = float_parameter_base<T>;
  using fixed_widget = take_widget_from_user;
};

template <class ValueStringWrapper>
struct get_parameter_base_class<int_parameter_ext<ValueStringWrapper>> {
  template <class T>
  using base         = int_parameter_base<T>;
  using fixed_widget = take_widget_from_user;
};

template <size_t N>
struct get_parameter_base_class<choice_parameter<N>> {
  template <class T>
  using base         = choice_parameter_base<T>;
  using fixed_widget = take_widget_from_user;
};

template <size_t N>
struct get_parameter_base_class<toggle_buttons_parameter<N>> {
  template <class T>
  using base         = int_parameter_base<T>;
  using fixed_widget = toggle_buttons<N>;
};

template <class T>
using parameter_define_base
  = detail::get_parameter_base_class<std::remove_const_t<T>>;
// -----------------------------------------------------------------------------
} // namespace detail
// -----------------------------------------------------------------------------
// Macro to define a unique type for a parameter.
//
// "type_name" is the name that this type will have.
//
// "num_elems" is the number of elements that this type will contain. Useful
//  to save map space and e.g. if you have many parameters that are equal. e.g.
//  on a synth with 2 filters you might have 2 cutoff frequency parameters.
//
// "common_v" is a constexpr instance of "parameter_common".
//
// "specific_v" is constexpr instance of any of the parameter types, e.g.
//    "float_parameter", "int_parameter"...
//
// "user_widget_type" is the desired Juce widget type. On some types it is
// ignored, as they only work on one widget type.

#define parameter_cpp_class_define( \
  type_name, num_elems, common_v, specific_v, user_widget_type) \
\
  namespace juce_id_prefixes { \
  constexpr auto type_name = BOOST_HANA_STRING (#type_name); \
  } \
\
  struct type_name : public detail::parameter_define_base< \
                       decltype (specific_v)>::base<type_name> { \
\
    static constexpr std::array<char const*, num_elems> juce_ids \
      = array_of_str_with_numeric_suffix ( \
        juce_id_prefixes::type_name, \
        std::make_index_sequence<num_elems> {}); \
    static constexpr parameter_common      common {common_v}; \
    static constexpr decltype (specific_v) type {specific_v}; \
\
    using fixed_widget \
      = detail::parameter_define_base<decltype (specific_v)>::fixed_widget; \
\
    using widget_type = std::conditional_t< \
      std::is_same<fixed_widget, detail::take_widget_from_user>::value, \
      user_widget_type, \
      fixed_widget>; \
  }
// -----------------------------------------------------------------------------
template <class... types>
class has_content_map;

template <class... Keys, class... Values>
class has_content_map<mp_list<Keys...>, mp_list<Values...>> {
public:
  //----------------------------------------------------------------------------
  template <class T>
  auto& p_get()
  {
    return _content[hana::type_c<T>];
  }
  //----------------------------------------------------------------------------
  template <class T>
  auto& p_get (T)
  {
    return p_get<T>();
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  auto p_get (uint index, mp_list<Ts...> lst)
  {
    using first_type = mp11::mp_front<mp_list<Ts...>>;
    using ret_type   = std::remove_reference_t<decltype (p_get<first_type>())>;

    ret_type* ret = nullptr;
    pinvoke_on_index (index, lst, [&] (auto key, auto& value) {
      ret = &value;
    });
    return ret;
  }
  //----------------------------------------------------------------------------
  template <class F>
  void pforeach (F func)
  {
    hana::for_each (_content, [&] (auto& v) {
      using Key = typename decltype (+hana::first (v))::type;
      func (Key {}, hana::second (v));
    });
  }
  //----------------------------------------------------------------------------
  template <class... Ts, class F>
  void pforeach (mp_list<Ts...> tl, F func)
  {
    mp11::mp_for_each<mp_list<Ts...>> ([&] (auto type) {
      func (type, _content[hana::type_c<decltype (type)>]);
    });
  }
  //----------------------------------------------------------------------------
  template <class F, class... Ts>
  void pinvoke_on_index (uint index, mp_list<Ts...>, F func)
  {
    uint i = 0;
    mp11::mp_for_each<mp_list<Ts...>> ([&] (auto type) {
      if (index == i) {
        func (type, _content[hana::type_c<decltype (type)>]);
      }
      ++i;
    });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <
    class Key,
    class T,
    std::enable_if_t<std::is_move_assignable<T>::value>* = nullptr>
  static auto construct_and_maybe_ptr_wrap()
  {
    return std::array<T, Key::juce_ids.size()> {};
  }
  //----------------------------------------------------------------------------
  template <
    class Key,
    class T,
    std::enable_if_t<!std::is_move_assignable<T>::value>* = nullptr>
  static auto construct_and_maybe_ptr_wrap()
  {
    std::array<std::unique_ptr<T>, Key::juce_ids.size()> ret;
    for (auto& elem : ret) {
      elem = std::make_unique<T>();
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  static constexpr auto make_map()
  {
    return hana::make_map (hana::make_pair (
      hana::type_c<Keys>, construct_and_maybe_ptr_wrap<Keys, Values>())...);
  }

  decltype (make_map()) _content {make_map()};
};
// -----------------------------------------------------------------------------
// class for the processor, takes a mp_list with the unique types
// describing the parameters (defined with the "parameter_cpp_class_define"
// macro) and creates a member compile time map containing the parameters
// themselves ,so "fetch_editor_values" can update the values. As every
// parameter is a type, the parameter values can be accessed bgy p_get<type>()
// too.

template <class... unique_types>
class has_processor_params;

template <class... types>
class has_processor_params<mp_list<types...>>
  : public has_content_map<
      mp_list<types...>,
      mp_list<typename decltype (types::type)::value_type...>> {
private:
  //----------------------------------------------------------------------------
  static constexpr auto get_default_visitor()
  {
    return [] (auto key, uint arr_idx, auto refresh_result_v) {};
  }

public:
  //----------------------------------------------------------------------------
  template <class T>
  struct refresh_result {
    using value_type = T;
    bool       changed() const { return previous != current; }
    value_type previous, current;
  };
  //----------------------------------------------------------------------------
  static juce::AudioProcessorValueTreeState::ParameterLayout make_apvts_layout()
  {
    juce::AudioProcessorValueTreeState::ParameterLayout ret;

    mp11::mp_for_each<mp_list<types...>> ([&ret] (auto param) {
      auto arr = param.make_parameter();
      ret.add (arr.begin(), arr.end());
    });
    // TODO: this is a static, it has to initialize the ptrs...

    return ret;
  }
  //----------------------------------------------------------------------------
  void init_aptvs_references (juce::AudioProcessorValueTreeState& params)
  {
    mp11::mp_for_each<mp_list<types...>> ([this, &params] (auto key) {
      auto& ptr_arr = _atomic_ptrs.p_get (key);
      for (uint i = 0; i < ptr_arr.size(); ++i) {
        ptr_arr[i] = params.getRawParameterValue (key.juce_ids[i]);
      }
    });
  }
  //----------------------------------------------------------------------------
  // fetches a single parameter from the juce::APVTS
  template <class T>
  auto p_refresh (uint param_idx)
  {
    return do_fetch (T {}, this->template p_get<T>(), param_idx);
  }
  //----------------------------------------------------------------------------
  // fetches a single parameter from the juce::APVTS
  template <class T>
  auto p_refresh (T, uint param_idx)
  {
    return p_refresh<T> (param_idx);
  }
  //----------------------------------------------------------------------------
  // fetches all array values on a single tag parameter type
  template <class T, class V = decltype (get_default_visitor())>
  void p_refresh_many (V visitor = get_default_visitor())
  {
    auto& varray = this->template p_get<T>();
    for (uint i = 0; i < varray.size(); ++i) {
      visitor (T {}, i, do_fetch (T {}, varray, i));
    }
  }
  //----------------------------------------------------------------------------
  // fetches all array values on a single tag parameter type
  template <class T, class V = decltype (get_default_visitor())>
  void p_refresh_many (T, V visitor = get_default_visitor())
  {
    p_refresh_many<T> (visitor);
  }
  //----------------------------------------------------------------------------
  // fetches all array values on a multiple tag parameter type
  template <class... Ts, class V = decltype (get_default_visitor())>
  void p_refresh_many (
    mp_list<Ts...> fetch_list,
    V              visitor = get_default_visitor())
  {
    this->template pforeach (
      fetch_list, [=, &visitor] (auto key, auto& varray) {
        for (uint i = 0; i < varray.size(); ++i) {
          visitor (key, i, do_fetch (key, varray, i));
        }
      });
  }
  //----------------------------------------------------------------------------
  // fetches all keys and values from the juce::APVTS.
  template <class V = decltype (get_default_visitor())>
  void p_refresh_all (V visitor = get_default_visitor())
  {
    this->template pforeach ([this, &visitor] (auto key, auto& varray) {
      for (uint i = 0; i < varray.size(); ++i) {
        visitor (key, i, do_fetch (key, varray, i));
      }
    });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T, class U, size_t N>
  auto do_fetch (T key, std::array<U, N>& values, uint idx)
  {
    auto  previous       = values[idx];
    auto& atomic_ptr_arr = _atomic_ptrs.p_get (key);
    values[idx] = (U) atomic_ptr_arr[idx]->load (std::memory_order_relaxed);
    return refresh_result<decltype (previous)> {previous, values[idx]};
  }
  //----------------------------------------------------------------------------
  has_content_map<
    mp_list<types...>,
    mp11::mp_fill<mp_list<types...>, std::atomic<float>*>>
    _atomic_ptrs;
  //----------------------------------------------------------------------------
};
// -----------------------------------------------------------------------------
// class for the processor, takes a mp_list with the unique types
// describing the parameters (defined with the "parameter_cpp_class_define"
// macro) and creates a member compile time map containing the widgets. As every
// parameter is a type, the widget arrays can be accessed by p_get<type>() too.
template <class... unique_types>
class editor_apvts_widgets;

template <class... types>
class editor_apvts_widgets<mp_list<types...>>
  : public has_content_map<
      mp_list<types...>,
      mp_list<typename types::widget_type...>> {
public:
  //----------------------------------------------------------------------------
  void init_widgets (
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    this->pforeach ([&, this] (auto key, auto& warray) {
      for (uint i = 0; i < warray.size(); ++i) {
        init_widget (key, key.type, *warray[i], parent, params, i);
      }
    });
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void init_widgets (
    mp_list<Ts...>,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    mp11::mp_for_each<mp_list<Ts...>> ([&] (auto key) {
      auto& warray = this->p_get (key);
      for (uint i = 0; i < warray.size(); ++i) {
        init_widget (key, key.type, *warray[i], parent, params, i);
      }
    });
  }
  //----------------------------------------------------------------------------
  template <class... Ts>
  void init_widgets (
    mp_list<Ts...>,
    uint                                array_index,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params)
  {
    mp11::mp_for_each<mp_list<Ts...>> ([&] (auto key) {
      auto& warray = this->p_get (key);
      init_widget (
        key, key.type, *warray[array_index], parent, params, array_index);
    });
  }
  //----------------------------------------------------------------------------
  void clear_widgets()
  {
    this->pforeach ([] (auto key, auto& warray) {
      for (auto& w : warray) {
        w->clear();
      }
    });
  }
  //----------------------------------------------------------------------------
private:
  template <class defs, class param_type>
  void init_widget (
    defs d,
    param_type,
    slider_ext&                         s,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params,
    uint                                juce_id_idx)
  {
    s.init (
      d.juce_ids[juce_id_idx],
      d.common.gui_text,
      d.type.units,
      (double) d.type.max,
      parent,
      params);
  }
  template <class defs, size_t N>
  void init_widget (
    defs d,
    choice_parameter<N>,
    slider_ext&                         s,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params,
    uint                                juce_id_idx)
  {
    s.init (
      d.juce_ids[juce_id_idx],
      d.common.gui_text,
      d.type.units,
      (double) (N - 1), // might have a reserved range afterwards
      parent,
      params);
  }
  //----------------------------------------------------------------------------
  template <class defs, class param_type, class T>
  void init_widget (
    defs d,
    param_type,
    button_ext<T>&                      b,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params,
    uint                                juce_id_idx)
  {
    b.init (d.juce_ids[juce_id_idx], d.common.gui_text, parent, params);
  }
  //----------------------------------------------------------------------------
  template <class defs, class param_type>
  void init_widget (
    defs d,
    param_type,
    combobox_ext&                       c,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params,
    uint                                juce_id_idx)
  {
    uint choices_now    = d.type.choices.size();
    uint choices_future = d.type.max;
    if (choices_future > choices_now) {
      choices_future -= choices_now;
    }
    else {
      choices_future = 0;
    }
    c.init (
      d.juce_ids[juce_id_idx],
      d.common.gui_text,
      d.type.choices.data(),
      choices_now,
      choices_future,
      d.type.defaultv,
      d.type.alphabetically_sorted,
      parent,
      params);
  }
  //----------------------------------------------------------------------------
  template <class defs, class param_type, size_t N>
  void init_widget (
    defs d,
    param_type,
    toggle_buttons<N>&                  tb,
    juce::Component&                    parent,
    juce::AudioProcessorValueTreeState& params,
    uint                                juce_id_idx)
  {
    tb.init (d.juce_ids[juce_id_idx], d.type.texts.data(), parent, params);
  }
  //----------------------------------------------------------------------------
};
// -----------------------------------------------------------------------------
} // namespace artv
