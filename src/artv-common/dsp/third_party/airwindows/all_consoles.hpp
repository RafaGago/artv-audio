#pragma once

#include <variant>

#include "artv-common/dsp/own/classes/plugin_context.hpp"

#include "artv-common/dsp/third_party/airwindows/console4.hpp"
#include "artv-common/dsp/third_party/airwindows/console5.hpp"
#include "artv-common/dsp/third_party/airwindows/console6.hpp"
#include "artv-common/dsp/third_party/airwindows/console7.hpp"

namespace artv { namespace airwindows {

//------------------------------------------------------------------------------
// adapter
template <size_t N>
class all_consoles {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::console;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    plugcontext = &pc;
    console_instances_get ([=] (auto& bus, auto& channels) {
      bus.reset (*plugcontext);
      for (auto& ch : channels) {
        ch.reset (*plugcontext);
      }
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    console_instances_get ([=] (auto& bus, auto& channels) {
      bus.process (outs, ins, samples);
    });
  }
  //----------------------------------------------------------------------------
  stereo_summing_processor& get_channel (uint channel)
  {
    stereo_summing_processor* ret {};
    console_instances_get ([=, &ret] (auto& bus, auto& channels) {
      assert (channel < channels.size());
      ret = &channels[channel];
    });
    return *ret;
  }
  //----------------------------------------------------------------------------
  struct drive_tag {};
  void set (drive_tag, float v)
  {
    console_instances_get ([=] (auto& bus, auto& channels) {
      using drive_tag =
        typename std::remove_reference_t<decltype (bus)>::drive_tag;

      bus.set (drive_tag {}, v);
      for (auto& ch : channels) {
        ch.set (drive_tag {}, v);
      }
    });
  }

  static constexpr auto get_parameter (drive_tag)
  {
    return float_param ("dB", -20.f, 20.f, 0.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct type_tag {};
  void set (type_tag, int v)
  {
    if (v == _console_type || v >= get_parameter (type_tag {}).choices.size()) {
      return;
    }
    _console_type = v;
    console_types_get ([=] (auto bus_tw, auto chnl_tw) {
      using bus_type      = typename decltype (bus_tw)::type;
      using channels_type = typename decltype (chnl_tw)::type;

      _bus      = bus_type {};
      _channels = channels_type {};
    });
    reset (*plugcontext);
  }

  static constexpr auto get_parameter (type_tag)
  {
    return choice_param (
      0, make_cstr_array ("Console4", "Console5", "Console6", "Console7"), 20);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<type_tag, drive_tag>;
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class Func>
  void console_types_get (Func f)
  {
    switch (_console_type) {
    case 0:
      f (
        type_wrapper<console4bus> {},
        type_wrapper<std::array<console4channel, N>> {});
      return;
    case 1:
      f (
        type_wrapper<console5bus> {},
        type_wrapper<std::array<console5channel, N>> {});
      return;
    case 2:
      f (
        type_wrapper<console6bus> {},
        type_wrapper<std::array<console6channel, N>> {});
      return;
    case 3:
      f (
        type_wrapper<console7bus> {},
        type_wrapper<std::array<console7channel, N>> {});
      return;
    default:
      assert (false);
      return;
    };
  }
  //----------------------------------------------------------------------------
  template <class Func>
  void console_instances_get (Func f)
  {
    console_types_get ([&, this] (auto bus_tw, auto chnl_tw) {
      using bus_type      = typename decltype (bus_tw)::type;
      using channels_type = typename decltype (chnl_tw)::type;

      f (std::get<bus_type> (_bus), std::get<channels_type> (_channels));
    });
  }
  //----------------------------------------------------------------------------
  using bus_variant
    = std::variant<console4bus, console5bus, console6bus, console7bus>;

  using channels_variant = std::variant<
    std::array<console4channel, N>,
    std::array<console5channel, N>,
    std::array<console6channel, N>,
    std::array<console7channel, N>>;

  bus_variant      _bus;
  channels_variant _channels;

  int             _console_type = 0;
  plugin_context* plugcontext   = nullptr;
};

}} // namespace artv::airwindows
