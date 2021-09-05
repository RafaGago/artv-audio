#pragma once

#include <array>
//#include <variant>
//#include "artv-common/misc/mp11.hpp"

namespace artv {
//------------------------------------------------------------------------------
class plugin_context;
//------------------------------------------------------------------------------
template <class T>
struct stereo_processor {
  virtual ~stereo_processor() {}
  virtual void reset (plugin_context& pc) = 0;
  virtual void process (crange<T*> outs, crange<T const*> ins, uint samples)
    = 0;

  // TODO: All of this is implementable on "stereo_processor_adapt", just not
  // done yet.
  //
  // virtual uint get_parameter_count () = 0;
  // virtual void set_parameter (unsigned id, std::variant<int, float> value) =
  // 0;
  // using param_desc = std::variant<float_parameter, int_parameter>;
  // virtual param_desc get_parameter_traits (unsigned id) = 0;
  // virtual dsp_type get_dsp_type() = 0;
};
//------------------------------------------------------------------------------
// to go back from template based to virtual, in case there is a need to reduce
// the code footprint in some stages.
template <class Impl, class T = float>
struct stereo_processor_adapt : public Impl, public stereo_processor<T> {
public:
  using stereo_processor_type = Impl;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc) override { Impl::reset (pc); };
  //----------------------------------------------------------------------------
  void process (crange<T*> outs, crange<T const*> ins, uint samples) override
  {
    Impl::process (outs, ins, samples);
  }
  //----------------------------------------------------------------------------
private:
  using Impl::process;
  using Impl::reset;
};
//------------------------------------------------------------------------------
} // namespace artv
