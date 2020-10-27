#pragma once

#include <Downsampler2xSse.h>
#include <Upsampler2xSse.h>
#include <array>

#include "artv-common/dsp/types.hpp"

namespace artv {

template <class HiirObj, size_t NCoeffs>
class hiir_wrapper;
//------------------------------------------------------------------------------
template <size_t NCoeffs>
class hiir_wrapper<hiir::Upsampler2xSse, NCoeffs>
  : public hiir::Upsampler2xSse<NCoeffs> {
public:
  using value_type               = float;
  static constexpr uint channels = 1;
  //----------------------------------------------------------------------------
  void reset (float) { clear_buffers(); }
  //----------------------------------------------------------------------------
  std::array<std::array<value_type, N>, channels> process_sample (
    contiguous_range<const value_type> in)
  {
    std::array < std::array<value_type, channels> ret;
    Upsampler2xSse::process_sample (ret[0][0], ret[0][1], in[0]);
    return ret;
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <size_t NCoeffs>
class hiir_wrapper<hiir::Downsampler2xSse, NCoeffs>
  : public hiir::Downsampler2xSse<NCoeffs> {
public:
  using value_type               = float;
  static constexpr uint channels = 1;
  //----------------------------------------------------------------------------
  void reset (float) { clear_buffers(); }
  //----------------------------------------------------------------------------
  std::array<value_type, channels> process_sample (
    contiguous_range<contiguous_range<const value_type>> in)
  {
    std::array<value_type, channels> ret
      = {Downsampler2xSse::process_sample (in[0][0], in[0][1])};
    return ret;
  }
  //----------------------------------------------------------------------------
}; // namespace artv
//------------------------------------------------------------------------------
} // namespace artv
