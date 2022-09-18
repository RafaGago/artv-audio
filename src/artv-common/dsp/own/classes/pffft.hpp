#pragma once

// Using pffft for doubles

#include <cassert>
#include <utility>

#include "pffft/pffft_double.h"

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace pffft {
//------------------------------------------------------------------------------
template <class T>
class fft;
//------------------------------------------------------------------------------
template <>
class fft<double> {
public:
  using value_type                                = double;
  static constexpr uint io_alignment              = 32;
  static constexpr bool io_can_alias              = true;
  static constexpr bool reorder_requires_buffer   = true;
  static constexpr bool transform_requires_buffer = true;
  //----------------------------------------------------------------------------
  fft() = default;
  ~fft() { reset (0); }
  //----------------------------------------------------------------------------
  fft (const fft&)            = delete;
  fft& operator= (const fft&) = delete;
  //----------------------------------------------------------------------------
  fft (fft&& other) { *this = std::forward<fft<value_type>> (other); };
  //----------------------------------------------------------------------------
  fft& operator= (fft&& other)
  {
    if (this == &other) {
      return *this;
    }
    _setup       = other._setup;
    other._setup = nullptr;
    return *this;
  };
  //----------------------------------------------------------------------------
  bool reset (uint blocksize, bool complex = true)
  {
    auto pfft_type = complex ? PFFFT_COMPLEX : PFFFT_REAL;

    if (blocksize != 0) {
      if (_setup) {
        reset (0);
      }
      assert (pffftd_is_valid_size (blocksize, pfft_type));
      _setup = pffftd_new_setup (blocksize, pfft_type);

      if (unlikely (!_setup)) {
        // probably a bad size.
        reset (0);
        return false;
      }
    }
    else {
      if (_setup) {
        pffftd_destroy_setup (_setup);
        _setup = nullptr;
      }
    }
    return true;
  }
  //----------------------------------------------------------------------------
  // Notice PFFT forward and backwards doesn't even leave the bins as an
  // unordered complex array, so in practice forward/backward is for little use
  // except for performing some functionality not yet wrapped.
  void forward (
    xspan<value_type>       out,
    const xspan<value_type> in,
    xspan<value_type>       work)
  {
    align_assert (out.data());
    align_assert (in.data());
    align_assert (work.data());
    pffftd_transform (_setup, &in[0], &out[0], &work[0], PFFFT_FORWARD);
  }
  //----------------------------------------------------------------------------
  // Notice PFFT forward and backwards doesn't even leave the bins as an
  // unordered complex array, so in practice forward/backward is for little use
  // except for performing some functionality not yet wrapped.
  void backward (
    xspan<value_type>       out,
    const xspan<value_type> in,
    xspan<value_type>       work)
  {
    align_assert (out.data());
    align_assert (in.data());
    align_assert (work.data());
    pffftd_transform (_setup, &in[0], &out[0], &work[0], PFFFT_BACKWARD);
  }
  //----------------------------------------------------------------------------
  void reorder_after_forward (xspan<value_type> out, xspan<value_type> work)
  {
    align_assert (out.data());
    align_assert (work.data());
    xspan_copy (work, out);
    pffftd_zreorder (_setup, &work[0], &out[0], PFFFT_FORWARD);
  }
  //----------------------------------------------------------------------------
  void reorder_before_backward (xspan<value_type> out, xspan<value_type> work)
  {
    align_assert (out.data());
    align_assert (work.data());
    xspan_copy (work, out);
    pffftd_zreorder (_setup, &work[0], &out[0], PFFFT_BACKWARD);
  }
  //----------------------------------------------------------------------------
  void forward_ordered (
    xspan<value_type>       out,
    const xspan<value_type> in,
    xspan<value_type>       work)
  {
    align_assert (out.data());
    align_assert (in.data());
    align_assert (work.data());
    pffftd_transform_ordered (_setup, &in[0], &out[0], &work[0], PFFFT_FORWARD);
  }
  //----------------------------------------------------------------------------
  void backward_ordered (
    xspan<value_type>       out,
    const xspan<value_type> in,
    xspan<value_type>       work)
  {
    align_assert (out.data());
    align_assert (in.data());
    align_assert (work.data());
    pffftd_transform_ordered (
      _setup, &in[0], &out[0], &work[0], PFFFT_BACKWARD);
  }
  //----------------------------------------------------------------------------
  void data_rescale (xspan<value_type> v, uint fft_size, bool complex)
  {
    uint elems = fft_size * (complex ? 2 : 1);
    assert (v.size() >= elems);
    auto f = (value_type) 1 / (value_type) fft_size;
    for (uint i = 0; i < elems; ++i) {
      v[i] *= f;
    }
  }
  //----------------------------------------------------------------------------
private:
  static void align_assert (void const* ptr)
  {
    assert ((((uintptr_t) ptr) & (io_alignment - 1)) == 0);
  }
  //----------------------------------------------------------------------------
  PFFFTD_Setup* _setup = nullptr;
};
//------------------------------------------------------------------------------
}} // namespace artv::pffft
