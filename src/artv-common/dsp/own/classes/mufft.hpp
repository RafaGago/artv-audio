#pragma once

// Mufft is very competitive on float vectors.

#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#include <xmmintrin.h>

#include "mufft/fft.h"

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace mufft {
//------------------------------------------------------------------------------
template <class T>
class fft;
//------------------------------------------------------------------------------
template <>
class fft<float> {
public:
  using value_type                                = float;
  static constexpr uint io_alignment              = 32;
  static constexpr bool io_can_alias              = false;
  static constexpr bool reorder_requires_buffer   = false;
  static constexpr bool transform_requires_buffer = false;
  //----------------------------------------------------------------------------
  fft() = default;
  ~fft() { reset (0); }
  //----------------------------------------------------------------------------
  fft (const fft&) = delete;
  fft& operator= (const fft&) = delete;
  //----------------------------------------------------------------------------
  fft (fft&& other) { *this = std::forward<fft<value_type>> (other); };
  //----------------------------------------------------------------------------
  fft& operator= (fft&& other)
  {
    if (this == &other) {
      return *this;
    }
    _fwd         = other._fwd;
    other._fwd   = nullptr;
    _bckwd       = other._bckwd;
    other._bckwd = nullptr;
    return *this;
  };
  //----------------------------------------------------------------------------
  bool reset (uint blocksize, bool complex = true)
  {
    if (blocksize != 0) {
      if (_fwd || _bckwd) {
        reset (0);
      }
      assert (
        is_pow2 (blocksize)
        && blocksize >= (io_alignment / sizeof (value_type)));
      if (complex) {
        _fwd   = mufft_create_plan_1d_c2c (blocksize, MUFFT_FORWARD, 0);
        _bckwd = mufft_create_plan_1d_c2c (blocksize, MUFFT_INVERSE, 0);
      }
      else {
        _fwd   = mufft_create_plan_1d_r2c (blocksize, 0);
        _bckwd = mufft_create_plan_1d_c2r (blocksize, 0);
      }
      if (unlikely (!_fwd || !_bckwd)) {
        reset (0);
        return false;
      }
    }
    else {
      if (_fwd) {
        mufft_free_plan_1d (_fwd);
        _fwd = nullptr;
      }
      if (_bckwd) {
        mufft_free_plan_1d (_bckwd);
        _bckwd = nullptr;
      }
    }
    return true;
  }
  //----------------------------------------------------------------------------
  void forward (crange<value_type> out, const crange<value_type> in)
  {
    align_assert (out.data());
    align_assert (in.data());
    assert (out.data() != in.data()); // mufft can't alias
    mufft_execute_plan_1d (_fwd, out.data(), in.data());
  }
  //----------------------------------------------------------------------------
  void backward (crange<value_type> out, const crange<value_type> in)
  {
    align_assert (out.data());
    align_assert (in.data());
    assert (out.data() != in.data()); // mufft can't alias
    mufft_execute_plan_1d (_bckwd, out.data(), in.data());
  }
  //----------------------------------------------------------------------------
  void reorder_after_forward (crange<value_type>) {}
  //----------------------------------------------------------------------------
  void reorder_before_backward (crange<value_type>) {}
  //----------------------------------------------------------------------------
  void forward_ordered (crange<value_type> out, const crange<value_type> in)
  {
    forward (out, in);
  }
  //----------------------------------------------------------------------------
  void backward_ordered (crange<value_type> out, const crange<value_type> in)
  {
    backward (out, in);
  }
  //----------------------------------------------------------------------------
  void data_rescale (crange<value_type> v, uint fft_size, bool complex)
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
  mufft_plan_1d* _fwd   = nullptr;
  mufft_plan_1d* _bckwd = nullptr;
};

}} // namespace artv::mufft
