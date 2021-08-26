#pragma once

// just a wrapper to be able to switch FFT implementations if needed...

// FFTS seems to be broken for power of 2 sizes...

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
  using value_type                   = float;
  static constexpr uint io_alignment = 32;

  fft() = default;
  ~fft() { reset (0); }

  fft (const fft&) = delete;
  fft& operator= (const fft&) = delete;

  fft (fft&& other) { *this = std::forward<fft<value_type>> (other); };
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

  bool reset (uint blocksize)
  {
    assert (
      is_pow2 (blocksize) && blocksize >= (io_alignment / sizeof (value_type))
      || blocksize == 0);

    if (blocksize != 0) {
      _fwd   = mufft_create_plan_1d_c2c (blocksize, MUFFT_FORWARD, 0);
      _bckwd = mufft_create_plan_1d_c2c (blocksize, MUFFT_INVERSE, 0);

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

  void forward_transform (value_type* out, const value_type* in)
  {
    align_assert (out);
    align_assert (in);
    mufft_execute_plan_1d (_fwd, out, in);
  }

  void backward_transform (value_type* out, const value_type* in)
  {
    align_assert (out);
    align_assert (in);
    mufft_execute_plan_1d (_bckwd, out, in);
  }

private:
  static void align_assert (void const* ptr)
  {
    assert ((((uintptr_t) ptr) & (io_alignment - 1)) == 0);
  }

  mufft_plan_1d* _fwd   = nullptr;
  mufft_plan_1d* _bckwd = nullptr;
};
//------------------------------------------------------------------------------
template <class T, size_t... blocksizes>
class initialized_ffts;

template <size_t... blocksizes>
class initialized_ffts<float, blocksizes...> {
public:
  using value_type                   = float;
  static constexpr uint io_alignment = fft<float>::io_alignment;

  initialized_ffts()
  {
    mp_foreach_idx (blocksize_list {}, [=] (auto index, auto bsz) {
      _ffts[index.value].reset (bsz.value);
    });
    _tmpbuff.resize (
      div_ceil<uint> (biggest_blocksize(), simdwrapper::n_elems) * 2);
  }

  void forward_transform (value_type* out, value_type* in, uint blocksize)
  {
    int idx = get_fft_idx (blocksize);
    assert (idx >= 0);
    _ffts[idx].forward_transform (out, in);
  }

  void forward_transform (value_type* io, uint blocksize)
  {
    int idx = get_fft_idx (blocksize);
    assert (idx >= 0);
    memcpy (_tmpbuff.data(), io, blocksize * 2 * sizeof (value_type));
    _ffts[idx].forward_transform (io, &_tmpbuff[0][0]);
  }

  void forward_permute (value_type* buffer, uint blocksize) {}

  void backward_transform (
    value_type*       out,
    const value_type* in,
    uint              blocksize)
  {
    int idx = get_fft_idx (blocksize);
    assert (idx >= 0);
    _ffts[idx].backward_transform (out, in);
  }

  void backward_transform (value_type* io, uint blocksize)
  {
    int idx = get_fft_idx (blocksize);
    assert (idx >= 0);
    memcpy (_tmpbuff.data(), io, blocksize * 2 * sizeof (value_type));
    _ffts[idx].backward_transform (io, &_tmpbuff[0][0]);
  }

  void backward_permute (value_type* buffer, uint blocksize) {}

private:
  static int get_fft_idx (uint blocksize)
  {
    int ret = -1;
    mp_foreach_idx (blocksize_list {}, [=, &ret] (auto index, auto bsz) {
      if (bsz.value == blocksize) {
        ret = index.value;
      }
    });
    return ret;
  }

  static uint biggest_blocksize()
  {
    uint maxbsz = 0;
    mp11::mp_for_each<blocksize_list> ([&] (auto bsz) {
      maxbsz = (bsz.value > maxbsz) ? bsz.value : maxbsz;
    });
    return maxbsz;
  }

  using simdwrapper
    = simd_mem<float, io_alignment / sizeof (float), io_alignment>;

  using blocksize_list = mp11::mp_list_c<int, blocksizes...>;

  std::array<fft<float>, sizeof...(blocksizes)> _ffts;
  std::vector<simdwrapper>                      _tmpbuff;
};
//------------------------------------------------------------------------------
}} // namespace artv::mufft
