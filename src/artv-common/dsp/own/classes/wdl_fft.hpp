#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>

#include "WDL/fft.h"

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace wdl {
//------------------------------------------------------------------------------
template <class T>
class initialized_ffts;

// thread safety: None
template <>
class initialized_ffts<double> {
public:
  static_assert (std::is_same<double, WDL_FFT_REAL>::value, "");

  using value_type                    = double;
  static constexpr uint min_blocksize = 16;
  static constexpr uint max_blocksize = 32768;

  initialized_ffts()
  {
    // TODO: maybe add a static mutex here...
    // can be called many times. Single threaded though...
    WDL_fft_init();
    _work_buffer.resize (max_blocksize * 2); // max fft
  }

  void forward_transform (value_type* out, const value_type* in, uint blocksize)
  {
    bsz_assert (blocksize);
    memcpy (out, in, blocksize * sizeof *out * 2 * (out != in));
    WDL_fft (((WDL_FFT_COMPLEX*) ((void*) out)), blocksize, false);
  }

  void forward_transform (value_type* io, uint blocksize)
  {
    bsz_assert (blocksize);
    WDL_fft (((WDL_FFT_COMPLEX*) ((void*) io)), blocksize, false);
  }

  void forward_permute (value_type* buffer, uint blocksize)
  {
    bsz_assert (blocksize);
    constexpr uint n_complex = sizeof (WDL_FFT_COMPLEX) / sizeof (value_type);

    memcpy (_work_buffer.data(), buffer, blocksize * sizeof (*buffer) * 2);
    int* tbl = WDL_fft_permute_tab (blocksize);

    for (uint i = 0; i < blocksize; ++i) {
      memcpy (
        &buffer[i * n_complex],
        &_work_buffer[tbl[i] * n_complex],
        sizeof *buffer * n_complex);
    }
  }

  void backward_transform (
    value_type*       out,
    const value_type* in,
    uint              blocksize)
  {
    bsz_assert (blocksize);
    memcpy (out, in, blocksize * sizeof *out * 2 * (out != in));
    WDL_fft (((WDL_FFT_COMPLEX*) ((void*) out)), blocksize, true);
  }

  void backward_transform (value_type* io, uint blocksize)
  {
    bsz_assert (blocksize);
    WDL_fft (((WDL_FFT_COMPLEX*) ((void*) io)), blocksize, true);
  }

  void backward_permute (value_type* buffer, uint blocksize)
  {
    bsz_assert (blocksize);
    constexpr uint n_complex = sizeof (WDL_FFT_COMPLEX) / sizeof (value_type);

    memcpy (_work_buffer.data(), buffer, blocksize * (sizeof *buffer) * 2);
    int* tbl = WDL_fft_permute_tab (blocksize);

    for (uint i = 0; i < blocksize; ++i) {
      memcpy (
        &buffer[tbl[i] * n_complex],
        &_work_buffer[i * n_complex],
        sizeof *buffer * n_complex);
    }
  }

private:
  static void bsz_assert (uint blocksize)
  {
    assert (
      is_pow2 (blocksize) && blocksize >= min_blocksize
      && blocksize <= max_blocksize);
  }
  std::vector<value_type> _work_buffer;
};
//------------------------------------------------------------------------------
}} // namespace artv::wdl
