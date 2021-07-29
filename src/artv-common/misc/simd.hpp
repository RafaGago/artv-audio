#pragma once

#include <array>
#include <utility>

#include <type_traits>

#include <xsimd/xsimd.hpp>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
//------------------------------------------------------------------------------
// This seems a bit redundant, but it was previously wrapping the JUCE's SIMD
// wrappers.

#if 0
template <class T>
class simd_register : public ::xsimd::simd_type<T> {
private:
  using parent = ::xsimd::simd_type<T>;

public:
  struct aligned {};

  static_assert (std::is_arithmetic_v<T>, "");

  using value_type                   = T;
  static constexpr size_t n_builtins = parent::size;
  static constexpr size_t alignment  = n_builtins * sizeof (value_type);
  //----------------------------------------------------------------------------
  template <size_t N>
  using array_rounded = std::array<T, round_ceil<size_t> (N, n_builtins)>;
  //----------------------------------------------------------------------------
  template <class... Ts>
  simd_register (Ts&&... args) : parent (std::forward<Ts> (args)...)
  {}
  //----------------------------------------------------------------------------
  simd_register (contiguous_range<value_type> r) { load (r); }
  simd_register (contiguous_range<value_type> r, aligned) { aligned_load (r); }
  //----------------------------------------------------------------------------
  simd_register (contiguous_range<const value_type> r) { load (r); }
  simd_register (contiguous_range<const value_type> r, aligned)
  {
    aligned_load (r);
  }
  //----------------------------------------------------------------------------
  simd_register (std::array<value_type, n_builtins> r) { load (r); }
  simd_register (std::array<value_type, n_builtins> r, aligned)
  {
    aligned_load (r);
  }
  //----------------------------------------------------------------------------
  template <int N, std::enable_if_t<N >= (n_builtins + 1)>* = nullptr>
  simd_register (std::array<const value_type, N> r)
  {
    load (r);
  }
  //----------------------------------------------------------------------------
  void load (contiguous_range<const value_type> src)
  {
    assert (src.size() >= n_builtins);
    for (uint i = 0; i < n_builtins; ++i) {
      this->operator[] (i) = src[i];
    }
  }
  //----------------------------------------------------------------------------
  void store (contiguous_range<value_type> dst)
  {
    assert (dst.size() >= n_builtins);
    for (uint i = 0; i < n_builtins; ++i) {
      dst[i] = this->operator[] (i);
    }
  }
  //----------------------------------------------------------------------------
  void aligned_load (contiguous_range<const value_type> src)
  {
    assert (src.size() >= n_builtins);
    *this = this->fromRawArray (src.data());
  }
  //----------------------------------------------------------------------------
  void aligned_store (contiguous_range<value_type> dst)
  {
    assert (dst.size() >= n_builtins);
    this->copyToRawArray (dst.data());
  }
  //----------------------------------------------------------------------------
  simd_register& operator= (contiguous_range<value_type> const r)
  {
    load (r);
    return *this;
  }
  //----------------------------------------------------------------------------
  using parent::operator[];
  using parent::operator+=;
  using parent::operator-=;
  using parent::operator*=;
  using parent::operator/=;
  using parent::operator=;
  using parent::operator&=;
  using parent::operator|=;
  using parent::operator^=;
  using parent::operator+;
  using parent::operator-;
  using parent::operator*;
  using parent::operator/;
  using parent::operator&;
  using parent::operator|;
  using parent::operator^;
  using parent::operator~;
  using parent::operator==;
  using parent::operator!=;
};
//------------------------------------------------------------------------------
#endif

template <class T, size_t size>
using simd_batch = xsimd::batch<T, size>;

template <class T, size_t instr_set_bytes>
using simd_reg = simd_batch<T, instr_set_bytes / sizeof (T)>;

static constexpr uint sse_bytes  = 16;
static constexpr uint avx_bytes  = 32;
static constexpr uint avx2_bytes = 64;

// This could come from flags at some point on the future
static constexpr uint default_simd_align = sse_bytes;

// SSE or equivalent
using simd_flt = simd_reg<float, 16>;
using simd_dbl = simd_reg<double, 16>;
using simd_u32 = simd_reg<u32, 16>;
using simd_s32 = simd_reg<s32, 16>;
using simd_u64 = simd_reg<u64, 16>;
using simd_s64 = simd_reg<s64, 16>;

// AVX256 or equivalent.
using simd_flt_x2 = simd_reg<float, 32>;
using simd_dbl_x2 = simd_reg<double, 32>;
using simd_u32_x2 = simd_reg<u32, 32>;
using simd_s32_x2 = simd_reg<s32, 32>;
using simd_u64_x2 = simd_reg<u64, 32>;

//------------------------------------------------------------------------------
// just an array rounded to a simd width.
template <class T, size_t N, size_t instr_set_bytes>
using simd_array
  = std::array<T, round_ceil<size_t> (N, (instr_set_bytes / sizeof (T)))>;
//------------------------------------------------------------------------------
template <class T, size_t N>
static simd_batch<T, N> sgn_no_zero (
  simd_batch<T, N> x,
  simd_batch<T, N> neg = simd_batch<T, N> {(T) -1.},
  simd_batch<T, N> pos = simd_batch<T, N> {(T) 1.})
{
  using batch = simd_batch<T, N>;
  return xsimd::select (x < batch {(T) 0.}, neg, pos);
}
//------------------------------------------------------------------------------
#define XSIMD_BROKEN_W_FAST_MATH 1
// https://github.com/xtensor-stack/xsimd/issues/515
// at least "log, pow, exp, tanh" are broken with ffast-math enabled

} // namespace artv
