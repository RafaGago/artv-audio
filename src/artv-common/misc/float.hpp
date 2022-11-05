#pragma once

#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/short_ints.hpp"

namespace artv {

template <class T, uint N_exp, uint N_frac>
struct float_layout {
  static_assert ((1 + N_exp + N_frac) == (sizeof (T) * 8), "No padding!");
  static_assert (N_exp > 0);
  static_assert (N_frac > 0);

  using value_type             = T;
  static constexpr uint n_exp  = N_exp;
  static constexpr uint n_frac = N_frac;

// remember, this codebase is restricted to clang  (and probably GCC)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  T frac : N_frac;
  T exp : N_exp;
  T sign : 1;
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  T sign : 1;
  T exp : N_exp;
  T frac : N_frac;
#else
#error "Endianess unsupported"
#endif
};
//------------------------------------------------------------------------------
template <uint N_exp = 5>
using f16_layout = float_layout<u16, N_exp, 15 - N_exp>;
using f32_layout = float_layout<u32, 8, 23>;
using f64_layout = float_layout<u64, 11, 52>;
//------------------------------------------------------------------------------
constexpr uint f16pack_dftz    = 1 << 0;
constexpr uint f16pack_nearest = 1 << 1;

// A homegrown storage format for float16, as __fp16 had linking problems and
// probably issues on Clang-cl at the time of writing this. At least this wheel
// reinvention allows to use the ranges better suited for audio and to remove
// unneeded features. It is probably slower than compiler generated code with
// hardware support.
//
// Assumes audio code: Denormals flushed to zero. No Inf/NaN. Inputs externally
// clamped to the allowed range.
//
// Used this as a reference.
// https://gist.github.com/rygorous/2156668
// https://fgiesen.wordpress.com/2012/03/28/half-to-float-done-quic/
//
// "N_exp": Number of exponent bits.
//
// "Max_exp:" is the maximum base 2 exponent the representation will reach.
//  Defaults to a symmetric range between positive and negative exponents. E.g
// exponent 0 = 2^0 = 1. Range [-1.99, 1.99].
//
// "flags":
// -"f16pack_nearest": Round to nearest instead of truncation. Notice that this
//  conversion doesn't support infinity, so enabling rounding can make values on
//  the extremes of the range to wrap.  This small headroom has to be considered
//  when/if externally clamping the range.
//
// -"f16pack_dftz": Flush denormalized float16 values to zero. Denormals allow
//  to win some extra low precision digits at the expense of CPU.
//
template <
  uint N_exp,
  int  Max_exp = (int) lsb_mask<uint> (N_exp) / 2 /*signed*/,
  uint flags   = f16pack_nearest>
class f16pack {
public:
  using half_layout                    = f16_layout<N_exp>;
  static constexpr int  max_exp        = Max_exp;
  static constexpr bool rounds_nearest = !!(flags & f16pack_nearest);
  static constexpr bool denormal_ftz   = !!(flags & f16pack_dftz);
  //----------------------------------------------------------------------------
  static constexpr u16 encode (float f)
  {
    f32_layout single;
    memcpy (&single, &f, sizeof f);

    half_layout half {};
    int         exp = (int) single.exp - rebias;
    assert (exp < exp_end);

    bool zero {};
    if constexpr (!denormal_ftz) {
      zero = std::abs (f) == 0.f;
    }
    if (exp > 0 || zero) {
      half.sign = single.sign;
      // zeroes on audio are very frequent, placed on the fast path.
      half.exp  = zero ? 0 : exp;
      half.frac = single.frac >> shift;
      u16 r {};
      memcpy (&r, &half, sizeof half);
      if constexpr (rounds_nearest) {
        uint round = !!(single.frac & bit<uint> (shift - 1));
        return r + round;
      }
      else {
        return r;
      }
    }
    else {
      if constexpr (denormal_ftz) {
        half.sign = single.sign;
        half.exp  = 0;
        half.frac = 0;
        u16 r {};
        memcpy (&r, &half, sizeof half);
        return r;
      }
      else {
        // store as denormal
        uint round = 0;
        half.sign  = single.sign;
        half.exp   = 0;
        half.frac  = 0;
        if ((-exp) <= half.n_frac) {
          uint frac = single.frac | bit<uint> (single.n_frac); // Hidden 1 bit
          uint fixed_shift = (shift + 1) - exp;
          half.frac        = frac >> fixed_shift;
          if constexpr (rounds_nearest) {
            round = !!(single.frac & bit<uint> (fixed_shift));
          }
        }
        u16 r {};
        memcpy (&r, &half, sizeof half);
        return r + round;
      }
    }
  }
  //----------------------------------------------------------------------------
  static constexpr float decode (u16 u)
  {
    f32_layout  single {};
    half_layout half {};

    memcpy (&half, &u, sizeof u);

    if (half.exp != 0 || u == 0) {
      // zeroes on audio are very frequent, placed on the fast path.
      single.frac = ((u32) half.frac) << shift;
      single.exp  = half.exp + ((u == 0) ? 0 : rebias);
      single.sign = half.sign;
      float r {};
      memcpy (&r, &single, sizeof single);
      return r;
    }
    else {
      if constexpr (denormal_ftz) {
        single.frac = 0;
        single.exp  = 0;
        single.sign = half.sign;

        float r {};
        memcpy (&r, &single, sizeof single);
        return r;
      }
      else {
#if 1 // using "last_bit_set" instruction to renormalize (set the virtual bit
      // number 24 as set).
        constexpr u32 hmask = lsb_mask<u32> (half.n_frac);
        uint          bit   = last_bit_set (half.frac);
        uint          exp   = half.n_frac - bit;
        single.frac         = (((u32) half.frac) << (exp + 1)) & hmask;
        single.exp          = rebias - exp;
        single.sign         = half.sign;

        float r {};
        memcpy (&r, &single, sizeof single);
        return r;
#else // using FPU renormalization
        constexpr u32 k_renorm = (rebias + 1) << single.n_frac;

        single.frac = ((u32) half.frac) << shift;
        single.exp  = half.exp + rebias + 1;

        float r {};
        memcpy (&r, &single, sizeof single);

        float renorm {};
        memcpy (&renorm, &k_renorm, sizeof renorm);
        r -= renorm; // fpu renormalization.

        memcpy (&single, &r, sizeof r);
        single.sign = half.sign;
        memcpy (&r, &single, sizeof single);
        return r;
#endif
      }
    }
  }
  //----------------------------------------------------------------------------
private:
  static constexpr uint bias_s = 127; // lsb_mask<uint> (f32_layout::n_exp - 1);
  static constexpr int  exp_end = (int) bit<uint> (half_layout::n_exp);
  static constexpr uint shift   = f32_layout::n_frac - half_layout::n_frac;
  static constexpr uint rebias  = bias_s - (exp_end - 1 - max_exp);

  static_assert (
    half_layout::n_exp <= f32_layout::n_exp,
    "N_exp has more bits than float's exponent");
  static_assert (max_exp <= 127, "Max_exp over float's range");
  static_assert (max_exp >= (-127 + (exp_end - 1)), "bottom range truncated");
};
//------------------------------------------------------------------------------
} // namespace artv
