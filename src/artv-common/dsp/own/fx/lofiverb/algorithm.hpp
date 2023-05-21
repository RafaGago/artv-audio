#pragma once

#include <array>
#include <type_traits>

#include "artv-common/dsp/own/fx/lofiverb/algorithm-engine.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "boost/mp11/algorithm.hpp"

namespace artv { namespace detail { namespace lofiverb {

// shared data structures and methods for algorithms, a bit free-form/chaotic...
class algorithm {
public:
  template <class T>
  struct smoothed_parameters {
    std::array<T, max_block_size>     decay;
    std::array<T, max_block_size>     character;
    std::array<T, max_block_size>     mod;
    std::array<float, max_block_size> stereo;
  };

  struct unsmoothed_parameters {
    float lf_amt;
    float hf_amt;
  };

  //----------------------------------------------------------------------------
  template <class T>
  using block_arr = std::array<T, max_block_size>;
  template <class T>
  using fb_block_arr = std::array<T, max_block_size + 1>;
  //----------------------------------------------------------------------------
  using fixpt_t = detail::lofiverb::fixpt_t;
  //----------------------------------------------------------------------------
  static constexpr auto blank     = detail::lofiverb::defaulted;
  static constexpr auto add_to    = detail::lofiverb::add_to;
  static constexpr auto overwrite = detail::lofiverb::overwrite;
  //----------------------------------------------------------------------------
  template <uint... Idxs>
  using sl = stage_list<Idxs...>;

protected:
  //----------------------------------------------------------------------------
  template <uint First, uint Last, class E, class T, class... Ts>
  void run_cascade (E& engine, xspan<T> io, Ts&&... args)
  {
    static_assert (Last >= First);
    mp11::mp_for_each<mp11::mp_iota_c<(Last + 1) - First>> ([&] (auto i) {
      engine.run (sl<First + i> {}, io, std::forward<Ts> (args)...);
    });
  }
  //----------------------------------------------------------------------------
  // just a convenience function for iterating block loops while not bloating
  // the code with more loop unroll hints than necessary
  template <class T, class F>
  static void span_visit (xspan<T> block, F&& visitor)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block.size(); ++i) {
      visitor (block[i], i);
    }
  }
  //----------------------------------------------------------------------------
  // just a convenience function for iterating block loops while not bloating
  // the code with more loop unroll hints than necessary
  template <class T, class U>
  static void span_mul (xspan<T> dst_lhs, U&& rhs)
  {
    span_visit (dst_lhs, [&] (auto& v, uint i) {
      if constexpr (is_array_subscriptable_v<U>) {
        v = (T) (v * rhs[i]);
      }
      else {
        v = (T) (v * rhs);
      }
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class U>
  static constexpr void span_add (xspan<T> dst, xspan<T const> lhs, U&& rhs)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < lhs.size(); ++i) {
      if constexpr (is_array_subscriptable_v<U>) {
        dst[i] = (T) (lhs[i] + rhs[i]);
      }
      else {
        dst[i] = (T) (lhs[i] + rhs);
      }
    }
  }

  template <class T, class U>
  static constexpr void span_add (xspan<T> dst_lhs, U&& rhs)
  {
    span_visit (dst_lhs, [&] (auto& v, uint i) {
      if constexpr (is_array_subscriptable_v<U>) {
        v = (T) (v + rhs[i]);
      }
      else {
        v = (T) (v + rhs);
      }
    });
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class V>
  static constexpr void span_add_with_factor (
    xspan<T>       dst,
    xspan<T const> lhs,
    U&&            rhs,
    V&&            factor)
  {
    constexpr bool u_subscriptable = is_array_subscriptable_v<U>;
    constexpr bool v_subscriptable = is_array_subscriptable_v<V>;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < lhs.size(); ++i) {
      if constexpr (u_subscriptable && v_subscriptable) {
        dst[i] = (T) (lhs[i] + rhs[i] * factor[i]);
      }
      else if constexpr (u_subscriptable && !v_subscriptable) {
        dst[i] = (T) (lhs[i] + rhs[i] * factor);
      }
      else if constexpr (!u_subscriptable && v_subscriptable) {
        dst[i] = (T) (lhs[i] + rhs * factor[i]);
      }
      else {
        dst[i] = (T) (lhs[i] + rhs * factor);
      }
    }
  }

  template <class T, class U, class V>
  static constexpr void span_add_with_factor (
    xspan<T> dst_lhs,
    U&&      rhs,
    V&&      factor)
  {
    span_add_with_factor (
      dst_lhs,
      dst_lhs.to_const(),
      std::forward<U> (rhs),
      std::forward<V> (factor));
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class V>
  static T clamp (T v, U min, V max)
  {
    if constexpr (std::is_same_v<T, float>) {
      return std::clamp (v, T {} + min, T {} + max);
    }
    else {
      return fixpt_clamp (v, T {} + min, T {} + max);
    }
  }
  //----------------------------------------------------------------------------
  // naive linear crossfade. 1 selects s2, 0 dst
  template <class T, class U, class V>
  static void crossfade (xspan<T> dst, U s2, V ctrl)
  {
    span_visit (dst, [&] (auto spl, uint i) {
      auto ctrl2 = (T) (1_r - ctrl[i]);
      dst[i]     = (T) (s2[i] * ctrl[i] + ctrl2 * spl);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  static T load_float (float v)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      return fixpt_t::from_float (v);
    }
    else {
      return v;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static float as_float (T v)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      return v.to_floatp();
    }
    else {
      return v;
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static auto tick_lfo (lfo<4>& lfo)
  {
    if constexpr (std::is_same_v<fixpt_t, T>) {
      auto ret = lfo.tick_sine_fixpt().spec_cast<fixpt_t>().value();
      // lfos return fixed points using vectors, we use scalars
      return vec_cast<fixpt_t::value_type> (ret);
    }
    else {
      return lfo.tick_sine();
    }
  }
  //----------------------------------------------------------------------------
  // reminder, unrolled to avoid unncessary quantization on fixed point
  template <class T>
  static void hadamard4 (std::array<T*, 4> io, uint block_size)
  {
    // these quantizations could use fraction saving too if moved to the
    // engine...
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block_size; ++i) {
      auto y1 = io[0][i] - io[1][i];
      auto y2 = io[0][i] + io[1][i];
      auto y3 = io[2][i] - io[3][i];
      auto y4 = io[2][i] + io[3][i];

      io[0][i] = (T) ((y1 - y3) * 0.5_r);
      io[1][i] = (T) ((y2 - y4) * 0.5_r);
      io[2][i] = (T) ((y1 + y3) * 0.5_r);
      io[3][i] = (T) ((y2 + y4) * 0.5_r);
    }
  }
  //----------------------------------------------------------------------------
  // reminder, unrolled to avoid unncessary quantization on fixed point
  template <class T>
  static void hadamard2 (std::array<T*, 2> io, uint block_size)
  {
    // these quantizations could use fraction saving too if moved to the
    // engine...
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block_size; ++i) {
      auto y1  = io[0][i] - io[1][i];
      auto y2  = io[0][i] + io[1][i];
      io[0][i] = (T) (y1 * 0.707106781_r);
      io[1][i] = (T) (y2 * 0.707106781_r);
    }
  }
  //----------------------------------------------------------------------------
  // 1 - (x - 1) * (x - 1) : a cheap curve that grows faster than linear,
  // similar to sqrt but compresses more the values closer to one instead of the
  // ones closer to zero
  template <class T>
  static auto fastgrowth (T x)
  {
    x = (T) (x - 1_r);
    return (T) (1_r - x * x);
  }
  //----------------------------------------------------------------------------
  // Equal power crossfade approx, 0 selects v1. 1 selects v2, results on v1
  template <class T, class F>
  static void ep_crossfade (
    xspan<T*>       v1, // in/out
    xspan<T const*> v2,
    T const*        cf,
    uint            block_size,
    F&&             cf_transform)
  {
    crossfade_impl (
      v1, v2, cf, block_size, std::forward<F> (cf_transform), 1.4186_r);
  }

  template <class T, class F>
  static void ep_crossfade (
    xspan<T*>       v1, // in/out
    xspan<T const*> v2,
    T const*        cf,
    uint            block_size)
  {
    crossfade_impl (
      v1, v2, cf, block_size, [] (auto x) { return x; }, 1.4186_r);
  }
  //----------------------------------------------------------------------------
  // Equal gain crossfade approx, 0 selects v1. 1 selects v2, results on v1
  template <class T, class F>
  static void eg_crossfade (
    xspan<T*>       v1, // in/out
    xspan<T const*> v2,
    T const*        cf,
    uint            block_size,
    F&&             cf_transform)
  {
    crossfade_impl (
      v1, v2, cf, block_size, std::forward<F> (cf_transform), 0.70912_r);
  }

  template <class T, class F>
  static void eg_crossfade (
    xspan<T*>       v1, // in/out
    xspan<T const*> v2,
    T const*        cf,
    uint            block_size)
  {
    crossfade_impl (
      v1, v2, cf, block_size, [] (auto x) { return x; }, 0.70912_r);
  }
  //----------------------------------------------------------------------------
protected:
  //----------------------------------------------------------------------------
  // https://signalsmith-audio.co.uk/writing/2021/cheap-energy-crossfade/
  // 0 selects v1. 1 selects v2
  template <class T, class R, class F>
  static void crossfade_impl (
    xspan<T*>       v1,
    xspan<T const*> v2,
    T const*        cf,
    uint            block_size,
    F&&             transform,
    R               constant)
  {
    using acum
      = std::conditional_t<std::is_floating_point_v<T>, T, fixpt_s<1, 4, 27>>;
    assert (v1.size() == v2.size());
    // not optimized, just eyebaled the values
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block_size; ++i) {
      acum x2 = transform (cf[i]);
      acum x1 = 1_r - x2;
      acum a  = x1 * x2;
      acum b  = a * (1_r + constant * a);
      acum c  = b + x1;
      acum d  = b + x2;
      c *= c;
      d *= d;
      for (uint j = 0; j < v1.size(); ++j) {
        v1[j][i] = (T) (v1[j][i] * c + v2[j][i] * d);
      }
    }
  }
  //----------------------------------------------------------------------------
};

}}} // namespace artv::detail::lofiverb
