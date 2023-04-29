#pragma once

#include <type_traits>

#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/primes_table.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------

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

  // fixpt_t truncates when dropping fractional bits (leaky).
  // fixpt_tr rounds to the nearest (never reaches full zero).
  //
  // Using fixpt_t as default, with fixpt_tr at some points compensate the
  // truncating leakage.
  //----------------------------------------------------------------------------
  template <class T>
  using block_arr = std::array<T, max_block_size>;
  template <class T>
  using fb_block_arr = std::array<T, max_block_size + 1>;
  //----------------------------------------------------------------------------
  using fixpt_t  = detail::lofiverb::fixpt_t;
  using fixpt_tr = detail::lofiverb::fixpt_tr;
  using float16  = detail::lofiverb::float16;
  //----------------------------------------------------------------------------
  static constexpr auto blank     = detail::lofiverb::defaulted;
  static constexpr auto add_to    = detail::lofiverb::add_to;
  static constexpr auto overwrite = detail::lofiverb::overwrite;
  //----------------------------------------------------------------------------
  template <uint... Idxs>
  using sl = stage_list<Idxs...>;

protected:
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
  static void span_mul (xspan<T> block, U val)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block.size(); ++i) {
      block[i] = (T) (block[i] * val);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void span_add (
    xspan<T>       dst,
    xspan<T const> lhs,
    xspan<T const> rhs)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < rhs.size(); ++i) {
      dst[i] = (T) (lhs[i] + rhs[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void span_add (xspan<T> dst_lhs, xspan<T const> rhs)
  {
    span_add<T> (dst_lhs, dst_lhs, rhs);
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
  template <class T, class U, class V, class W>
  static void crossfade (xspan<T> dst, U s2, V ctrl, W one)
  {
    span_visit (dst, [&] (auto spl, uint i) {
      auto ctrl2 = (T) (one - ctrl[i]);
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
      // At this point "ret" has fixed traits, but different fixed point
      // conversions configured
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
  template <class T, class U>
  static auto fastgrowth (T x, U one)
  {
    x = (T) (x - one);
    return (T) (one - x * x);
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
//------------------------------------------------------------------------------
struct debug : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<debug, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f1 = 0.43f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.75f, 1.f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // left channel output conditioning part 1
      make_ap (71), // 0
      make_ap (116), // 1
      make_ap (172), // 2

      // right channel output conditioning part 2
      make_ap (72), // 3
      make_ap (117), // 4
      make_ap (175), // 5

      // left channel output conditioning part 2 (chorus)
      make_ap (1987, 0, 123), // 6
      make_lp(), // 7
      // right channel output conditioning part 2 (chorus)
      make_ap (2177, 0, 138), // 8
      make_lp(), // 9

      make_ap (277), // 10 // L last AP
      make_ap (274), // 11 // R last AP

      // m conditioning-diffusors
      make_lp (0.5f), // 12
      make_hp (0.87f), // 13
      make_ap (23, 0.8), // 14
      make_ap (130, 0.8), // 15
      make_ap (217, 0.8), // 16

      // s conditioning-diffusors
      make_lp (0.5f), // 17
      make_hp (0.87f), // 18
      make_ap (22, 0.8), // 19
      make_ap (131, 0.8), // 20
      make_ap (219, 0.8), // 21

      make_quantizer(), // 22
      make_quantizer(), // 23
      make_quantizer(), // 24
      make_quantizer(), // 25

      // block a iteration 1
      make_ap (153, 0.2), // 26
      make_ap (89, -0.04), // 27
      make_ap (60, 0.04), // 28
      make_delay (201), // 29

      // block b iteration 1
      make_delay (185), // 39

      // block c iteration 1
      make_ap (149, 0.2), // 31
      make_ap (83, 0.04), // 32
      make_ap (59, -0.04), // 33
      make_delay (225), // 34

      // block d iteration 1
      make_ap (167, -0.2), // 35
      make_ap (97, -0.04), // 36
      make_ap (67, 0.04), // 37
      make_delay (221), // 38

      make_quantizer(), // 39
      make_quantizer(), // 40
      make_quantizer(), // 41
      make_quantizer(), // 42

      // crossovers
      make_crossover2(), // 43
      make_crossover2(), // 44
      make_crossover2(), // 45
      make_crossover2(), // 46

      // block a iteration 2
      make_ap (119, 0.6), // 47
      make_ap (67, 0.6), // 48
      make_ap (47, -0.6), // 49
      make_block_delay (171), // feedback point // 50

      // block b iteration 2
      make_ap (114, -0.1), // nested 3x // 51
      make_ap (66, -0.04), // 52
      make_ap (47, -0.04), // 53
      make_ap (9), // 54
      make_block_delay (185), // feedback point // 55

      // block c iteration 2
      make_ap (116, -0.1), // nested 3x // 56
      make_ap (65, -0.04), // 57
      make_ap (46, -0.04), // 58
      make_ap (3), //  59
      make_block_delay (179), // feedback point // 60

      // block d iteration 2
      make_ap (121, -0.1), // nested 3x // 61
      make_ap (69, 0.04), // 62
      make_ap (47, 0.04), // 63
      make_block_delay (175) // feedback point // 64
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr_fb a_mem, d_mem;
    arr    b_mem, c_mem, tmp1, tmp2, tmp3;
    arr    lfo1, lfo2;

    xspan a {a_mem.data(), io.size() + 1};
    xspan b {b_mem.data(), io.size()};
    xspan c {c_mem.data(), io.size()};
    xspan d {d_mem.data(), io.size() + 1};

    // fetch a block from the output (a and d)
    rev.fetch_block (sl<50> {}, a, 1); // feedback, fetching block + 1 samples
    rev.fetch_block (sl<64> {}, d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // invert character first
      par.character[i] = (T) (rev.one - par.character[i]);
      tmp1[i]          = (T) (0.4_r + 0.3_r * par.character[i]);
      tmp2[i]          = (T) (-0.6_r * par.character[i]);
      tmp3[i]          = (T) (0.6_r * par.character[i]);
      // decay fixup
      auto d       = rev.one - par.decay[i];
      par.decay[i] = (T) (0.985_r - d * d * 0.21_r);
      // lfo
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = T {lfo[1]};
      // lfo3 = -lfo1, lfo4 = -lfo2
    }

    // output preconditioning
    rev.run (sl<0> {}, b, blank, tmp1);
    rev.run (sl<1> {}, b, blank, tmp2);
    rev.run (sl<2> {}, b, blank, tmp2);

    rev.run (sl<3> {}, c, blank, tmp1);
    rev.run (sl<4> {}, c, blank, tmp2);
    rev.run (sl<5> {}, c, blank, tmp2);

    // tmp1 and tmp2 free
    xspan cho1 {tmp1.data(), io.size()};
    xspan cho2 {tmp2.data(), io.size()};

    rev.run (sl<6> {}, cho1, b.to_const(), lfo1, tmp3);
    rev.run (sl<7> {}, cho1);
    rev.run (
      sl<8> {}, cho2, c.to_const(), [&] (uint i) { return -lfo1[i]; }, tmp3);
    rev.run (sl<9> {}, cho2);

    auto v1 = make_array (&b[0], &c[0]);
    auto v2 = make_array<T const*> (&cho1[0], &cho2[0]);
    ep_crossfade<T> (v1, v2, &par.mod[0], io.size(), [=] (auto v) {
      return fastgrowth ((T) (v * 0.5_r), rev.one);
    });

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      b[i] = (T) (b[i] + cho1[i]);
      c[i] = (T) (c[i] + cho2[i]);
    }

    rev.run (sl<10> {}, b, blank, tmp3);
    rev.run (sl<11> {}, c, blank, tmp3); // tmp3 free

    xspan i1 {tmp1.data(), io.size()};
    xspan i2 {tmp2.data(), io.size()};
    // output write + preparations
    span_visit (io, [&] (auto& spls, uint i) {
      i1[i]   = (T) ((spls[0] + spls[1]) * 0.25_r); // m
      i2[i]   = (T) ((spls[0] - spls[1]) * 0.25_r); // s
      spls[0] = b[i];
      spls[1] = c[i];
      b[i]    = (T) (0.8_r + 0.1_r * par.character[i]); // character 1 on b
      c[i]    = (T) (0.8_r * par.character[i]); // character 2 on c
    });

    // input preconditioning
    rev.run (sl<12> {}, i1);
    rev.run (sl<13> {}, i1);
    rev.run (sl<14> {}, i1, blank, b);
    rev.run (sl<15> {}, i1, blank, c);
    rev.run (sl<16> {}, i1, blank, c);

    rev.run (sl<17> {}, i2);
    rev.run (sl<18> {}, i2);
    rev.run (sl<19> {}, i2, blank, b);
    rev.run (sl<20> {}, i2, blank, c);
    rev.run (sl<21> {}, i2, blank, c);

    // fetch b and d feedback (a and d are already ready) 1 block only
    rev.fetch_block (sl<55> {}, b, 1); // feedback
    rev.fetch_block (sl<60> {}, c, 1); // feedback

    // quantized decay
    rev.run (sl<22> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<23> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<24> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<25> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // a single hadamard iteration can duplicate the range on one
      // of the channels, so halving once more (TODO needed?).
      b[i] = (T) (b[i] + tmp1[i] * 0.5_r);
      d[i] = (T) (d[i] + tmp2[i] * 0.5_r);
    }
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // channel a block 1
    rev.run (sl<26, 27, 28> {}, a);
    rev.run (sl<29> {}, a);

    // channel b block 1
    rev.run (sl<30> {}, b);

    // channel c block 1
    rev.run (sl<31, 32, 33> {}, c);
    rev.run (sl<34> {}, c);

    // channel d block 1
    rev.run (sl<35, 36, 37> {}, d);
    rev.run (sl<38> {}, d);

    // quantized decay
    rev.run (sl<39> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<40> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<41> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<42> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // block 2
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // swaps.
    a = xspan {d_mem.data(), io.size()};
    b = xspan {c_mem.data(), io.size()};
    c = xspan {a_mem.data(), io.size()};
    d = xspan {b_mem.data(), io.size()};

    // damp -----------------------------------
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.9f + upar.lf_amt * 0.09f);
    T fhi = load_float<T> (0.61f - upar.hf_amt * upar.hf_amt * 0.12f);
    T ghi = load_float<T> (0.6f + upar.hf_amt * 0.18f);

    rev.run (sl<43> {}, a, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<44> {}, b, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<45> {}, c, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<46> {}, d, flo, glo, fhi, rev.one, ghi);

    // channel a block 2
    rev.run (sl<47> {}, a);
    rev.run (sl<48> {}, a);
    rev.run (sl<49> {}, a);
    rev.push (sl<50> {}, a.to_const());

    // channel b block 2
    rev.run (sl<51, 52, 53> {}, b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      T m     = fastgrowth (par.mod[i], rev.one);
      tmp1[i] = (T) (0.4_r + 0.4_r * lfo2[i] * m);
      tmp2[i] = (T) (0.6_r + 0.4_r * -lfo2[i] * m); // for block c
    }
    rev.run (sl<54> {}, b, tmp1);
    rev.push (sl<55> {}, b.to_const());

    // channel c block 2
    rev.run (sl<56, 57, 58> {}, c);
    rev.run (sl<59> {}, c, tmp2);
    rev.push (sl<60> {}, c.to_const());

    // channel d block 2
    rev.run (sl<61, 62, 63> {}, d);
    rev.push (sl<64> {}, d.to_const());
  }
};
//------------------------------------------------------------------------------
class abyss : public algorithm {
public:
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<abyss, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f_er   = 0.3f + mod * 0.3f;
    auto f_late = 0.1f + mod * 1.2f;
    lfo.set_freq (f32_x4 {f_er, f_er, f_late, f_late}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp (0.3), // 0
      // diffusors
      make_ap (147, -0.707), // 1
      make_ap (183, 0.707), // 2
      make_ap (389, -0.6), // 3
      make_ap (401, 0.6), // 4
      // er (not ER at the end, as this is now another reverb...)
      make_ap (1367, 0.35, 71 + 70), // 5
      make_crossover2(), // 6
      make_quantizer(), // 7
      make_ap (1787, 0.5, 261), // 8
      make_block_delay (max_block_size), // 9 to allow block processing
      // loop1
      make_ap (977, 0.5 /*overridden*/, 51), // 10
      make_delay (2819), // 11
      make_crossover2(), // 12
      make_quantizer(), // 13
      make_ap (863, -0.5 /*overridden*/), // 14
      make_delay (1021), // 15
      make_quantizer(), // 16
      make_ap (1453, 0.618), // 17
      make_block_delay (
        787), // 18 delay (allows block processing) (> blocksz + 1)
      // loop2
      make_ap (947, 0.5 /*overridden*/, 67), // 19
      make_delay (3191), // 20
      make_crossover2(), // 21
      make_quantizer(), // 22
      make_ap (887, -0.5 /*overridden*/), // 23
      make_delay (1049), // 24
      make_quantizer(), // 25
      make_ap (1367, 0.618), // 26
      make_hp (0.98), // 27
      make_block_delay (
        647)); // 28 delay (allows block processing) (> blocksz + 1)
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (T) ((io[i][0] + io[i][1]) * 0.5_r);
      auto mod   = (T) (0.25_r + (1_r - par.mod[i]) * 0.75_r);
      // ER + late lfo
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = (T) (T {lfo[1]} * (0.5_r + par.character[i] * 0.5_r));
      lfo3[i]  = (T) (T {lfo[2]} * mod);
      lfo4[i]  = (T) (T {lfo[3]} * mod);

      // decay fixup
      auto decay   = (T) (rev.one - par.decay[i]);
      decay        = (T) (rev.one - decay * decay);
      par.decay[i] = (T) (0.6_r + decay * 0.38_r);
    }
    // damp -----------------------------------
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.55f + upar.lf_amt * 0.4f);
    T fhi = load_float<T> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.6f + upar.hf_amt * 0.25f);

    rev.run (sl<0> {}, late_in);
    // diffusion -----------------------------
    rev.run (sl<1> {}, late_in);
    rev.run (sl<2> {}, late_in);
    rev.run (sl<3> {}, late_in);
    rev.run (sl<4> {}, late_in);

    // ER (first reverb, not exactly ER at the end...) -----------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    rev.fetch_block (sl<9> {}, er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = (T) (late_in[i] * 0.5_r + er2[i]);
    });
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run (sl<5> {}, er1, xspan {lfo2});
    rev.run (sl<6> {}, er1, flo, glo, fhi, rev.one, ghi);
    xspan_memcpy (er1b, er1);
    rev.run (sl<7> {}, er1b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<8> {}, er1b, xspan {lfo1});
    rev.push (sl<9> {}, er1b.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto e_l = (T) (-er1[i] * 0.66_r - er2[i] * 0.34_r);
      auto e_r = (T) (-er1[i] * 0.66_r + er2[i] * 0.34_r);
      er1[i]   = e_l;
      er2[i]   = e_r;
    }

    // Late -----------------------------
    arr    late_arr;
    arr_fb l_arr;
    arr_fb r_arr;
    arr    g_arr;

    auto late = xspan {late_arr.data(), io.size()};
    auto g    = xspan {g_arr.data(), io.size()};
    auto l    = xspan {l_arr.data(), io.size() + 1}; // +1: Feedback on head
    auto r    = xspan {r_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    rev.fetch_block (sl<18> {}, l, 1); // feedback, fetching block + 1 samples
    rev.fetch_block (sl<28> {}, r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (rev.one - er_amt) + er_sig * er_amt);
      g[i] = (T) (0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r));
    }
    r.cut_head (1); // drop feedback sample from previous block

    rev.run (sl<10> {}, late, xspan {lfo3}, g);
    rev.run (sl<11> {}, late);
    rev.run (sl<12> {}, late, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<13> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    rev.run (sl<14> {}, late, blank, [g] (uint i) { return -g[i]; });
    rev.run (sl<15> {}, late);
    rev.run (sl<16> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    rev.run (sl<17> {}, late);
    rev.push (sl<18> {}, late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (rev.one - er_amt) + er_sig * er_amt);
    }
    l.cut_head (1); // drop feedback sample from previous block
    rev.run (sl<19> {}, late, xspan {lfo4}, g);
    rev.run (sl<20> {}, late);
    rev.run (sl<21> {}, late, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<22> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    rev.run (sl<23> {}, late, blank, [g] (uint i) { return -g[i]; });
    rev.run (sl<24> {}, late);
    rev.run (sl<25> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    rev.run (sl<26> {}, late);
    rev.run (sl<27> {}, late);
    rev.push (sl<28> {}, late.to_const()); // feedback point

    // Mixdown
    auto v1 = make_array (&l[0], &er1[0]);
    auto v2 = make_array<T const*> (&r[0], &er2[0]);
    ep_crossfade<T> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return fastgrowth ((T) (v * 0.1_r), rev.one);
    });

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct plate1 : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<plate1, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f1 = 0.43f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.75f, 1.f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // left channel output conditioning part 1
      make_ap (71), // 0
      make_ap (116), // 1
      make_ap (172), // 2

      // right channel output conditioning part 2
      make_ap (72), // 3
      make_ap (117), // 4
      make_ap (175), // 5

      // left channel output conditioning part 2 (chorus)
      make_ap (1987, 0, 123), // 6
      make_lp(), // 7
      // right channel output conditioning part 2 (chorus)
      make_ap (2177, 0, 138), // 8
      make_lp(), // 9

      make_ap (277), // 10 // L last AP
      make_ap (274), // 11 // R last AP

      // m conditioning-diffusors
      make_lp (0.5f), // 12
      make_hp (0.87f), // 13
      make_ap (23, 0.8), // 14
      make_ap (130, 0.8), // 15
      make_ap (217, 0.8), // 16

      // s conditioning-diffusors
      make_lp (0.5f), // 17
      make_hp (0.87f), // 18
      make_ap (22, 0.8), // 19
      make_ap (131, 0.8), // 20
      make_ap (219, 0.8), // 21

      make_quantizer(), // 22
      make_quantizer(), // 23
      make_quantizer(), // 24
      make_quantizer(), // 25

      // block a iteration 1
      make_ap (153, 0.2), // 26
      make_ap (89, -0.04), // 27
      make_ap (60, 0.04), // 28
      make_delay (201), // 29

      // block b iteration 1
      make_delay (185), // 39

      // block c iteration 1
      make_ap (149, 0.2), // 31
      make_ap (83, 0.04), // 32
      make_ap (59, -0.04), // 33
      make_delay (225), // 34

      // block d iteration 1
      make_ap (167, -0.2), // 35
      make_ap (97, -0.04), // 36
      make_ap (67, 0.04), // 37
      make_delay (221), // 38

      make_quantizer(), // 39
      make_quantizer(), // 40
      make_quantizer(), // 41
      make_quantizer(), // 42

      // crossovers
      make_crossover2(), // 43
      make_crossover2(), // 44
      make_crossover2(), // 45
      make_crossover2(), // 46

      // block a iteration 2
      make_ap (119, 0.6), // 47
      make_ap (67, 0.6), // 48
      make_ap (47, -0.6), // 49
      make_block_delay (171), // feedback point // 50

      // block b iteration 2
      make_ap (114, -0.1), // nested 3x // 51
      make_ap (66, -0.04), // 52
      make_ap (47, -0.04), // 53
      make_ap (9), // 54
      make_block_delay (185), // feedback point // 55

      // block c iteration 2
      make_ap (116, -0.1), // nested 3x // 56
      make_ap (65, -0.04), // 57
      make_ap (46, -0.04), // 58
      make_ap (3), //  59
      make_block_delay (179), // feedback point // 60

      // block d iteration 2
      make_ap (121, -0.1), // nested 3x // 61
      make_ap (69, 0.04), // 62
      make_ap (47, 0.04), // 63
      make_block_delay (175) // feedback point // 64
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr_fb a_mem, d_mem;
    arr    b_mem, c_mem, tmp1, tmp2, tmp3;
    arr    lfo1, lfo2;

    xspan a {a_mem.data(), io.size() + 1};
    xspan b {b_mem.data(), io.size()};
    xspan c {c_mem.data(), io.size()};
    xspan d {d_mem.data(), io.size() + 1};

    // fetch a block from the output (a and d)
    rev.fetch_block (sl<50> {}, a, 1); // feedback, fetching block + 1 samples
    rev.fetch_block (sl<64> {}, d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // invert character first
      par.character[i] = (T) (rev.one - par.character[i]);
      tmp1[i]          = (T) (0.4_r + 0.3_r * par.character[i]);
      tmp2[i]          = (T) (-0.6_r * par.character[i]);
      tmp3[i]          = (T) (0.6_r * par.character[i]);
      // decay fixup
      auto d       = rev.one - par.decay[i];
      par.decay[i] = (T) (0.985_r - d * d * 0.21_r);
      // lfo
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = T {lfo[1]};
      // lfo3 = -lfo1, lfo4 = -lfo2
    }

    // output preconditioning
    rev.run (sl<0> {}, b, blank, tmp1);
    rev.run (sl<1> {}, b, blank, tmp2);
    rev.run (sl<2> {}, b, blank, tmp2);

    rev.run (sl<3> {}, c, blank, tmp1);
    rev.run (sl<4> {}, c, blank, tmp2);
    rev.run (sl<5> {}, c, blank, tmp2);

    // tmp1 and tmp2 free
    xspan cho1 {tmp1.data(), io.size()};
    xspan cho2 {tmp2.data(), io.size()};

    rev.run (sl<6> {}, cho1, b.to_const(), lfo1, tmp3);
    rev.run (sl<7> {}, cho1);
    rev.run (
      sl<8> {}, cho2, c.to_const(), [&] (uint i) { return -lfo1[i]; }, tmp3);
    rev.run (sl<9> {}, cho2);
    auto v1 = make_array (&b[0], &c[0]);
    auto v2 = make_array<T const*> (&cho1[0], &cho2[0]);
    ep_crossfade<T> (v1, v2, &par.mod[0], io.size(), [=] (auto v) {
      return fastgrowth ((T) (v * 0.5_r), rev.one);
    });

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      b[i] = (T) (b[i] + cho1[i]);
      c[i] = (T) (c[i] + cho2[i]);
    }

    rev.run (sl<10> {}, b, blank, tmp3);
    rev.run (sl<11> {}, c, blank, tmp3); // tmp3 free

    xspan i1 {tmp1.data(), io.size()};
    xspan i2 {tmp2.data(), io.size()};
    // output write + preparations
    span_visit (io, [&] (auto& spls, uint i) {
      i1[i]   = (T) ((spls[0] + spls[1]) * 0.25_r); // m
      i2[i]   = (T) ((spls[0] - spls[1]) * 0.25_r); // s
      spls[0] = b[i];
      spls[1] = c[i];
      b[i]    = (T) (0.8_r + 0.1_r * par.character[i]); // character 1 on b
      c[i]    = (T) (0.8_r * par.character[i]); // character 2 on c
    });

    // input preconditioning
    rev.run (sl<12> {}, i1);
    rev.run (sl<13> {}, i1);
    rev.run (sl<14> {}, i1, blank, b);
    rev.run (sl<15> {}, i1, blank, c);
    rev.run (sl<16> {}, i1, blank, c);

    rev.run (sl<17> {}, i2);
    rev.run (sl<18> {}, i2);
    rev.run (sl<19> {}, i2, blank, b);
    rev.run (sl<20> {}, i2, blank, c);
    rev.run (sl<21> {}, i2, blank, c);

    // fetch b and d feedback (a and d are already ready) 1 block only
    rev.fetch_block (sl<55> {}, b, 1); // feedback
    rev.fetch_block (sl<60> {}, c, 1); // feedback

    // quantized decay
    rev.run (sl<22> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<23> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<24> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<25> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // a single hadamard iteration can duplicate the range on one
      // of the channels, so halving once more (TODO needed?).
      b[i] = (T) (b[i] + tmp1[i] * 0.5_r);
      d[i] = (T) (d[i] + tmp2[i] * 0.5_r);
    }
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // channel a block 1
    rev.run (sl<26, 27, 28> {}, a);
    rev.run (sl<29> {}, a);

    // channel b block 1
    rev.run (sl<30> {}, b);

    // channel c block 1
    rev.run (sl<31, 32, 33> {}, c);
    rev.run (sl<34> {}, c);

    // channel d block 1
    rev.run (sl<35, 36, 37> {}, d);
    rev.run (sl<38> {}, d);

    // quantized decay
    rev.run (sl<39> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<40> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<41> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run (sl<42> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // block 2
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // swaps.
    a = xspan {d_mem.data(), io.size()};
    b = xspan {c_mem.data(), io.size()};
    c = xspan {a_mem.data(), io.size()};
    d = xspan {b_mem.data(), io.size()};

    // damp -----------------------------------
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.9f + upar.lf_amt * 0.09f);
    T fhi = load_float<T> (0.61f - upar.hf_amt * upar.hf_amt * 0.12f);
    T ghi = load_float<T> (0.6f + upar.hf_amt * 0.18f);

    rev.run (sl<43> {}, a, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<44> {}, b, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<45> {}, c, flo, glo, fhi, rev.one, ghi);
    rev.run (sl<46> {}, d, flo, glo, fhi, rev.one, ghi);

    // channel a block 2
    rev.run (sl<47> {}, a);
    rev.run (sl<48> {}, a);
    rev.run (sl<49> {}, a);
    rev.push (sl<50> {}, a.to_const());

    // channel b block 2
    rev.run (sl<51, 52, 53> {}, b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      T m     = fastgrowth (par.mod[i], rev.one);
      tmp1[i] = (T) (0.4_r + 0.4_r * lfo2[i] * m);
      tmp2[i] = (T) (0.6_r + 0.4_r * -lfo2[i] * m); // for block c
    }
    rev.run (sl<54> {}, b, tmp1);
    rev.push (sl<55> {}, b.to_const());

    // channel c block 2
    rev.run (sl<56, 57, 58> {}, c);
    rev.run (sl<59> {}, c, tmp2);
    rev.push (sl<60> {}, c.to_const());

    // channel d block 2
    rev.run (sl<61, 62, 63> {}, d);
    rev.push (sl<64> {}, d.to_const());
  }
};
//------------------------------------------------------------------------------
static constexpr std::array<u8, 32> get_room_ffwd_table()
{
  std::array<u8, 32> r {};
  // clamping at the block size
  r[0] = 32;
  assert (max_block_size == 32);
  r[1] = 33;
  r[2] = 34;
  for (uint i = 3; i < r.size(); ++i) {
    r[i] = primes_table_raw[8 + i];
  }
  return r;
}
//------------------------------------------------------------------------------
class room : public algorithm {
public:
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<room, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.43f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 1.f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp(),
      make_hp(),
      make_lp(),
      make_hp(), // 3
      make_variable_delay (ffwd_table[0], ffwd_table[16 + 0], true),
      make_variable_delay (ffwd_table[1], ffwd_table[16 + 1], true),
      make_variable_delay (ffwd_table[2], ffwd_table[16 + 2], true),
      make_variable_delay (ffwd_table[3], ffwd_table[16 + 3], true),
      make_variable_delay (ffwd_table[4], ffwd_table[16 + 4], true),
      make_variable_delay (ffwd_table[5], ffwd_table[16 + 5], true),
      make_variable_delay (ffwd_table[6], ffwd_table[16 + 6], true),
      make_variable_delay (ffwd_table[7], ffwd_table[16 + 7], true),
      make_variable_delay (ffwd_table[8], ffwd_table[16 + 8], true),
      make_variable_delay (ffwd_table[9], ffwd_table[16 + 9], true),
      make_variable_delay (ffwd_table[10], ffwd_table[16 + 10], true),
      make_variable_delay (ffwd_table[11], ffwd_table[16 + 11], true),
      make_variable_delay (ffwd_table[12], ffwd_table[16 + 12], true),
      make_variable_delay (ffwd_table[13], ffwd_table[16 + 13], true),
      make_variable_delay (ffwd_table[14], ffwd_table[16 + 14], true),
      make_variable_delay (ffwd_table[15], ffwd_table[16 + 15], true), // 19
      make_ap (913, 0.f, 413, interpolation::linear), // 20
      make_ap (491, 0.f, 213, interpolation::linear), // 21
      make_crossover2(), // 22
      make_ap (491, 0.f, 213, interpolation::linear), // 23
      make_ap (133, 0.f), // 24
      make_ap (191, 0.f), // 25
      make_ap (211, 0.f), // 26
      make_ap (311, 0.f, 73)); // 27
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr = block_arr<T>;

    arr   l_m, r_m;
    xspan l {l_m.data(), io.size()};
    xspan r {r_m.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      l[i] = (T) (spl[0] * 0.25_r);
      r[i] = (T) (spl[1] * 0.25_r);
    });

    T lo = load_float<T> (0.75f - upar.hf_amt * 0.2f);
    T hi = load_float<T> (0.9f + upar.lf_amt * 0.075f);

    rev.run (sl<0> {}, l, lo);
    rev.run (sl<1> {}, l, hi);
    rev.run (sl<2> {}, r, lo);
    rev.run (sl<3> {}, r, hi);

    // Mystran's feedforward diffusor
    mp11::mp_for_each<mp11::mp_from_sequence<std::make_index_sequence<16>>> (
      [&] (auto idx) {
        rev.run (sl<idx + 4> {}, (idx & 1) ? r : l, [&] (uint i) {
          return ffwd_table[idx + (uint) (par.decay[i] * 16_r)];
        });
        hadamard2 (make_array (l.data(), r.data()), l.size());
      });

    arr   m_m, ltdec, t1, t2, t3;
    xspan m {m_m.data(), io.size()};

    rev.fetch_block (sl<20> {}, xspan {t1.data(), io.size()}, 33);
    rev.fetch_block (sl<20> {}, xspan {t2.data(), io.size()}, 157);

    span_visit (m, [&] (auto& spl, uint i) {
      spl      = (T) ((l[i] + r[i]) * 0.25_r);
      auto inv = rev.one - par.decay[i];
      l[i]     = (T) (l[i] + t1[i] * inv);
      r[i]     = (T) (r[i] - t2[i] * inv);
      ltdec[i] = (T) (par.decay[i] * par.decay[i]);
      t1[i]    = (T) (0.39_r * ltdec[i]);
      t2[i]    = (T) (-0.28_r * ltdec[i]);
      t3[i]    = (T) (0.19_r * ltdec[i]);
    });
    // triple nested AP with crossover
    rev.run (
      sl<20, 21, 22, 23> {},
      m,
      par.character,
      t1,
      par.character,
      t2,
      load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f),
      load_float<T> (0.55f - upar.hf_amt * upar.hf_amt * 0.15f),
      load_float<T> (0.85f + upar.lf_amt * 0.13f),
      rev.one,
      load_float<T> (0.35f + upar.hf_amt * 0.14f),
      par.character,
      t3);

    span_visit (xspan {ltdec.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (T) (0.3_r + v * 0.3_r);
      t2[i] = (T) (-0.25_r - v * 0.3_r);
      t3[i] = (T) (0.25_r + v * 0.3_r);
    });
    rev.run (sl<24> {}, m, blank, t1);
    rev.run (sl<25> {}, m, blank, t2);
    rev.run (sl<26> {}, m, blank, t3);

    span_visit (xspan {par.mod.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (T) (T {tick_lfo<T> (lfo_obj)[0]} * v);
      t2[i] = (T) (0.365_r + v * 0.2_r);
    });
    rev.run (sl<27> {}, xspan {t3.data(), io.size()}, m.to_const(), t1, t2);
    span_visit (io, [&] (auto& spl, uint i) {
      spl[0] = (T) (l[i] + par.decay[i] * t3[i]);
      spl[1] = (T) (r[i] - par.decay[i] * m[i]);
    });
  }
  //----------------------------------------------------------------------------
  static constexpr std::array<u8, 32> ffwd_table {get_room_ffwd_table()};
};
//------------------------------------------------------------------------------
struct midifex49 : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<midifex49, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.2f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // diffusors
      make_ap (321, 0.5), // 0 PreAP
      make_ap (431, 0.5), // 1 PreAP
      make_ap (968, 0.5), // 2 PreAP
      make_ap (1620, 0.5), // 3 PreAP

      make_delay (21), // 4 L1
      make_block_delay (1010, 0), // 5 R1

      make_delay (1624), // 6 Loop
      make_ap (1992, 0.5, 17), // 7 Loop

      make_block_delay (1891, 0), // 8 L2
      make_block_delay (890, 0), // 9 R2

      make_delay (2110), // 10 Loop
      make_quantizer(), // 11
      make_ap (2371, 0.5), // 12 Loop nested allpass 1
      make_ap (1378, 0.2), // 13 Loop nested allpass 2

      make_block_delay (2003, 0), // 14 L3
      make_block_delay (671, 0), // 15 R3

      make_delay (2157 - max_block_size), // 16 Loop
      make_quantizer(), // 17
      make_crossover2(), // 18
      make_ap (2712, 0.5, 22), // 19 Loop nested allpass 1
      make_ap (1783, 0.2), // 20 Loop nested allpass 2

      make_block_delay (max_block_size) // 22 delay block (== blocksz + 1)
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr = block_arr<T>;

    arr   tmp_arr;
    xspan tmp {tmp_arr.data(), io.size()};
    arr   sig_arr;
    xspan sig {sig_arr.data(), io.size()};
    arr   l_arr;
    xspan l {l_arr.data(), io.size()};
    arr   r_arr;
    xspan r {r_arr.data(), io.size()};
    arr   lfo1_arr;
    xspan lfo1 {lfo1_arr.data(), io.size()};
    arr   lfo2_arr;
    xspan lfo2 {lfo2_arr.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      sig[i] = (T) (((spl[0] + spl[1]) * 0.25_r)); // gain = 0.5 + feedback = 1
      // LFO
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      // decay fixup
      auto decay   = (T) (rev.one - par.decay[i]);
      decay        = (T) (rev.one - decay * decay);
      par.decay[i] = (T) (0.1_r + decay * 0.8375_r);
    });

    rev.run (sl<0> {}, sig);
    rev.run (sl<1> {}, sig);
    rev.run (sl<2> {}, sig);
    rev.run (sl<3> {}, sig);

    // feedback handling, fetching the block with a negative offset of 1
    rev.fetch_block (sl<21> {}, tmp, 1);
    span_visit (sig, [&] (auto& v, uint i) { v = (T) (v + tmp[i]); });

    // 1st output point for L and R signal
    xspan_memcpy (l, sig); // delay LT a block -> might overlap, requires copy
    rev.run (sl<4> {}, l);
    rev.fetch_block (sl<5> {}, r, 0); // delay GT a block will never overlap
    rev.push (sl<5> {}, sig.to_const());

    // continuing the loop
    rev.run (sl<6> {}, sig);
    rev.run (sl<7> {}, sig, lfo1, blank);

    // 2nd output point for L and R signal
    rev.fetch_block (sl<8> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push (sl<8> {}, sig.to_const());
    rev.fetch_block (sl<9> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push (sl<9> {}, sig.to_const());

    // continuing the loop
    rev.run (sl<10> {}, sig);
    rev.run (sl<11> {}, sig, [&] (auto v, uint i) { return v * par.decay[i]; });
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.14_r);
    });
    rev.run (sl<12, 13> {}, sig, blank, blank, blank, tmp);

    // 3rd output point for L and R signal
    rev.fetch_block (sl<14> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push (sl<14> {}, sig.to_const());
    rev.fetch_block (sl<15> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push (sl<15> {}, sig.to_const());

    // continuing the loop
    rev.run (sl<16> {}, sig);
    rev.run (sl<17> {}, sig, [&] (auto v, uint i) { return v * par.decay[i]; });

    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + upar.lf_amt * 0.6f);
    T fhi = load_float<T> (0.8f - upar.hf_amt * upar.hf_amt * 0.1f);
    T ghi = load_float<T> (0.05f + upar.hf_amt * 0.6f);

    rev.run (sl<18> {}, sig, flo, glo, fhi, rev.one, ghi);
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.2_r);
    });
    rev.run (sl<19, 20> {}, sig, lfo2, blank, blank, tmp);
    // push to delay feedback
    rev.push (sl<21> {}, sig.to_const());

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
};

struct midifex50 : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<midifex50, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.3f + mod * 0.1f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // outputs first
      make_quantizer(), // 0
      make_ap (147, 0.5), // 1 L diff
      make_ap (43, 0.5), // 2 L diff
      make_ap (55, 0.5), // 3 L diff
      make_ap (249, 0.5), // 4 R diff
      make_ap (48, 0.5), // 5 R diff
      make_ap (21, 0.5), // 6 R diff
      // loop
      make_ap (13, 0.5), // 7
      make_ap (83, 0.5), // 8
      make_ap (116, 0.5), // 9
      make_ap (239, 0.5), // 10
      make_ap (339, 0.5), // 11
      make_ap (481, 0.5), // 12
      make_ap (555, 0.5), // 13
      make_ap (823, 0.5), // 14
      make_ap (999, 0.5), // 15
      make_ap (1100, 0.5), // 16
      make_ap (1347, 0.5), // 17
      make_ap (1563, 0.5), // 18
      make_ap (1841 - 32, 0.5), // 19
      make_block_delay (32), // 20 C1
      make_ap (2001 - 32, 0.5, 67), // 21
      make_block_delay (32), // 22 C2
      make_ap (2083 - 32, 0.5, 127), // 23
      make_crossover2(), // 24
      make_block_delay (max_block_size) // 25 (FB point) (== blocksz + 1)
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr_fb tmp1;
    xspan  loop {tmp1.data(), io.size() + 1};
    arr    l_arr;
    xspan  l {l_arr.data(), io.size()};
    arr    r_arr;
    xspan  r {r_arr.data(), io.size()};

    // xspan_memcpy (loop, io);
    rev.fetch_block (sl<25> {}, loop, 1);
    xspan_memcpy (r, loop.advanced (1));
    xspan_memcpy (l, loop.advanced (1));

    loop.cut_tail (1); // feedback samples
    rev.run (sl<0> {}, loop, [&] (T fb_spl, uint i) {
      // decay fixup
      auto decay = (T) (rev.one - par.decay[i]);
      decay      = (T) (rev.one - decay * decay);
      decay      = (T) - (0.3_r + decay * 0.45_r);
      return (fb_spl * decay) + ((io[i][0] + io[i][1]) * 0.25_r); // gain = 1
    }); // feedback + input summing with quantizer

    rev.run (sl<1> {}, l);
    rev.run (sl<2> {}, l);
    rev.run (sl<3> {}, l);

    rev.run (sl<4> {}, r);
    rev.run (sl<5> {}, r);
    rev.run (sl<6> {}, r);

    rev.run (sl<7> {}, loop);
    rev.run (sl<8> {}, loop);
    rev.run (sl<9> {}, loop);
    rev.run (sl<10> {}, loop);
    rev.run (sl<11> {}, loop);
    rev.run (sl<12> {}, loop);
    rev.run (sl<13> {}, loop);
    rev.run (sl<14> {}, loop);
    rev.run (sl<15> {}, loop);
    rev.run (sl<16> {}, loop);
    rev.run (sl<17> {}, loop);
    rev.run (sl<18> {}, loop);
    rev.run (sl<19> {}, loop);

    arr tmp2;
    arr tmp3;

    xspan c1 {tmp2.data(), loop.size()};
    xspan lfo2 {tmp3.data(), loop.size()};
    rev.fetch_block (sl<20> {}, c1, 0);
    rev.push (sl<20> {}, loop.to_const());
    auto lfo1 = loop; // data inserted already (lfo1 -> tmp1)
    loop      = c1; // avoid a copy. loop -> tmp2
    span_visit (l, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = rev.one - c;
      v      = (T) (k * v + loop[i] * c); // L done
      // unrelated but done here to skip one iteration
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
    });
    rev.run (sl<21> {}, loop, xspan {lfo1}, blank); // tmp1 free
    xspan c2 {tmp1.data(), loop.size()};
    rev.fetch_block (sl<22> {}, c2, 0);
    rev.push (sl<22> {}, loop.to_const());
    loop = c2; // avoid copy. tmp2 (loop) free.
    span_visit (r, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = rev.one - c;
      v      = (T) (k * v + loop[i] * c); // R done
    });
    // outputs done, dump now that they have been recently touched
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
    rev.run (sl<23> {}, loop, xspan {lfo2}, blank); // tmp3 free

    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + upar.lf_amt * 0.6f);
    T fhi = load_float<T> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.35f + upar.hf_amt * 0.3f);

    rev.run (sl<24> {}, loop, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<25> {}, loop.to_const());
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct dre2000a : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<dre2000a, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.75f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.12), // 1
      make_hp (0.99), // 2
      make_parallel_delay (1, 42, -g, 586, g, 1099, -g), // 3 (L)
      make_parallel_delay (1, 105, g, 490, -g, 1290, g), // 4 (R)

      make_comb (5719 - 1, 0.f, 43), // 5
      make_crossover2(), // 6
      make_parallel_delay (
        1, 640, -g, 1494, g, 2199, -g, 3122, g, 4135, -g, 4952, g), // 7

      make_comb (5779 - 1, 0.f, 44), // 8
      make_crossover2(), // 9
      make_parallel_delay (
        1,
        902,
        -g,
        1830,
        g,
        2528,
        -g,
        3641,
        g,
        4670,
        -g,
        5432,
        g), // 10

      make_comb (5905 - 1, 0.f, 42), // 11
      make_crossover2(), // 12
      make_parallel_delay (
        1,
        979,
        -g,
        2031,
        g,
        3267,
        -g,
        4104,
        g,
        5103,
        -g,
        5886,
        g), // 13
      // L
      make_ap (224, -0.7), // 14
      make_delay (311, 311), // 15
      make_ap (227, -0.7), // 16
      make_ap (1343, -0.7), // 17
      // R
      make_ap (224, -0.7), // 18
      make_delay (277, 400), // 19
      make_ap (91, -0.7), // 20
      make_ap (1182, -0.7)); // 21
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run (sl<1> {}, in);
    rev.run (sl<2> {}, l, in.to_const());
    xspan_memdump (r.data(), l);
    rev.run (sl<3> {}, l, overwrite);
    rev.run (sl<4> {}, r, overwrite);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60 (sl<5, 8, 11> {}, 0.9f + dec2 * 19.1f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.8f + upar.lf_amt * 0.18f);
    T fhi = load_float<T> (0.7f - upar.hf_amt * upar.hf_amt * 0.05f);
    T ghi = load_float<T> (dec2 * 0.25f + upar.hf_amt * 0.45f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block (sl<5> {}, comb_fb, lfo1, -gains[0]);
    rev.run (sl<6> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<5> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<7> {}, comb_fb.to_const(), overwrite, tank);

    rev.fetch_block (sl<8> {}, comb_fb, lfo2, -gains[1]);
    rev.run (sl<9> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<8> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<10> {}, comb_fb.to_const(), add_to, tank);

    rev.fetch_block (sl<11> {}, comb_fb, lfo3, -gains[2]);
    rev.run (sl<12> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<11> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<13> {}, comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      l[i]     = (T) (l[i] + tank[i]);
      r[i]     = (T) (r[i] + tank[i]);
      c        = (T) ((c - rev.one * 0.5_r) * 2_r);
      c        = (T) ((rev.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    auto stank = xspan {tank.data(), io.size()}.to_const();
    rev.run (sl<14> {}, l, stank);
    xspan er {tmp2.data(), io.size()};
    rev.run (sl<15> {}, er, in.to_const(), par.character);
    crossfade (l, er, eramt, rev.one);
    rev.run (sl<16> {}, l);
    rev.run (sl<17> {}, l);

    rev.run (sl<18> {}, r, stank);
    rev.run (sl<19> {}, er, in.to_const(), par.character);
    crossfade (r, er, eramt, rev.one);
    rev.run (sl<20> {}, r);
    rev.run (sl<21> {}, r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

struct dre2000b : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<dre2000b, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.15f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.15), // 1
      make_hp (0.985), // 2

      make_comb (3821 - 1, 0.f, 54), // 3
      make_crossover2(), // 4
      make_parallel_delay (
        1, 429, -g, 1000, g, 1472, -g, 2088, g, 2765, -g, 3311, g), // 5

      make_comb (4036 - 1, 0.f, 53), // 6
      make_crossover2(), // 7
      make_parallel_delay (
        1,
        616,
        -g,
        1225,
        g,
        1691,
        -g,
        2434,
        g,
        3122,
        -g,
        3631,
        g), // 8

      make_comb (4059 - 1, 0.f, 44.f), // 9
      make_crossover2(), // 10
      make_parallel_delay (
        1,
        657,
        -g,
        1359,
        g,
        2184,
        -g,
        2744,
        g,
        3411,
        -g,
        3934,
        g), // 11

      make_ap (282, -0.7), // 12
      make_ap (343, -0.7), // 13
      // L
      make_delay (311, 311), // 14
      make_ap (233, -0.7), // 15
      make_ap (273, -0.7), // 16
      make_ap (534, -0.7), // 17
      // R
      make_delay (277, 400), // 18
      make_ap (194, -0.7), // 19
      make_ap (426, -0.7), // 20
      make_ap (566, -0.7)); // 21
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run (sl<1> {}, in);
    rev.run (sl<2> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60 (sl<3, 6, 9> {}, 0.25f + dec2 * 10.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + upar.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.4f);
    T ghi = load_float<T> (0.4f + dec2 * 0.4f + upar.hf_amt * 0.15f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block (sl<3> {}, comb_fb, lfo1, -gains[0]);
    rev.run (sl<4> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<3> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<5> {}, comb_fb.to_const(), overwrite, tank);

    rev.fetch_block (sl<6> {}, comb_fb, lfo2, -gains[1]);
    rev.run (sl<7> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run (sl<8> {}, comb_fb.to_const(), add_to, tank);

    rev.fetch_block (sl<9> {}, comb_fb, lfo3, -gains[2]);
    rev.run (sl<10> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run (sl<11> {}, comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - rev.one * 0.5_r) * 2_r);
      c        = (T) ((rev.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    xspan stank {tank.data(), io.size()};
    rev.run (sl<12> {}, stank);
    rev.run (sl<13> {}, l, stank.to_const());
    xspan_memdump (r.data(), l);

    xspan er {tmp2.data(), io.size()};
    rev.run (sl<14> {}, er, in.to_const(), par.character);
    crossfade (l, er, eramt, rev.one);
    rev.run (sl<15> {}, l);
    rev.run (sl<16> {}, l);
    rev.run (sl<17> {}, l);

    rev.run (sl<18> {}, er, in.to_const(), par.character);
    crossfade (r, er, eramt, rev.one);
    rev.run (sl<19> {}, r);
    rev.run (sl<20> {}, r);
    rev.run (sl<21> {}, r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

struct dre2000c : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<dre2000c, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.75f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.18), // 1
      make_hp (0.985), // 2
      make_ap (115, -0.7f), // 3
      make_ap (160, -0.7f), // 4
      make_ap (231, -0.7f), // 5

      make_comb (3794 - 1, 0.f, 23), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 20 + 13, -g, 959, g, 1817, -g, 2855, g), // 8 R/L

      make_comb (3838 - 1, 0.f, 24), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 339, -g, 1309, g, 2271, -g, 3221, g), // 11 R/L

      make_comb (3861 - 1, 0.f, 25), // 12
      make_crossover2(), // 13
      make_parallel_delay (2, 130, -g, 1104, g, 2065, -g, 3391, g), // 14 R/L

      make_comb (3894 - 1, 0.f, 26), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 499, -g, 1445, g, 2071, -g, 2885, g), // 17 R/L
      // L
      make_ap (140, -0.7), // 18
      make_ap (160, -0.7), // 19
      make_delay (212, 237), // 20
      // R
      make_ap (160, -0.7), // 21
      make_ap (227, -0.7), // 22
      make_delay (213, 126) // 23
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run (sl<1> {}, in);
    rev.run (sl<2> {}, in);
    rev.run (sl<3> {}, in);
    rev.run (sl<4> {}, in);
    rev.run (sl<5> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60 (sl<6, 9, 12, 15> {}, 0.25f + dec2 * 20.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + upar.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.3f + dec2 * 0.25f + upar.hf_amt * 0.4499f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block (sl<6> {}, comb_fb, lfo1, gains[0]);
    rev.run (sl<7> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<8> {}, comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block (sl<9> {}, comb_fb, lfo2, gains[1]);
    rev.run (sl<10> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<11> {}, comb_fb.to_const(), add_to, r, l);

    rev.fetch_block (sl<12> {}, comb_fb, lfo3, gains[2]);
    rev.run (sl<13> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<12> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<14> {}, comb_fb.to_const(), add_to, r, l);

    rev.fetch_block (sl<15> {}, comb_fb, lfo4, gains[3]);
    rev.run (sl<16> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<15> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<17> {}, comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - rev.one * 0.5_r) * 2_r);
      c        = (T) ((rev.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run (sl<18> {}, l);
    rev.run (sl<19> {}, l);
    rev.run (sl<20> {}, lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt, rev.one);

    rev.run (sl<21> {}, r);
    rev.run (sl<22> {}, r);
    rev.run (sl<23> {}, in, par.character);
    crossfade (r, in, eramt, rev.one);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct dre2000d : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<dre2000d, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.35f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.185), // 1
      make_hp (0.99), // 2
      make_ap (115, -0.7f), // 3
      make_ap (214, -0.7f), // 4
      make_ap (365, -0.7f), // 5

      make_comb (5652 - 1, 0.f, 57), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 745, -g, 2164, g, 2922, -g, 5005, g), // 8 R/L

      make_comb (5679 - 1, 0.f), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 625, -g, 2080, g, 3523, -g, 4948, g), // 11 R/L

      make_comb (5689 - 1, 0.f), // 12
      make_crossover2(), // 13
      make_parallel_delay (
        2, 11 + 21, -g, 1435, g, 2769, -g, 4279, g), // 14 R/L

      make_comb (5711 - 1, 0.f, 57), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 59, -g, 1557, g, 3539, -g, 4304, g), // 17 R/L
      // L
      make_ap (214, -0.7, 30), // 18
      make_ap (469, -0.7), // 19
      make_delay (234, 337), // 20
      // R
      make_ap (245, -0.7, 29), // 21
      make_ap (426, -0.7), // 22
      make_delay (434, 247) // 23
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<T> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1, tmp2;
    xspan        in {in_mem.data(), io.size()};
    xspan        l {l_mem.data(), io.size()};
    xspan        r {r_mem.data(), io.size()};

    rev.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run (sl<1> {}, in);
    rev.run (sl<2> {}, in);
    rev.run (sl<3> {}, in);
    rev.run (sl<4> {}, in);
    rev.run (sl<5> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60 (sl<6, 9, 12, 15> {}, 0.25f + dec2 * 5.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.80f + upar.lf_amt * 0.199f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.4f + dec2 * 0.1f + upar.hf_amt * 0.49f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block (sl<6> {}, comb_fb, lfo1, gains[0]);
    rev.run (sl<7> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<8> {}, comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block (sl<9> {}, comb_fb, blank, gains[1]);
    rev.run (sl<10> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<11> {}, comb_fb.to_const(), add_to, r, l);

    rev.fetch_block (sl<12> {}, comb_fb, blank, gains[2]);
    rev.run (sl<13> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<12> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<14> {}, comb_fb.to_const(), add_to, r, l);

    rev.fetch_block (sl<15> {}, comb_fb, lfo4, gains[3]);
    rev.run (sl<16> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<15> {}, comb_fb, comb_fb.to_const(), in.to_const());
    rev.run (sl<17> {}, comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - rev.one * 0.5_r) * 2_r);
      c        = (T) ((rev.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run (sl<18> {}, l, lfo2);
    rev.run (sl<19> {}, l);
    rev.run (sl<20> {}, lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt, rev.one);

    rev.run (sl<21> {}, r, lfo3);
    rev.run (sl<22> {}, r);
    rev.run (sl<23> {}, in, par.character);
    crossfade (r, in, eramt, rev.one);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct rev5_l_hall : public algorithm {
  //----------------------------------------------------------------------------
  template <delay::data_type Dt>
  using engine = detail::lofiverb::engine<rev5_l_hall, Dt, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.35f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.18), // 1
      make_hp (0.987), // 2

      make_quantizer(), // 3
      make_ap (1257, 0.75f), // 4
      make_ap (978, 0.75f), // 5

      // rev main
      make_comb (2076, 0., 43), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 949, g, 837, g), // 8
      make_comb (2894, 0., 43), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 1330, g, 1412, g), // 11
      make_comb (3295, 0.), // 12
      make_crossover2(), // 13
      make_parallel_delay (2, 2393, g, 2704, g), // 14
      make_comb (3919, 0.), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 3263, g, 3011, g), // 17
      make_comb (4570, 0.), // 18
      make_crossover2(), // 19
      make_parallel_delay (2, 3667, g, 3667, g), // 20
      // sub reverb cond
      make_variable_delay (533, 533 + 2411), // 21
      make_ap (1212, 0.75f), // 22
      make_ap (960, 0.75f), // 23
      // sub reverb
      make_comb (3169, 0., 44), // 24
      make_crossover2(), // 25
      make_parallel_delay (2, 1757, g, 1644, g), // 26
      make_comb (3753, 0., 44), // 27
      make_crossover2(), // 28
      make_parallel_delay (2, 1900, g, 1981, g), // 29
      make_comb (4280, 0.), // 30
      make_crossover2(), // 31
      make_parallel_delay (2, 2838, g, 3148, g), // 32
      make_comb (4491, 0.), // 33
      make_crossover2(), // 34
      make_parallel_delay (2, 3798, g, 3545, g), // 35
      make_comb (5091, 0.), // 36
      make_crossover2(), // 37
      make_parallel_delay (2, 4298, g, 4298, g), // 38
      make_ap (682, 0.75f), // 39
      make_ap (830, 0.75f), // 40
      make_ap (695, 0.75f), // 41
      make_ap (844, 0.75f), // 42
      make_parallel_delay (
        2,
        18 + 14,
        0.95f,
        18 + 14,
        0.95f,
        895,
        0.95f,
        695,
        0.95f,
        1238,
        0.75f,
        1487,
        0.75f,
        1870,
        0.75f,
        1936,
        0.75f,
        2724,
        0.5f,
        2635,
        0.5f,
        3365,
        0.25f,
        3609,
        0.25f) // 43 L/R
    );
  }
  //----------------------------------------------------------------------------
  template <class T, delay::data_type Dt>
  static void process_block (
    engine<Dt>&                  rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<T> in_mem, main_mem, sub_mem, l_mem, r_mem, lfo1, lfo2, lfo3,
      lfo4, tmp1, tmp2;
    xspan main {main_mem.data(), io.size()};
    xspan sub {sub_mem.data(), io.size()};
    xspan in {in_mem.data(), io.size()};
    xspan l {l_mem.data(), io.size()};
    xspan r {r_mem.data(), io.size()};

    rev.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run (sl<1> {}, in);
    rev.run (sl<2> {}, in);
    xspan_memdump (main.data(), in);
    xspan_memdump (sub.data(), in);
    rev.run (sl<3> {}, main, [&] (auto v, uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      tmp2[i] = (T) (0.35_r + 0.7_r * clamp (c, 0._r, 0.4999_r)); // cg2
      return v * 0.25_r;
    });
    rev.run (sl<4> {}, main, blank, tmp1);
    rev.run (sl<5> {}, main, blank, tmp2);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains = rev.get_gain_for_rt60 (
      sl<6, 9, 12, 15, 18> {}, 1.f + dec2 * 10.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.9f + upar.lf_amt * 0.09f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.1f);
    T ghi = load_float<T> (0.75f + upar.hf_amt * 0.2f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block (sl<6> {}, comb_fb, lfo1, gains[0]);
    rev.run (sl<7> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.85_r));
    rev.push (sl<6> {}, comb_fb, comb_fb.to_const(), main.to_const());
    rev.run (sl<8> {}, comb_fb.to_const(), overwrite, l, r);

    rev.fetch_block (sl<9> {}, comb_fb, lfo3, gains[1]);
    rev.run (sl<10> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.9_r));
    rev.push (sl<9> {}, comb_fb, comb_fb.to_const(), main.to_const());
    rev.run (sl<11> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<12> {}, comb_fb, blank, gains[2]);
    rev.run (sl<13> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.95_r));
    rev.push (sl<12> {}, comb_fb, comb_fb.to_const(), main.to_const());
    rev.run (sl<14> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<15> {}, comb_fb, blank, gains[3]);
    rev.run (sl<16> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.97_r));
    rev.push (sl<15> {}, comb_fb, comb_fb.to_const(), main.to_const());
    rev.run (sl<17> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<18> {}, comb_fb, blank, gains[4]);
    rev.run (sl<19> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<18> {}, comb_fb, comb_fb.to_const(), main.to_const());
    rev.run (sl<20> {}, comb_fb.to_const(), add_to, l, r);

    rev.run (sl<21> {}, sub, par.character);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      // tmp2 still not overwritten...
      // tmp2[i] = (T) (0.35_r + 0.999_r * clamp (c, 0._r, 0.333_r)); //
      // cg2
    }
    rev.run (sl<22> {}, main, tmp2);
    rev.run (sl<23> {}, main, tmp1);
    gains = rev.get_gain_for_rt60 (
      sl<24, 27, 30, 33, 36> {}, 1.f + dec2 * 6.f, srate);
    flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    glo = load_float<T> (0.9f + upar.lf_amt * 0.04f);
    fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    ghi = load_float<T> (0.7f + upar.hf_amt * 0.25f);

    // sub  reverb
    rev.fetch_block (sl<24> {}, comb_fb, lfo2, gains[0]);
    rev.run (sl<25> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.9_r));
    rev.push (sl<24> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run (sl<26> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<27> {}, comb_fb, lfo4, gains[1]);
    rev.run (sl<28> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.93_r));
    rev.push (sl<27> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run (sl<29> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<30> {}, comb_fb, blank, gains[2]);
    rev.run (sl<31> {}, comb_fb, flo, glo, fhi, rev.one, (T) (ghi * 0.97_r));
    rev.push (sl<30> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run (sl<32> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<33> {}, comb_fb, blank, gains[3]);
    rev.run (sl<34> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<33> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run (sl<35> {}, comb_fb.to_const(), add_to, l, r);

    rev.fetch_block (sl<36> {}, comb_fb, blank, gains[4]);
    rev.run (sl<37> {}, comb_fb, flo, glo, fhi, rev.one, ghi);
    rev.push (sl<36> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run (sl<38> {}, comb_fb.to_const(), add_to, l, r);

    rev.run (sl<39> {}, l, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo1[i]);
    });
    rev.run (sl<40> {}, l, blank, [&] (uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.325_r + 0.999_r * clamp (c, 0._r, 0.333_r)); // cg1
      return (T) (tmp1[i] + 0.05_r * lfo2[i]);
    });
    rev.run (sl<41> {}, r, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo3[i]);
    });
    rev.run (sl<42> {}, r, blank, [&] (uint i) {
      return (T) (tmp1[i] + 0.05_r * lfo4[i]);
    });

    rev.run (sl<43> {}, in.to_const(), add_to, l, r); // ER
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

}}} // namespace artv::detail::lofiverb
