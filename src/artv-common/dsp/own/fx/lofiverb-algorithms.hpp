#pragma once

#include <type_traits>

#include "artv-common/dsp/own/fx/lofiverb-engine.hpp"
#include "artv-common/dsp/own/parts/oscillators/lfo.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/primes_table.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
static constexpr uint max_block_size = 32;

// this is using 16 bits fixed-point arithmetic, positive values can't
// represent one, so instead of correcting everywhere the parameters are
// scaled instead to never reach 1.
static constexpr auto fixpt_max_flt = fixpt_t::max_float();

struct algorithm {

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
  static constexpr auto one()
  {
    if constexpr (std::is_floating_point_v<T>) {
      return fixpt_max_flt;
    }
    else {
      return fixpt_t::max();
    }
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
  // 1 selects s2, 0 dst
  template <class T, class U, class V>
  static void crossfade (xspan<T> dst, U s2, V ctrl)
  {
    span_visit (dst, [&] (auto spl, uint i) {
      auto ctrl2 = (T) (one<T>() - ctrl[i]);
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
};

class debug : private algorithm {
public:
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<debug, max_block_size>;
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
      make_ap (147, 0.5), // 0
      make_ap (183, 0.4), // 1
      make_ap (389, 0.3), // 2
      make_block_delay (max_block_size) // 3
    );
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void process_block (
    engine&                      rev,
    lfo<4>&                      lfo_obj,
    xspan<std::array<T, 2>>      io,
    smoothed_parameters<T>&      par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<T>;
    using arr_fb = fb_block_arr<T>;

    arr    in_arr;
    arr_fb loop_arr;

    auto in   = xspan {in_arr.data(), io.size()};
    auto loop = xspan {loop_arr.data(), io.size() + 1};

    span_visit (in, [&] (auto& v, uint i) {
      v = (T) ((io[i][0] + io[i][1]) * 0.25_r);
    });
    rev.fetch_block<3> (loop, 1);
    span_visit (io, [&] (auto& v, uint i) { v[0] = v[1] = loop[i + 1]; });
    loop.cut_tail (1);
    span_visit (loop, [&] (auto& v, uint i) {
      v = (T) (in[i] + v * (0.7_r + par.decay[i] * 0.28_r));
    });
    rev.run<0, 1, 2> (loop);
    rev.push<3> (loop.to_const()); // feedback point
  }
  //----------------------------------------------------------------------------
};

class abyss : private algorithm {
public:
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<abyss, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.6_r + decay * 0.38_r);
    }
    // damp -----------------------------------
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.55f + upar.lf_amt * 0.4f);
    T fhi = load_float<T> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.6f + upar.hf_amt * 0.25f);

    rev.run<0> (late_in);
    // diffusion -----------------------------
    rev.run<1> (late_in);
    rev.run<2> (late_in);
    rev.run<3> (late_in);
    rev.run<4> (late_in);

    // ER (first reverb, not exactly ER at the end...) -----------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    rev.fetch_block<9> (er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = (T) (late_in[i] * 0.5_r + er2[i]);
    });
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run<5> (er1, xspan {lfo2});
    rev.run<6> (er1, flo, glo, fhi, one<T>(), ghi);
    xspan_memcpy (er1b, er1);
    rev.run<7> (er1b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<8> (er1b, xspan {lfo1});
    rev.push<9> (er1b.to_const()); // feedback point

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
    rev.fetch_block<18> (l, 1); // feedback, fetching block + 1 samples
    rev.fetch_block<28> (r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
      g[i] = (T) (0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r));
    }
    r.cut_head (1); // drop feedback sample from previous block

    rev.run<10> (late, xspan {lfo3}, g);
    rev.run<11> (late);
    rev.run<12> (late, flo, glo, fhi, one<T>(), ghi);
    rev.run<13> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<14> (late, blank, [g] (uint i) { return -g[i]; });
    rev.run<15> (late);
    rev.run<16> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<17> (late);
    rev.push<18> (late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (T) (loopsig * (one<T>() - er_amt) + er_sig * er_amt);
    }
    l.cut_head (1); // drop feedback sample from previous block
    rev.run<19> (late, xspan {lfo4}, g);
    rev.run<20> (late);
    rev.run<21> (late, flo, glo, fhi, one<T>(), ghi);
    rev.run<22> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<23> (late, blank, [g] (uint i) { return -g[i]; });
    rev.run<24> (late);
    rev.run<25> (late, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<26> (late);
    rev.run<27> (late);
    rev.push<28> (late.to_const()); // feedback point

    // Mixdown
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto er_amt     = par.character[i] * 0.2_r;
      auto e_l        = (-er1[i] * 0.66_r - er2[i] * 0.34_r) * er_amt;
      auto e_r        = (-er1[i] * 0.66_r + er2[i] * 0.34_r) * er_amt;
      auto direct_amt = one<T>() - er_amt;
      io[i][0]        = (T) (l[i] * direct_amt + e_l);
      io[i][1]        = (T) (r[i] * direct_amt + e_r);
    }
  }
  //----------------------------------------------------------------------------
};

struct small_space : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<small_space, max_block_size>;
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f1 = 0.1f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
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
      // left channel output conditioning
      make_ap (71), // 0
      make_ap (116), // 1
      make_ap (172), // 2
      make_ap (277), // 3

      // right channel output conditioning
      make_ap (72), // 4
      make_ap (117), // 5
      make_ap (175), // 6
      make_ap (274), // 7

      // m diffusors
      make_ap (23, 0.8), // 8
      make_ap (130, 0.8), // 9
      make_ap (217, 0.8), // 10
      // s diffusors
      make_ap (22, 0.8), // 11
      make_ap (131, 0.8), // 12
      make_ap (219, 0.8), // 13

      make_quantizer(), // 14
      make_quantizer(), // 15
      make_quantizer(), // 16
      make_quantizer(), // 17
      // block a iteration 1
      make_ap (153), // nested 3x (mod g) // 18
      make_ap (89, -0.04), // 19
      make_ap (60, 0.04), // 20
      make_delay (201), // 21

      // block b iteration 1
      make_ap (139), // nested 3x (mod g) // 22
      make_ap (79, -0.04), // 23
      make_ap (53, 0.04), // 24
      make_delay (185), // mod delay // 25

      // block c iteration 1
      make_ap (149), // nested 3x (mod g) // 26
      make_ap (83, 0.04), // 27
      make_ap (59, -0.04), // 28
      make_delay (193), // 29

      // block d iteration 1
      make_ap (167), // nested 3x (mod g) // 30
      make_ap (97, -0.04), // 31
      make_ap (67, 0.04), // 32
      make_delay (221), // mod delay // 33

      make_quantizer(), // 34
      make_quantizer(), // 35
      make_quantizer(), // 36
      make_quantizer(), // 37

      // block a iteration 2
      make_crossover2(), // 38
      make_ap (113, 0.1), // nested 3x // 39
      make_ap (67, 0.04), // 40
      make_ap (47, -0.04), // nested end // 41
      make_delay (122), // 42
      make_ap (119, 0.6), // 43
      make_ap (67, 0.6), // 44
      make_ap (47, -0.6), // 45
      make_block_delay (32), // feedback point // 46

      // block b iteration 2
      make_crossover2(), // 47
      make_ap (114, -0.1), // nested 3x // 48
      make_ap (66, -0.04), // 49
      make_ap (47, -0.04), // 50
      make_ap (9, 0, 2), // single. modulated both gain and time // 51
      make_block_delay (149), // feedback point // 52

      // block c iteration 2
      make_crossover2(), // 53
      make_ap (116, -0.1), // nested 3x // 54
      make_ap (65, -0.04), // 55
      make_ap (46, -0.04), // 56
      make_ap (9, 0, 3), // single. modulated both gain and time // 57
      make_block_delay (151), // feedback point // 58

      // block d iteration 2
      make_crossover2(), // 59
      make_ap (121, -0.1), // nested 3x // 60
      make_ap (69, 0.04), // 61
      make_ap (47, 0.04), // 62
      make_block_delay (157) // feedback point // 63
    );
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void process_block (
    engine&                      rev,
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
    rev.fetch_block<46> (a, 1); // feedback, fetching block + 1 samples
    rev.fetch_block<63> (d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // invert character first
      par.character[i] = (T) (one<T>() - par.character[i]);
      tmp1[i]          = (T) (0.4_r + 0.3_r * par.character[i]);
      tmp2[i]          = (T) (-0.6_r * par.character[i]);
      tmp3[i]          = (T) (0.6_r * par.character[i]);
      // decay fixup
      auto d       = one<T>() - par.decay[i];
      par.decay[i] = (T) (0.9_r - d * d * 0.2_r);
    }

    // output preconditioning
    rev.run<0> (b, 0, tmp1);
    rev.run<1> (b, 0, tmp2);
    rev.run<2> (b, 0, tmp2);
    rev.run<3> (b, 0, tmp3);

    rev.run<4> (c, 0, tmp1);
    rev.run<5> (c, 0, tmp2);
    rev.run<6> (c, 0, tmp2);
    rev.run<7> (c, 0, tmp3);

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
    rev.run<8> (i1, 0, b);
    rev.run<9> (i1, 0, c);
    rev.run<10> (i1, 0, c);

    rev.run<11> (i2, 0, b);
    rev.run<12> (i2, 0, c);
    rev.run<13> (i2, 0, c);

    // fetch b and d feedback (a and d are already ready) 1 block only
    rev.fetch_block<52> (b, 1); // feedback
    rev.fetch_block<58> (c, 1); // feedback

    // quantized decay
    rev.run<14> (a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<15> (b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<16> (c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<17> (d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // a single hadamard iteration can duplicate the range on one
      // of the channels, so halving once more.
      b[i] = (T) (b[i] + tmp1[i] * 0.5_r);
      d[i] = (T) (d[i] + tmp2[i] * 0.5_r);
    }
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = T {lfo[0]};
      lfo2[i]  = T {lfo[1]};
      // lfo3 = -lfo1, lfo4 = -lfo2
      // mod fixup
      par.mod[i] = (T) (one<T>() - par.mod[i]);
      par.mod[i] = (T) (one<T>() - (par.mod[i] * par.mod[i]));
      // preparing modlations for block a
      tmp1[i] = (T) (0.1_r + 0.1_r * lfo1[i]);
    }
    // channel a block 1
    rev.run<18, 19, 20> (a, blank, tmp1);
    rev.run<21> (a);

    // channel b block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (-0.1_r - 0.1_r * -lfo1[i]); // lfo3 = -lfo1
    }
    rev.run<22, 23, 24> (b, 0, tmp1);
    rev.run<25> (b);

    // channel c block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (0.1_r + 0.1_r * -lfo2[i]); // lfo4 = -lfo2
    }
    rev.run<26, 27, 28> (c, 0, tmp1);
    rev.run<29> (c);

    // channel d block 1
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (0.1_r - 0.1_r * lfo2[i]);
    }
    rev.run<30, 31, 32> (d, 0, tmp1);
    rev.run<33> (d);

    // quantized decay
    rev.run<34> (a, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<35> (b, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<36> (c, [&] (auto v, uint i) { return v * par.decay[i]; });
    rev.run<37> (d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // block 2
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // swaps.
    a = xspan {d_mem.data(), io.size()};
    b = xspan {c_mem.data(), io.size()};
    c = xspan {a_mem.data(), io.size()};
    d = xspan {b_mem.data(), io.size()};

    // damp -----------------------------------
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.87f + upar.lf_amt * 0.1299f);
    T fhi = load_float<T> (0.75f - upar.hf_amt * upar.hf_amt * 0.55f);
    T ghi = load_float<T> (0.7f + upar.hf_amt * 0.24f);

    // channel a block 2
    rev.run<38> (a, flo, glo, fhi, one<T>(), ghi);
    rev.run<39, 40, 41> (a);
    rev.run<42> (a);
    rev.run<43> (a);
    rev.run<44> (a);
    rev.run<45> (a);
    rev.push<46> (a.to_const());

    // channel b block 2
    rev.run<47> (b, flo, glo, fhi, one<T>(), ghi);
    rev.run<48, 49, 50> (b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (lfo1[i] * par.mod[i]);
      tmp2[i] = (T) (0.4_r + 0.4_r * lfo2[i] * par.mod[i]);
    }
    rev.run<51> (b, tmp1, tmp2);
    rev.push<52> (b.to_const());

    // channel c block 2
    rev.run<53> (c, flo, glo, fhi, one<T>(), ghi);
    rev.run<54, 55, 56> (c);
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = (T) (-tmp1[i]); // lfo3 = -lfo1.
      tmp2[i] = (T) (0.6_r - 0.4_r * -lfo2[i] * par.mod[i]); // lfo4 = -lfo2
    }
    rev.run<57> (c, tmp1, tmp2);
    rev.push<58> (c.to_const());

    // channel d block 2
    rev.run<59> (d, flo, glo, fhi, one<T>(), ghi);
    rev.run<60, 61, 62> (d);
    rev.push<63> (d.to_const());
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
class room : private algorithm {
public:
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<room, max_block_size>;
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
      make_ap (211, 0.f), // 22
      make_crossover2(), // 23
      make_ap (133, 0.f), // 24
      make_ap (191, 0.f), // 25
      make_ap (211, 0.f), // 26
      make_ap (311, 0.f, 73)); // 27
  }
  //----------------------------------------------------------------------------
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (l, lo);
    rev.run<1> (l, hi);
    rev.run<2> (r, lo);
    rev.run<3> (r, hi);

    // Mystran's diffusor
    mp11::mp_for_each<mp11::mp_from_sequence<std::make_index_sequence<16>>> (
      [&] (auto idx) {
        rev.run<idx + 4> ((idx & 1) ? l : r, [&] (uint i) {
          return ffwd_table[idx + (uint) (par.decay[i] * 16_r)];
        });
        hadamard2 (make_array (l.data(), r.data()), l.size());
      });

    arr   m_m, ltdec, t1, t2, t3;
    xspan m {m_m.data(), io.size()};

    rev.fetch_block<20> (xspan {t1.data(), io.size()}, 33);
    rev.fetch_block<20> (xspan {t2.data(), io.size()}, 157);

    span_visit (m, [&] (auto& spl, uint i) {
      spl      = (T) ((l[i] + r[i]) * 0.25_r);
      auto inv = one<T>() - par.decay[i];
#if 1
      l[i] = r[i] = T {0};
      spl         = (T) ((io[i][0] + io[i][1]) * 0.25_r);
#endif
#if 0
      l[i]     = (T) (l[i] + t1[i] * inv);
      r[i]     = (T) (r[i] - t2[i] * inv);
#endif
      ltdec[i] = (T) (par.decay[i] * par.decay[i]);
      t1[i]    = (T) (0.39_r * ltdec[i]);
      t2[i]    = (T) (-0.28_r * ltdec[i]);
      t3[i]    = (T) (0.19_r * ltdec[i]);
    });
    // triple nested AP with crossover
    rev.run<20, 21, 22, 23> (
      m,
      par.character,
      t1,
      par.character,
      t2,
      blank,
      t3,
      load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f),
      load_float<T> (0.85f + upar.lf_amt * 0.13f),
      load_float<T> (0.55f - upar.hf_amt * upar.hf_amt * 0.15f),
      one<T>(),
      load_float<T> (0.35f + upar.hf_amt * 0.14f));

#if 0
    span_visit (xspan {ltdec.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (T) (0.3_r + v * 0.3_r);
      t2[i] = (T) (-0.25_r - v * 0.3_r);
      t3[i] = (T) (0.25_r + v * 0.3_r);
    });
    rev.run<24> (m, blank, t1);
    rev.run<25> (m, blank, t2);
    rev.run<26> (m, blank, t3);

    span_visit (xspan {par.mod.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (T) (T {tick_lfo<T> (lfo_obj)[0]} * v);
      t2[i] = (T) (0.365_r + v * 0.2_r);
    });
    rev.run<27> (xspan {t3.data(), io.size()}, m.to_const(), t1, t2);
#else
    xspan_memcpy<T> (t3, m);
#endif
    span_visit (io, [&] (auto& spl, uint i) {
      spl[0] = (T) (l[i] + par.decay[i] * t3[i]);
      spl[1] = (T) (r[i] - par.decay[i] * m[i]);
    });
  }
  //----------------------------------------------------------------------------
  static constexpr auto ffwd_table {get_room_ffwd_table()};
};
//------------------------------------------------------------------------------
struct midifex49 : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<midifex49, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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
      auto decay   = (T) (one<T>() - par.decay[i]);
      decay        = (T) (one<T>() - decay * decay);
      par.decay[i] = (T) (0.1_r + decay * 0.8375_r);
    });

    rev.run<0> (sig);
    rev.run<1> (sig);
    rev.run<2> (sig);
    rev.run<3> (sig);

    // feedback handling, fetching the block with a negative offset of 1
    rev.fetch_block<21> (tmp, 1);
    span_visit (sig, [&] (auto& v, uint i) { v = (T) (v + tmp[i]); });

    // 1st output point for L and R signal
    xspan_memcpy (l, sig); // delay LT a block -> might overlap, requires copy
    rev.run<4> (l);
    rev.fetch_block<5> (r, 0); // delay GT a block will never overlap
    rev.push<5> (sig.to_const());

    // continuing the loop
    rev.run<6> (sig);
    rev.run<7> (sig, lfo1, blank);

    // 2nd output point for L and R signal
    rev.fetch_block<8> (tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<8> (sig.to_const());
    rev.fetch_block<9> (tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) ((v + tmp[i]) * (2_r / 3_r));
    });
    rev.push<9> (sig.to_const());

    // continuing the loop
    rev.run<10> (sig);
    rev.run<11> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.14_r);
    });
    rev.run<12, 13> (sig, blank, blank, blank, tmp);

    // 3rd output point for L and R signal
    rev.fetch_block<14> (tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<14> (sig.to_const());
    rev.fetch_block<15> (tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (T) (v + tmp[i] * (1_r / 3_r));
    });
    rev.push<15> (sig.to_const());

    // continuing the loop
    rev.run<16> (sig);
    rev.run<17> (sig, [&] (auto v, uint i) { return v * par.decay[i]; });

    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + upar.lf_amt * 0.6f);
    T fhi = load_float<T> (0.8f - upar.hf_amt * upar.hf_amt * 0.1f);
    T ghi = load_float<T> (0.05f + upar.hf_amt * 0.6f);

    rev.run<18> (sig, flo, glo, fhi, one<T>(), ghi);
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (T) (par.character[i] * 0.2_r);
    });
    rev.run<19, 20> (sig, lfo2, blank, blank, tmp);
    // push to delay feedback
    rev.push<21> (sig.to_const());

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
};

struct midifex50 : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<midifex50, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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
    rev.fetch_block<25> (loop, 1);
    xspan_memcpy (r, loop.advanced (1));
    xspan_memcpy (l, loop.advanced (1));

    loop.cut_tail (1); // feedback samples
    rev.run<0> (loop, [&] (T fb_spl, uint i) {
      // decay fixup
      auto decay = (T) (one<T>() - par.decay[i]);
      decay      = (T) (one<T>() - decay * decay);
      decay      = (T) - (0.3_r + decay * 0.45_r);
      return (fb_spl * decay) + ((io[i][0] + io[i][1]) * 0.25_r); // gain = 1
    }); // feedback + input summing with quantizer

    rev.run<1> (l);
    rev.run<2> (l);
    rev.run<3> (l);

    rev.run<4> (r);
    rev.run<5> (r);
    rev.run<6> (r);

    rev.run<7> (loop);
    rev.run<8> (loop);
    rev.run<9> (loop);
    rev.run<10> (loop);
    rev.run<11> (loop);
    rev.run<12> (loop);
    rev.run<13> (loop);
    rev.run<14> (loop);
    rev.run<15> (loop);
    rev.run<16> (loop);
    rev.run<17> (loop);
    rev.run<18> (loop);
    rev.run<19> (loop);

    arr tmp2;
    arr tmp3;

    xspan c1 {tmp2.data(), loop.size()};
    xspan lfo2 {tmp3.data(), loop.size()};
    rev.fetch_block<20> (c1, 0);
    rev.push<20> (loop.to_const());
    auto lfo1 = loop; // data inserted already (lfo1 -> tmp1)
    loop      = c1; // avoid a copy. loop -> tmp2
    span_visit (l, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = one<T>() - c;
      v      = (T) (k * v + loop[i] * c); // L done
      // unrelated but done here to skip one iteration
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
    });
    rev.run<21> (loop, xspan {lfo1}, blank); // tmp1 free
    xspan c2 {tmp1.data(), loop.size()};
    rev.fetch_block<22> (c2, 0);
    rev.push<22> (loop.to_const());
    loop = c2; // avoid copy. tmp2 (loop) free.
    span_visit (r, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = one<T>() - c;
      v      = (T) (k * v + loop[i] * c); // R done
    });
    // outputs done, dump now that they have been recently touched
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
    rev.run<23> (loop, xspan {lfo2}, blank); // tmp3 free

    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.35f + upar.lf_amt * 0.6f);
    T fhi = load_float<T> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.35f + upar.hf_amt * 0.3f);

    rev.run<24> (loop, flo, glo, fhi, one<T>(), ghi);
    rev.push<25> (loop.to_const());
  }
  //----------------------------------------------------------------------------
};

struct dre2000a : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<dre2000a, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (l, in.to_const());
    xspan_memdump (r.data(), l);
    rev.run<3> (l, overwrite);
    rev.run<4> (r, overwrite);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 5, 9, 13> (0.9f + dec2 * 19.1f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.8f + upar.lf_amt * 0.18f);
    T fhi = load_float<T> (0.7f - upar.hf_amt * upar.hf_amt * 0.05f);
    T ghi = load_float<T> (dec2 * 0.25f + upar.hf_amt * 0.45f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<5> (comb_fb, lfo1, -gains[0]);
    rev.run<6> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<5> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<7> (comb_fb.to_const(), overwrite, tank);

    rev.fetch_block<8> (comb_fb, lfo2, -gains[1]);
    rev.run<9> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<8> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<10> (comb_fb.to_const(), add_to, tank);

    rev.fetch_block<11> (comb_fb, lfo3, -gains[2]);
    rev.run<12> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<11> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<13> (comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      l[i]     = (T) (l[i] + tank[i]);
      r[i]     = (T) (r[i] + tank[i]);
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    auto stank = xspan {tank.data(), io.size()}.to_const();
    rev.run<14> (l, stank);
    xspan er {tmp2.data(), io.size()};
    rev.run<15> (er, in.to_const(), par.character);
    crossfade (l, er, eramt);
    rev.run<16> (l);
    rev.run<17> (l);

    rev.run<18> (r, stank);
    rev.run<19> (er, in.to_const(), par.character);
    crossfade (r, er, eramt);
    rev.run<20> (r);
    rev.run<21> (r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

struct dre2000b : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<dre2000b, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 3, 7, 11> (0.25f + dec2 * 10.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + upar.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.4f);
    T ghi = load_float<T> (0.4f + dec2 * 0.4f + upar.hf_amt * 0.15f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<3> (comb_fb, lfo1, -gains[0]);
    rev.run<4> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<3> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<5> (comb_fb.to_const(), overwrite, tank);

    rev.fetch_block<6> (comb_fb, lfo2, -gains[1]);
    rev.run<7> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<6> (comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run<8> (comb_fb.to_const(), add_to, tank);

    rev.fetch_block<9> (comb_fb, lfo3, -gains[2]);
    rev.run<10> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<9> (comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (T) (s + in[i]); });
    rev.run<11> (comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    xspan stank {tank.data(), io.size()};
    rev.run<12> (stank);
    rev.run<13> (l, stank.to_const());
    xspan_memdump (r.data(), l);

    xspan er {tmp2.data(), io.size()};
    rev.run<14> (er, in.to_const(), par.character);
    crossfade (l, er, eramt);
    rev.run<15> (l);
    rev.run<16> (l);
    rev.run<17> (l);

    rev.run<18> (er, in.to_const(), par.character);
    crossfade (r, er, eramt);
    rev.run<19> (r);
    rev.run<20> (r);
    rev.run<21> (r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

struct dre2000c : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<dre2000c, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    rev.run<3> (in);
    rev.run<4> (in);
    rev.run<5> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18> (0.25f + dec2 * 20.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.88f + upar.lf_amt * 0.1f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.3f + dec2 * 0.25f + upar.hf_amt * 0.4499f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<6> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<8> (comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block<9> (comb_fb, lfo2, gains[1]);
    rev.run<10> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<9> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<11> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<12> (comb_fb, lfo3, gains[2]);
    rev.run<13> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<12> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<14> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<15> (comb_fb, lfo4, gains[3]);
    rev.run<16> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<15> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run<18> (l);
    rev.run<19> (l);
    rev.run<20, T> (lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt);

    rev.run<21> (r);
    rev.run<22> (r);
    rev.run<23> (in, par.character);
    crossfade (r, in, eramt);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct dre2000d : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<dre2000d, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    rev.run<3> (in);
    rev.run<4> (in);
    rev.run<5> (in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18> (0.25f + dec2 * 5.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.80f + upar.lf_amt * 0.199f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    T ghi = load_float<T> (0.4f + dec2 * 0.1f + upar.hf_amt * 0.49f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<6> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<8> (comb_fb.to_const(), overwrite, r, l);

    rev.fetch_block<9> (comb_fb, blank, gains[1]);
    rev.run<10> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<9> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<11> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<12> (comb_fb, blank, gains[2]);
    rev.run<13> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<12> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<14> (comb_fb.to_const(), add_to, r, l);

    rev.fetch_block<15> (comb_fb, lfo4, gains[3]);
    rev.run<16> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<15> (comb_fb, comb_fb.to_const(), in.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (T) ((c - one<T>() * 0.5_r) * 2_r);
      c        = (T) ((one<T>() - c * c) * 0.4_r);
      eramt[i] = c;
    });

    rev.run<18> (l, lfo2);
    rev.run<19> (l);
    rev.run<20, T> (lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt);

    rev.run<21> (r, lfo3);
    rev.run<22> (r);
    rev.run<23> (in, par.character);
    crossfade (r, in, eramt);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
struct rev5_l_hall : private algorithm {
  //----------------------------------------------------------------------------
  using engine = detail::lofiverb::engine<rev5_l_hall, max_block_size>;
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
  template <class T>
  static void process_block (
    engine&                      rev,
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

    rev.run<0> (in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<T> (lfo_obj);
      lfo1[i]  = (T) (T {lfo[0]} * par.mod[i]);
      lfo2[i]  = (T) (T {lfo[1]} * par.mod[i]);
      lfo3[i]  = (T) (T {lfo[2]} * par.mod[i]);
      lfo4[i]  = (T) (T {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    rev.run<1> (in);
    rev.run<2> (in);
    xspan_memdump (main.data(), in);
    xspan_memdump (sub.data(), in);
    rev.run<3> (main, [&] (auto v, uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      tmp2[i] = (T) (0.35_r + 0.7_r * clamp (c, 0._r, 0.4999_r)); // cg2
      return v * 0.25_r;
    });
    rev.run<4> (main, blank, tmp1);
    rev.run<5> (main, blank, tmp2);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = rev.get_gain_for_rt60<T, 6, 10, 14, 18, 22> (1.f + dec2 * 10.f, srate);
    T flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    T glo = load_float<T> (0.9f + upar.lf_amt * 0.09f);
    T fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.1f);
    T ghi = load_float<T> (0.75f + upar.hf_amt * 0.2f);

    xspan comb_fb {tmp1.data(), io.size()};
    rev.fetch_block<6> (comb_fb, lfo1, gains[0]);
    rev.run<7> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.85_r));
    rev.push<6> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<8> (comb_fb.to_const(), overwrite, l, r);

    rev.fetch_block<9> (comb_fb, lfo3, gains[1]);
    rev.run<10> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.9_r));
    rev.push<9> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<11> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<12> (comb_fb, blank, gains[2]);
    rev.run<13> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.95_r));
    rev.push<12> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<14> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<15> (comb_fb, blank, gains[3]);
    rev.run<16> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.97_r));
    rev.push<15> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<17> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<18> (comb_fb, blank, gains[4]);
    rev.run<19> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<18> (comb_fb, comb_fb.to_const(), main.to_const());
    rev.run<20> (comb_fb.to_const(), add_to, l, r);

    rev.run<21> (sub, par.character);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.5_r + 0.2_r * c); // cg3
      // tmp2 still not overwritten...
      // tmp2[i] = (T) (0.35_r + 0.999_r * clamp (c, 0._r, 0.333_r)); //
      // cg2
    }
    rev.run<22> (main, tmp2);
    rev.run<23> (main, tmp1);
    gains
      = rev.get_gain_for_rt60<T, 29, 33, 37, 41, 45> (1.f + dec2 * 6.f, srate);
    flo = load_float<T> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    glo = load_float<T> (0.9f + upar.lf_amt * 0.04f);
    fhi = load_float<T> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    ghi = load_float<T> (0.7f + upar.hf_amt * 0.25f);

    // sub  reverb
    rev.fetch_block<24> (comb_fb, lfo2, gains[0]);
    rev.run<25> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.9_r));
    rev.push<24> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<26> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<27> (comb_fb, lfo4, gains[1]);
    rev.run<28> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.93_r));
    rev.push<27> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<29> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<30> (comb_fb, blank, gains[2]);
    rev.run<31> (comb_fb, flo, glo, fhi, one<T>(), (T) (ghi * 0.97_r));
    rev.push<30> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<32> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<33> (comb_fb, blank, gains[3]);
    rev.run<34> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<33> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<35> (comb_fb.to_const(), add_to, l, r);

    rev.fetch_block<36> (comb_fb, blank, gains[4]);
    rev.run<37> (comb_fb, flo, glo, fhi, one<T>(), ghi);
    rev.push<36> (comb_fb, comb_fb.to_const(), sub.to_const());
    rev.run<38> (comb_fb.to_const(), add_to, l, r);

    rev.run<39> (l, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo1[i]);
    });
    rev.run<40> (l, blank, [&] (uint i) {
      auto c  = par.character[i];
      tmp1[i] = (T) (0.325_r + 0.999_r * clamp (c, 0._r, 0.333_r)); // cg1
      return (T) (tmp1[i] + 0.05_r * lfo2[i]);
    });
    rev.run<41> (r, blank, [&] (uint i) {
      return (T) (0.65_r + 0.1_r * lfo3[i]);
    });
    rev.run<42> (r, blank, [&] (uint i) {
      return (T) (tmp1[i] + 0.05_r * lfo4[i]);
    });

    rev.run<43> (in.to_const(), add_to, l, r); // ER
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
};

}}} // namespace artv::detail::lofiverb
