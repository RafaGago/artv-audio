#pragma once

#include <type_traits>

#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/misc/primes_table.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
template <delay::data_type Dt>
class abyss : public algorithm {
private:
  using engine = detail::lofiverb::algo_engine<abyss<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (10500, 16800, 21000, 25200, 31500, 40320, 50400);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  void reset (xspan<u8> mem, float t_spl)
  {
    _eng.reset_memory (mem);
    _lfo.reset();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (1000.f),
      vec_set<f32_x2> (0.42f),
      vec_set<f32_x2> (3.85f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f_er   = 0.3f + mod * 0.3f;
    auto f_late = 0.1f + mod * 1.2f;
    _lfo.set_freq (f32_x4 {f_er, f_er, f_late, f_late}, t_spl);
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
      make_ap (1787, 0.5, 261), // 7
      make_block_delay (max_block_size), // 8 to allow block processing
      // loop1
      make_ap (977, 0.5 /*overridden*/, 51), // 9
      make_delay (2819), // 10
      make_crossover2(), // 11
      make_ap (863, -0.5 /*overridden*/), // 12
      make_delay (1021), // 13
      make_ap (1453, 0.618), // 14
      make_block_delay (
        787), // 15 delay (allows block processing) (> blocksz + 1)
      // loop2
      make_ap (947, 0.5 /*overridden*/, 67), // 16
      make_delay (3191), // 17
      make_crossover2(), // 18
      make_ap (887, -0.5 /*overridden*/), // 19
      make_delay (1049), // 20
      make_ap (1367, 0.618), // 21
      make_hp (0.98), // 22
      make_block_delay (
        647)); // 23 delay (allows block processing) (> blocksz + 1)
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<sample>;
    using arr_fb = fb_block_arr<sample>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (io[i][0] + io[i][1]) * 0.5_r;
      auto mod   = 0.25_r + (1_r - par.mod[i]) * 0.75_r;
      // ER + late lfo
      auto lfo = tick_lfo<sample> (_lfo);
      lfo1[i]  = sample {lfo[0]};
      lfo2[i]  = sample {lfo[1]} * (0.5_r + par.character[i] * 0.5_r);
      lfo3[i]  = sample {lfo[2]} * mod;
      lfo4[i]  = sample {lfo[3]} * mod;
      // decay fixup
      auto decay   = 1_r - par.decay[i];
      decay        = 1_r - decay * decay;
      par.decay[i] = 0.6_r + decay * 0.38_r;
    }
    // damp -----------------------------------
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.55f + upar.lf_amt * 0.4f);
    sample fhi = load_float<sample> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    sample ghi = load_float<sample> (0.6f + upar.hf_amt * 0.25f);

    _eng.run (sl<0> {}, late_in);
    // diffusion -----------------------------
    _eng.run (sl<1> {}, late_in);
    _eng.run (sl<2> {}, late_in);
    _eng.run (sl<3> {}, late_in);
    _eng.run (sl<4> {}, late_in);

    // ER (first reverb, not exactly ER at the end...) -----------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    _eng.fetch_block (sl<8> {}, er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = late_in[i] * 0.5_r + er2[i];
    });
    er2.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<5> {}, er1, xspan {lfo2});
    _eng.run (sl<6> {}, er1, flo, glo, fhi, 1_r, ghi);
    xspan_memcpy (er1b, er1);
    span_mul (er1b, par.decay);
    _eng.run (sl<7> {}, er1b, xspan {lfo1});
    _eng.push (sl<8> {}, er1b.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto e_l = -er1[i] * 0.66_r - er2[i] * 0.34_r;
      auto e_r = -er1[i] * 0.66_r + er2[i] * 0.34_r;
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
    _eng.fetch_block (sl<15> {}, l, 1); // feedback, fetching block + 1 samples
    _eng.fetch_block (sl<23> {}, r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = loopsig * (1_r - er_amt) + er_sig * er_amt;
      g[i]         = 0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r);
    }
    r.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<9> {}, late, xspan {lfo3}, g);
    _eng.run (sl<10> {}, late);
    _eng.run (sl<11> {}, late, flo, glo, fhi, 1_r, ghi);
    span_mul (late, par.decay);
    _eng.run (sl<12> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<13> {}, late);
    span_mul (late, par.decay);
    _eng.run (sl<14> {}, late);
    _eng.push (sl<15> {}, late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = loopsig * (1_r - er_amt) + er_sig * er_amt;
    }
    l.cut_head (1); // drop feedback sample from previous block
    _eng.run (sl<16> {}, late, xspan {lfo4}, g);
    _eng.run (sl<17> {}, late);
    _eng.run (sl<18> {}, late, flo, glo, fhi, 1_r, ghi);
    span_mul (late, par.decay);
    _eng.run (sl<19> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<20> {}, late);
    span_mul (late, par.decay);
    _eng.run (sl<21> {}, late);
    _eng.run (sl<22> {}, late);
    _eng.push (sl<23> {}, late.to_const()); // feedback point

    // Mixdown
    auto v1 = make_array (&l[0], &er1[0]);
    auto v2 = make_array<sample const*> (&r[0], &er2[0]);
    ep_crossfade<sample> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return fastgrowth (v * 0.06_r);
    });

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)));
    });
  }
  //----------------------------------------------------------------------------
private:
  engine                              _eng;
  lfo<4>                              _lfo;
  part_class_array<andy::svf, f32_x2> _eq;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class plate1 : public algorithm {
private:
  using engine = detail::lofiverb::algo_engine<plate1<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (10500, 16800, 21000, 25200, 31500, 40320, 50400);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  void reset (xspan<u8> mem, float t_spl)
  {
    _eng.reset_memory (mem);
    _lfo.reset();
    _lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (300.f),
      vec_set<f32_x2> (0.55f),
      vec_set<f32_x2> (1.75f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 0.43f + mod * 0.2f;
    _lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
  }
  //----------------------------------------------------------------------------
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
      make_lp (0.1f), // 12
      make_hp (0.96), // 13
      make_ap (23, 0.8), // 14
      make_ap (130, 0.8), // 15
      make_ap (217, 0.8), // 16

      // s conditioning-diffusors
      make_lp (0.1f), // 17
      make_hp (0.96f), // 18
      make_ap (22, 0.8), // 19
      make_ap (131, 0.8), // 20
      make_ap (219, 0.8), // 21

      // block a iteration 1
      make_ap (153, 0.2), // 22
      make_ap (89, -0.04), // 23
      make_ap (60, 0.04), // 24
      make_delay (201), // 25

      // block b iteration 1
      make_delay (185), // 26

      // block c iteration 1
      make_ap (149, 0.2), // 27
      make_ap (83, 0.04), // 28
      make_ap (59, -0.04), // 29
      make_delay (225), // 30

      // block d iteration 1
      make_ap (167, -0.2), // 31
      make_ap (97, -0.04), // 32
      make_ap (67, 0.04), // 33
      make_delay (221), // 34

      // crossovers
      make_crossover2(), // 35
      make_crossover2(), // 36
      make_crossover2(), // 37
      make_crossover2(), // 38

      // block a iteration 2
      make_ap (119, 0.6), // 39
      make_ap (67, 0.6), // 49
      make_ap (47, -0.6), // 41
      make_block_delay (181), // feedback point // 42

      // block b iteration 2
      make_ap (114, -0.1), // nested 3x // 43
      make_ap (66, -0.04), // 44
      make_ap (47, -0.04), // 45
      make_ap (9), // 46
      make_delay (50, 11), // 47
      make_block_delay (185 - 50), // feedback point // 48

      // block c iteration 2
      make_ap (116, -0.1), // nested 3x // 49
      make_ap (65, -0.04), // 50
      make_ap (46, -0.04), // 51
      make_ap (3), //  52
      make_delay (50, 11), // 53
      make_block_delay (178 - 50), // feedback point // 54

      // block d iteration 2
      make_ap (121, -0.1), // nested 3x // 55
      make_ap (69, 0.04), // 56
      make_ap (47, 0.04), // 57
      make_block_delay (175) // feedback point // 58
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<sample>;
    using arr_fb = fb_block_arr<sample>;

    arr_fb a_mem, d_mem;
    arr    b_mem, c_mem, tmp1, tmp2, tmp3;
    arr    lfo1, lfo2;

    xspan a {a_mem.data(), io.size() + 1};
    xspan b {b_mem.data(), io.size()};
    xspan c {c_mem.data(), io.size()};
    xspan d {d_mem.data(), io.size() + 1};

    // fetch a block from the output (a and d)
    _eng.fetch_block (sl<42> {}, a, 1); // feedback, fetching block + 1 samples
    _eng.fetch_block (sl<58> {}, d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp1[i] = 0.4_r + 0.3_r * par.character[i];
      tmp2[i] = -0.6_r * par.character[i];
      tmp3[i] = 0.65_r * par.character[i];
      // decay fixup
      auto d       = 1_r - par.decay[i];
      par.decay[i] = 0.985_r - d * d * 0.21_r;
      // lfo
      auto lfo = tick_lfo<sample> (_lfo);
      lfo1[i]  = sample {lfo[0]};
      lfo2[i]  = sample {lfo[1]};
      // lfo3 = -lfo1, lfo4 = -lfo2
    }

    // output preconditioning
    _eng.run (sl<0> {}, b, blank, tmp1);
    _eng.run (sl<1> {}, b, blank, tmp2);
    _eng.run (sl<2> {}, b, blank, tmp2);

    _eng.run (sl<3> {}, c, blank, tmp1);
    _eng.run (sl<4> {}, c, blank, tmp2);
    _eng.run (sl<5> {}, c, blank, tmp2);

    // tmp1 and tmp2 free
    xspan cho1 {tmp1.data(), io.size()};
    xspan cho2 {tmp2.data(), io.size()};

    _eng.run (sl<6> {}, cho1, b.to_const(), lfo1, tmp3);
    _eng.run (sl<7> {}, cho1);
    _eng.run (
      sl<8> {}, cho2, c.to_const(), [&] (uint i) { return -lfo1[i]; }, tmp3);
    _eng.run (sl<9> {}, cho2);
    auto v1 = make_array (&b[0], &c[0]);
    auto v2 = make_array<sample const*> (&cho1[0], &cho2[0]);
    ep_crossfade<sample> (v1, v2, &par.mod[0], io.size(), [=] (auto v) {
      return fastgrowth (v * 0.35_r);
    });

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      tmp3[i] = 0.6_r * par.character[i];
    }

    _eng.run (sl<10> {}, b, blank, tmp3);
    _eng.run (sl<11> {}, c, blank, tmp3); // tmp3 free

    xspan i1 {tmp1.data(), io.size()};
    xspan i2 {tmp2.data(), io.size()};
    // output write + preparations
    span_visit (io, [&] (auto& spls, uint i) {
      i1[i]   = (spls[0] + spls[1]) * 0.25_r; // m
      i2[i]   = (spls[0] - spls[1]) * 0.25_r; // s
      spls[0] = b[i];
      spls[1] = c[i];
      b[i]    = 0.4_r + 0.1_r * par.character[i]; // character 1 on b
      c[i]    = 0.8_r * par.character[i]; // character 2 on c
    });

    // input preconditioning
    _eng.run (sl<12> {}, i1);
    _eng.run (sl<13> {}, i1);
    _eng.run (sl<14> {}, i1, blank, b);
    _eng.run (sl<15> {}, i1, blank, c);
    _eng.run (sl<16> {}, i1, blank, c);

    _eng.run (sl<17> {}, i2);
    _eng.run (sl<18> {}, i2);
    _eng.run (sl<19> {}, i2, blank, b);
    _eng.run (sl<20> {}, i2, blank, c);
    _eng.run (sl<21> {}, i2, blank, c);

    // fetch b and d feedback (a and d are already ready) 1 block only
    _eng.fetch_block (sl<48> {}, b, 1); // feedback
    _eng.fetch_block (sl<54> {}, c, 1); // feedback

    span_mul (a, par.decay);
    span_mul (b, par.decay);
    span_mul (c, par.decay);
    span_mul (d, par.decay);

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      b[i] += i1[i] * 0.5_r;
      d[i] += i2[i] * 0.5_r;
    }
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // channel a block 1
    _eng.run (sl<22, 23, 24> {}, a);
    _eng.run (sl<25> {}, a);

    // channel b block 1
    _eng.run (sl<26> {}, b);

    // channel c block 1
    _eng.run (sl<27, 28, 29> {}, c);
    _eng.run (sl<30> {}, c);

    // channel d block 1
    _eng.run (sl<31, 32, 33> {}, d);
    _eng.run (sl<34> {}, d);

    // quantized decay
    span_mul (a, par.decay);
    span_mul (b, par.decay);
    span_mul (c, par.decay);
    span_mul (d, par.decay);

    // block 2
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // swaps.
    a = xspan {d_mem.data(), io.size()};
    b = xspan {c_mem.data(), io.size()};
    c = xspan {a_mem.data(), io.size()};
    d = xspan {b_mem.data(), io.size()};

    // damp -----------------------------------
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.9f + upar.lf_amt * 0.09f);
    sample fhi = load_float<sample> (0.61f - upar.hf_amt * upar.hf_amt * 0.12f);
    sample ghi = load_float<sample> (0.6f + upar.hf_amt * 0.28f);

    _eng.run (sl<35> {}, a, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<36> {}, b, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<37> {}, c, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<38> {}, d, flo, glo, fhi, 1_r, ghi);

    // channel a block 2
    _eng.run (sl<39> {}, a);
    _eng.run (sl<40> {}, a);
    _eng.run (sl<41> {}, a);
    _eng.push (sl<42> {}, a.to_const());

    // channel b block 2
    _eng.run (sl<43, 44, 45> {}, b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      sample m = fastgrowth (par.mod[i]);
      tmp1[i]  = 0.4_r + 0.4_r * lfo2[i] * m;
      tmp2[i]  = 0.6_r + 0.4_r * -lfo2[i] * m; // for block c
    }
    _eng.run (sl<46> {}, b, tmp1);
    _eng.run (sl<47> {}, b, tmp1);
    _eng.push (sl<48> {}, b.to_const());

    // channel c block 2
    _eng.run (sl<49, 50, 51> {}, c);
    _eng.run (sl<52> {}, c, tmp2);
    _eng.run (sl<53> {}, c, tmp2);
    _eng.push (sl<54> {}, c.to_const());

    // channel d block 2
    _eng.run (sl<55, 56, 57> {}, d);
    _eng.push (sl<58> {}, d.to_const());
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)));
    });
  }
  //----------------------------------------------------------------------------
private:
  engine                              _eng;
  lfo<4>                              _lfo;
  part_class_array<andy::svf, f32_x2> _eq;
};
//------------------------------------------------------------------------------
static constexpr std::array<u8, 32> get_ambience_ffwd_table()
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
template <delay::data_type Dt>
class ambience : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<ambience<Dt>, Dt, max_block_size>;
  static constexpr std::array<u8, 32> ffwd_table {get_ambience_ffwd_table()};

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (10500, 16800, 21000, 23400, 33600, 42000, 50400);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  void reset (xspan<u8> mem, float t_spl)
  {
    _eng.reset_memory (mem);
    _lfo.reset();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 1.f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (400.f),
      vec_set<f32_x2> (0.38f),
      vec_set<f32_x2> (6.15f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f = 0.43f + mod * 0.2f;
    _lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
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
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr = block_arr<sample>;

    arr   l_m, r_m;
    xspan l {l_m.data(), io.size()};
    xspan r {r_m.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      l[i] = spl[0] * 0.25_r;
      r[i] = spl[1] * 0.25_r;
    });

    sample lo = load_float<sample> (0.75f - upar.hf_amt * 0.2f);
    sample hi = load_float<sample> (0.9f + upar.lf_amt * 0.075f);

    _eng.run (sl<0> {}, l, lo);
    _eng.run (sl<1> {}, l, hi);
    _eng.run (sl<2> {}, r, lo);
    _eng.run (sl<3> {}, r, hi);

    // Mystran's feedforward diffusor
    mp11::mp_for_each<mp11::mp_from_sequence<std::make_index_sequence<16>>> (
      [&] (auto idx) {
        _eng.run (sl<idx + 4> {}, (idx & 1) ? r : l, [&] (uint i) {
          return ffwd_table[idx + (uint) (par.decay[i] * 16_r)];
        });
        hadamard2 (make_array (l.data(), r.data()), l.size());
      });

    arr   m_m, ltdec, t1, t2, t3;
    xspan m {m_m.data(), io.size()};

    _eng.fetch (sl<20> {}, xspan {t1.data(), io.size()}, 33);
    _eng.fetch (sl<20> {}, xspan {t2.data(), io.size()}, 157);

    span_visit (m, [&] (auto& spl, uint i) {
      spl      = (l[i] + r[i]) * 0.25_r;
      auto inv = 1_r - par.decay[i];
      l[i] += t1[i] * inv;
      r[i] -= t2[i] * inv;
      ltdec[i] = par.decay[i] * par.decay[i];
      t1[i]    = 0.39_r * ltdec[i];
      t2[i]    = -0.28_r * ltdec[i];
      t3[i]    = 0.19_r * ltdec[i];
    });
    // triple nested AP with crossover
    _eng.run (
      sl<20, 21, 22, 23> {},
      m,
      par.character,
      t1,
      par.character,
      t2,
      load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f),
      load_float<sample> (0.55f - upar.hf_amt * upar.hf_amt * 0.15f),
      load_float<sample> (0.85f + upar.lf_amt * 0.13f),
      1_r,
      load_float<sample> (0.35f + upar.hf_amt * 0.14f),
      par.character,
      t3);

    span_visit (xspan {ltdec.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = 0.3_r + v * 0.3_r;
      t2[i] = -0.25_r - v * 0.3_r;
      t3[i] = 0.25_r + v * 0.3_r;
    });
    _eng.run (sl<24> {}, m, blank, t1);
    _eng.run (sl<25> {}, m, blank, t2);
    _eng.run (sl<26> {}, m, blank, t3);

    span_visit (xspan {par.mod.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = sample {tick_lfo<sample> (_lfo)[0]} * v;
      t2[i] = 0.365_r + v * 0.2_r;
    });
    _eng.run (sl<27> {}, xspan {t3.data(), io.size()}, m.to_const(), t1, t2);
    span_visit (io, [&] (auto& spl, uint i) {
      spl[0] = l[i] + par.decay[i] * t3[i];
      spl[1] = r[i] - par.decay[i] * m[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)));
    });
  }
  //----------------------------------------------------------------------------
private:
  engine                              _eng;
  lfo<4>                              _lfo;
  part_class_array<andy::svf, f32_x2> _eq;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class arena : public algorithm {
private:
  using engine = detail::lofiverb::algo_engine<arena<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (10500, 16800, 21000, 23400, 33600, 42000, 50400);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  void reset (xspan<u8> mem, float t_spl)
  {
    _eng.reset_memory (mem);
    _lfo.reset();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (360.f),
      vec_set<f32_x2> (0.45f),
      vec_set<f32_x2> (2.6f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.73f - mod * 0.63f;
    auto f2 = 1.53f - mod * 0.43f;
    _lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp (0.26f), // 0
      make_hp (0.975f), // 1
      make_lp (0.26f), // 2
      make_hp (0.975f), // 3

      make_block_delay (1447), // feedback point -> h1 4

      make_ap (376, 0, 27), // 5
      make_delay (1007), // 6 -> a1

      make_crossover2(), // 7
      make_ap (363), // 8
      make_delay (1107), // 9 -> b1

      make_ap (414, 0, 23), // 10
      make_delay (1207), // 11 -> c1

      make_crossover2(), // 12
      make_ap (477), // 13
      make_delay (1307), // 14 -> d1

      make_ap (420, 0, 21), // 15
      make_delay (1407), // 16 -> e1

      make_crossover2(), // 17
      make_ap (252), // 18
      make_delay (1507), // 19 -> f1

      make_ap (413, 0, 22), // 20
      make_delay (1607), // 21 -> g1

      make_crossover2(), // 22
      make_ap (813) // 23
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<sample>;
    using arr_fb = fb_block_arr<sample>;

    arr loop_mem;
    arr l, r, l_in, r_in, r_cp, m, s, k1, k2, tmp_mem, loop2_mem;
    arr lfo1, lfo2;

    xspan loop {loop_mem.data(), io.size()};
    xspan loop2 {loop2_mem.data(), io.size()};
    xspan tmp {tmp_mem.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      l_in[i] = spl[0] * 0.25_r;
      r_in[i] = spl[1] * 0.25_r;
    });

    _eng.run (sl<0> {}, xspan {l_in.data(), io.size()});
    _eng.run (sl<1> {}, xspan {l_in.data(), io.size()});
    _eng.run (sl<2> {}, xspan {r_in.data(), io.size()});
    _eng.run (sl<3> {}, xspan {r_in.data(), io.size()});
    _eng.fetch_block (sl<4> {}, loop, 1); // feedback signal

    span_visit (loop, [&] (auto& v, uint i) {
      m[i]         = (l_in[i] + r_in[i]) * 0.5_r;
      s[i]         = (l_in[i] - r_in[i]) * 0.5_r;
      k1[i]        = 0.4_r + par.character[i] * 0.2_r;
      k2[i]        = 0.4_r + par.character[i] * 0.15_r;
      auto lfo     = tick_lfo<sample> (_lfo);
      lfo1[i]      = sample {lfo[0]} * par.mod[i];
      lfo2[i]      = sample {lfo[2]} * par.mod[i];
      v            = v + l_in[i] * 0.75_r + m[i] * 0.25_r;
      par.decay[i] = 0.05_r + fastgrowth (par.decay[i]) * 0.95_r;
      loop[i] *= par.decay[i];
    });

    auto flo = load_float<sample> (0.90 + upar.lf_amt * upar.lf_amt * 0.05);
    auto glo = load_float<sample> (0.75 + upar.lf_amt * 0.23);
    auto fhi = load_float<sample> (0.9 - upar.hf_amt * upar.hf_amt * 0.4);
    auto ghi = load_float<sample> (0.7 + upar.hf_amt * 0.23);
    // A
    _eng.run (sl<5> {}, loop, lfo1, k1);
    _eng.fetch (sl<6> {}, tmp, 400 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = loop[i];
      r[i] = v * 0.4_r;
    });
    _eng.run (sl<6> {}, loop);
    // B
    _eng.run (sl<7> {}, loop, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<8> {}, loop, blank, k1);
    _eng.fetch (sl<9> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] += v * 0.4_r;
      r[i] += loop[i];
    });
    _eng.run (sl<9> {}, loop);
    span_mul (loop, par.decay);
    // C
    span_visit (loop, [&] (auto& v, uint i) { v += s[i]; });
    _eng.run (sl<10> {}, loop, lfo2, [&] (uint i) { return -k2[i]; });
    _eng.fetch (sl<11> {}, tmp, 1001 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] -= loop[i];
      r[i] -= v * 0.42_r;
    });
    _eng.run (sl<11> {}, loop);
    // D
    _eng.run (sl<12> {}, loop, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<13> {}, loop, blank, k2);
    _eng.fetch (sl<14> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] += v * 0.42_r;
      r[i] -= loop[i] * 0.49_r;
    });
    _eng.run (sl<14> {}, loop);
    span_mul (loop, par.decay);
    // E
    span_visit (loop, [&] (auto& v, uint i) {
      v = v + r_in[i] * 0.75_r + m[i] * 0.25_r;
    });
    _eng.run (
      sl<15> {}, loop, [&] (uint i) { return -lfo1[i]; }, k1);
    _eng.fetch (sl<16> {}, tmp, 801 - 37 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] += loop[i] * 0.49_r;
      r[i] -= v;
    });
    _eng.run (sl<16> {}, loop);
    // F
    _eng.run (sl<17> {}, loop, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<18> {}, loop, blank, k1);
    _eng.fetch (sl<19> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] += v;
      r[i] += loop[i];
    });
    _eng.run (sl<19> {}, loop);
    span_mul (loop, par.decay);
    // G
    span_visit (loop, [&] (auto& v, uint i) { v += s[i]; });
    _eng.run (
      sl<20> {}, loop, [&] (uint i) { return -lfo2[i]; }, k2);
    _eng.fetch (sl<21> {}, tmp, 1001 - 27 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] -= loop[i];
      r[i] -= v;
    });
    _eng.run (sl<21> {}, loop);
    // H
    _eng.run (sl<22> {}, loop, flo, glo, fhi, 1_r, ghi);
    _eng.run (sl<23> {}, loop, blank, k2);
    _eng.fetch (sl<4> {}, tmp, 1001 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] += loop[i] * 0.2_r;
      r[i] += v * 0.2_r;
    });
    _eng.push (sl<4> {}, loop.to_const());

    span_visit (io, [&] (auto& spl, uint i) {
      constexpr auto sqrt8_recip = 0.3535534_r; // 1/sqrt(8), equal power mixing
      spl[0]                     = l[i] * sqrt8_recip;
      spl[1]                     = r[i] * sqrt8_recip;
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)));
    });
  }
  //----------------------------------------------------------------------------
private:
  engine                              _eng;
  lfo<4>                              _lfo;
  part_class_array<andy::svf, f32_x2> _eq;
};
//------------------------------------------------------------------------------
}}} // namespace artv::detail::lofiverb
