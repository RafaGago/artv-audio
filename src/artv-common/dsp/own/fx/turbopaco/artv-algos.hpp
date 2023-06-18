#pragma once

#include <type_traits>

#include "artv-common/dsp/own/fx/turbopaco/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/misc/primes_table.hpp"

namespace artv { namespace detail { namespace tpaco {

//------------------------------------------------------------------------------
template <delay::data_type Dt>
class abyss : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<abyss<Dt>, Dt, max_block_size>;

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
  using engine = detail::tpaco::algo_engine<plate1<Dt>, Dt, max_block_size>;

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
  using engine = detail::tpaco::algo_engine<ambience<Dt>, Dt, max_block_size>;
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
      make_ap (913, 0.f, 413, interpolation::zero_order_hold), // 20
      make_ap (491, 0.f, 213, interpolation::zero_order_hold), // 21
      make_crossover2(), // 22
      make_ap (491, 0.f, 213, interpolation::zero_order_hold), // 23
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
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 2.f);
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
class palace : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<palace<Dt>, Dt, max_block_size>;

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
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 4.f);
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
  using engine = detail::tpaco::algo_engine<arena<Dt>, Dt, max_block_size>;

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
      vec_set<f32_x2> (305.f),
      vec_set<f32_x2> (0.64f),
      vec_set<f32_x2> (3.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.73f - mod * 0.63f;
    auto f2 = 1.53f - mod * 0.43f;
    _lfo.set_freq (f32_x4 {f1, f1, f2, f2}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.0625f;
    return make_array<stage_data> (
      make_lp (0.16), // 0
      make_hp (0.99), // 1
      make_ap (123, 0.), // 3
      make_lp (0.16), // 4
      make_hp (0.99), // 5
      make_ap (124, 0.), // 6

      make_comb (1821, 0.f), // 6
      make_crossover2(), // 7
      make_ap (912, 0.), // 8
      make_ap (736, 0.), // 9
      make_ap (515, 0.1), // 10
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        114, g, 0,
        112, -g, 0,
        342, g, 13,
        377, g, 16,
        412, g, 77,
        410, -g, 88,
        1112 - 46, g, 89,
        1112 - 45, g, 73
        ), // 11
      // clang-format on
      make_comb (1827 + 387, 0.f), // 12
      make_crossover2(), // 13
      make_ap (772, 0.), // 14
      make_ap (666, 0.), // 15
      make_ap (535, 0.1), // 16
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        232, -g, 0,
        230, g, 0,
        666, g, -77,
        667, g, -88,
        812, g, -21,
        831, -g,-15,
        1312, -g, -89,
        1312, g, -73
        ), // 17
      // clang-format on
      make_comb (1827 + 38 + 56, 0.f), // 18
      make_crossover2(), // 19
      make_ap (672, 0.), // 20
      make_ap (526, 0.), // 21
      make_ap (425, 0.1), // 22
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        192, -g, 0,
        196, g, 0,
        376, -g, 77,
        377, g, 88,
        777, -g, 89,
        776, g, 73,
        1472, g, 17,
        1471, g, 23
        ), // 23
      // clang-format on
      make_comb (1827 + 38 + 56 + 37, 0.f), // 24
      make_crossover2(), // 25
      make_ap (572, 0.), // 26
      make_ap (426, 0.), // 27
      make_ap (325, 0.1), // 28
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        377, -g, -22,
        376, g, -17,
        642, g, -89,
        677, -g, -73,
        992, g, -77,
        992, g, -98,
        1612, -g, 0,
        1612, g, 0
        ) // 29
      // clang-format on
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
    uint                         srate)
  {
    block_arr<sample> l_in, r_in, m_in, l, r, lfo1, lfo1m, lfo2, lfo2m, k1, k2,
      k4;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo      = tick_lfo<sample> (_lfo);
      lfo1[i]       = sample {lfo[0]};
      lfo2[i]       = sample {lfo[2]};
      lfo1m[i]      = lfo1[i] * par.mod[i];
      lfo2m[i]      = lfo2[i] * par.mod[i];
      l_in[i]       = io[i][0] * 0.125_r * 0.125_r;
      r_in[i]       = io[i][1] * 0.125_r * 0.125_r;
      k1[i]         = 0.12_r + par.decay[i] * 0.2_r;
      k2[i]         = 0.12_r + par.decay[i] * 0.1_r;
      auto fg_decay = fastgrowth (par.decay[i]);
      k4[i]         = 0.1_r + fg_decay * 0.1_r + par.character[i] * 0.3_r;
    }
    run_cascade<0, 1> (_eng, xspan {l_in.data(), io.size()});
    _eng.run (sl<2> {}, xspan {l_in.data(), io.size()}, blank, k4);
    run_cascade<3, 4> (_eng, xspan {r_in.data(), io.size()});
    _eng.run (sl<5> {}, xspan {r_in.data(), io.size()}, blank, k4);

    float dec   = as_float (par.decay[0]);
    auto  gains = _eng.get_gain_for_rt60 (
      sl<6, 12, 18, 24> {}, 0.35f + dec * dec * 4.75f, srate);
    sample flo  = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo  = load_float<sample> (0.85f + upar.lf_amt * 0.145f);
    sample fhi1 = load_float<sample> (0.9f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample fhi2 = load_float<sample> (0.93f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample ghi  = load_float<sample> (0.5f + upar.hf_amt * 0.5f);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      m_in[i] = (l_in[i] + r_in[i]) * 0.5_r;
    }

    block_arr<sample> combmem;
    xspan             comb {combmem.data(), io.size()};
    auto              in = xspan {m_in.data(), io.size()}.to_const();

    _eng.fetch (sl<6> {}, comb, blank, gains[0]);
    _eng.run (sl<7> {}, comb, flo, glo, fhi1, 1_r, ghi);
    _eng.run (sl<8, 9, 10> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<6> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<11> {},
      comb.to_const(),
      overwrite,
      l,
      r,
      blank,
      blank,
      lfo1,
      lfo1,
      lfo1m,
      lfo1m,
      lfo1m,
      lfo1m);

    _eng.fetch (sl<12> {}, comb, blank, gains[1]);
    _eng.run (sl<13> {}, comb, flo, glo, fhi2, 1_r, ghi);
    _eng.run (sl<14, 15, 16> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<12> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<17> {},
      comb.to_const(),
      add_to,
      l,
      r,
      blank,
      blank,
      lfo1m,
      lfo1m,
      lfo1,
      lfo1,
      lfo1m,
      lfo1m);

    _eng.fetch (sl<18> {}, comb, 615); // ER L
    span_add_with_factor (xspan {l.data(), io.size()}, comb, -0.1_r);

    in = xspan {l_in.data(), io.size()}.to_const();
    _eng.fetch (sl<18> {}, comb, blank, gains[2]);
    _eng.run (sl<19> {}, comb, flo, glo, fhi1, 1_r, ghi);
    _eng.run (sl<20, 21, 22> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<18> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<23> {},
      comb.to_const(),
      add_to,
      l,
      r,
      blank,
      blank,
      lfo2m,
      lfo2m,
      lfo2m,
      lfo2m,
      lfo2,
      lfo2);

    _eng.fetch (sl<24> {}, comb, 873); // ER R
    span_add_with_factor (xspan {r.data(), io.size()}, comb, 0.1_r);

    in = xspan {r_in.data(), io.size()}.to_const();
    _eng.fetch (sl<24> {}, comb, blank, gains[3]);
    _eng.run (sl<25> {}, comb, flo, glo, fhi2, 1_r, ghi);
    _eng.run (sl<26, 27, 28> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<24> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<29> {},
      comb.to_const(),
      add_to,
      l,
      r,
      lfo2,
      lfo2,
      lfo2m,
      lfo2m,
      lfo2m,
      lfo2m,
      blank,
      blank);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 64.f);
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
class hall : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<hall<Dt>, Dt, max_block_size>;

  struct stage {
    uint  d1, d2, dd, tap_l, tap_r;
    float fac_l, fac_r, k1, k2;
  };

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
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.f, 0.f, 0.f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (360.f),
      vec_set<f32_x2> (0.48f),
      vec_set<f32_x2> (3.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.27f - mod * 0.53f;
    auto f2 = 1.51f - mod * 0.43f;
    auto f3 = 1.74f - mod * 0.37f;
    auto f4 = 1.84f - mod * 0.31f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f4}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    // just to avoid copy-paste issues from the prototype
    stage b1 {}, b2 {}, b3 {}, b4 {}, b5 {}, b6 {}, b7 {}, b8 {}, b9 {}, b10 {};
    b1.d1     = 571;
    b1.d2     = 359;
    b1.dd     = 857;
    b1.fac_l  = -0.40780916679832;
    b1.fac_r  = 0.59219083320168;
    b1.k1     = -0.33;
    b1.k2     = 0.14;
    b1.tap_l  = 329;
    b1.tap_r  = 343;
    b2.d1     = 467;
    b2.d2     = 257;
    b2.dd     = 707;
    b2.fac_l  = -0.49429981212465;
    b2.fac_r  = 0.50570018787535;
    b2.k1     = -0.32;
    b2.k2     = 0.14;
    b2.tap_l  = 316;
    b2.tap_r  = 286;
    b3.d1     = 577;
    b3.d2     = 383;
    b3.dd     = 941;
    b3.fac_l  = 0.25340987584866;
    b3.fac_r  = -0.74659012415134;
    b3.k1     = -0.33;
    b3.k2     = 0.14;
    b3.tap_l  = 492;
    b3.tap_r  = 456;
    b4.d1     = 449;
    b4.d2     = 283;
    b4.dd     = 757;
    b4.fac_l  = -0.77724914480402;
    b4.fac_r  = 0.22275085519598;
    b4.k1     = 0.32;
    b4.k2     = 0.14;
    b4.tap_l  = 287;
    b4.tap_r  = 295;
    b5.d1     = 619;
    b5.d2     = 353;
    b5.dd     = 991 + 100;
    b5.fac_l  = -0.4534214228516;
    b5.fac_r  = 0.5465785771484;
    b5.k1     = -0.33;
    b5.k2     = 0.13;
    b5.tap_l  = 353;
    b5.tap_r  = 388;
    b6.d1     = 443;
    b6.d2     = 241;
    b6.dd     = 691 + 100;
    b6.fac_l  = 0.38615093528902;
    b6.fac_r  = -0.61384906471098;
    b6.k1     = 0.32;
    b6.k2     = -0.14;
    b6.tap_l  = 285;
    b6.tap_r  = 264;
    b7.d1     = 599;
    b7.d2     = 373;
    b7.dd     = 919;
    b7.fac_l  = -0.93732491320402;
    b7.fac_r  = 0.06267508679598;
    b7.k1     = -0.33;
    b7.k2     = 0.14;
    b7.tap_l  = 424;
    b7.tap_r  = 459;
    b8.d1     = 467;
    b8.d2     = 269;
    b8.dd     = 709 + 84;
    b8.fac_l  = -0.44557774496395;
    b8.fac_r  = 0.554422;
    b8.k1     = -0.32;
    b8.k2     = 0.14;
    b8.tap_l  = 302;
    b8.tap_r  = 265;
    b9.d1     = 601;
    b9.d2     = 383;
    b9.dd     = 919;
    b9.fac_l  = 0.18701;
    b9.fac_r  = 0.8129;
    b9.k1     = -0.33;
    b9.k2     = -0.14;
    b9.tap_l  = 332;
    b9.tap_r  = 354;
    b10.d1    = 467;
    b10.d2    = 277;
    b10.dd    = 701;
    b10.fac_l = 0.15538;
    b10.fac_r = 0.8446;
    b10.k1    = -0.32;
    b10.k2    = -0.14;
    b10.tap_l = 357;
    b10.tap_r = 360;

    // on the prototype fetch is done after pushing, on here it's done before
    b1.tap_l -= 1;
    b2.tap_l -= 1;
    b3.tap_l -= 1;
    b4.tap_l -= 1;
    b5.tap_l -= 1;
    b6.tap_l -= 1;
    b7.tap_l -= 1;
    b8.tap_l -= 1;
    b9.tap_l -= 1;
    b10.tap_l -= 1;

    b1.tap_r -= 1;
    b2.tap_r -= 1;
    b3.tap_r -= 1;
    b4.tap_r -= 1;
    b5.tap_r -= 1;
    b6.tap_r -= 1;
    b7.tap_r -= 1;
    b8.tap_r -= 1;
    b9.tap_r -= 1;
    b10.tap_r -= 1;

    return make_array<stage_data> (
      make_lp (0.15), // 0
      make_hp (0.99), // 1
      make_ap (73, 0.5), // 2
      make_lp (0.15), // 3
      make_hp (0.99), // 4
      make_ap (77, -0.5), // 5

      make_ap (b1.d1, b1.k1, 23), // 6
      make_ap (b1.d2, b1.k2), // 7
      make_parallel_delay (
        2, b1.dd, b1.tap_l, b1.fac_l, b1.tap_r, b1.fac_r), // 8
      make_crossover2(), // 9

      make_ap (b2.d1, b2.k1), // 10
      make_ap (b2.d2, b2.k2), // 11
      make_parallel_delay (
        2, b2.dd, b2.tap_l, b2.fac_l, b2.tap_r, b2.fac_r), // 12

      make_ap (b3.d1, b3.k1, 17), // 13
      make_ap (b3.d2, b3.k2), // 14
      make_parallel_delay (
        2, b3.dd, b3.tap_l, b3.fac_l, b3.tap_r, b3.fac_r), // 15

      make_ap (b4.d1, b4.k1), // 16
      make_ap (b4.d2, b4.k2), // 17
      make_parallel_delay (
        2, b4.dd, b4.tap_l, b4.fac_l, b4.tap_r, b4.fac_r), // 18

      make_ap (b5.d1, b5.k1, 18), // 19
      make_ap (b5.d2, b5.k2), // 20
      make_parallel_delay (
        2, b5.dd, b5.tap_l, b5.fac_l, b5.tap_r, b5.fac_r), // 21
      make_crossover2(), // 22

      make_ap (b6.d1, b6.k1), // 23
      make_ap (b6.d2, b6.k2), // 24
      make_parallel_delay (
        2, b6.dd, b6.tap_l, b6.fac_l, b6.tap_r, b6.fac_r), // 25

      make_ap (b7.d1, b7.k1, 11), // 26
      make_ap (b7.d2, b7.k2), // 27
      make_parallel_delay (
        2, b7.dd, b7.tap_l, b7.fac_l, b7.tap_r, b7.fac_r), // 28

      make_ap (b8.d1, b8.k1), // 29
      make_ap (b8.d2, b8.k2), // 30
      make_parallel_delay (
        2, b8.dd, b8.tap_l, b8.fac_l, b8.tap_r, b8.fac_r), // 31

      make_ap (b5.d1, b5.k1, -12), // 32
      make_ap (b5.d2, b5.k2), // 33
      make_parallel_delay (
        2, b5.dd, b5.tap_l, b5.fac_l, b5.tap_r, b5.fac_r), // 34
      make_crossover2(), // 35

      make_ap (b6.d1, b6.k1), // 36
      make_ap (b6.d2, b6.k2), // 37
      make_parallel_delay (
        2, b6.dd, b6.tap_l, b6.fac_l, b6.tap_r, b6.fac_r) // 38

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
    uint                         srate)
  {
    block_arr<sample> l_in, r_in, m_in, loopm, l1, r1, l2, r2, lfo1, lfo2, lfo3,
      lfo4;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo     = tick_lfo<sample> (_lfo);
      lfo1[i]      = sample {lfo[0]} * par.mod[i];
      lfo2[i]      = sample {lfo[1]} * par.mod[i];
      lfo3[i]      = sample {lfo[2]} * par.mod[i];
      lfo4[i]      = sample {lfo[3]} * par.mod[i];
      l_in[i]      = io[i][0] * 0.125_r * 0.25_r;
      r_in[i]      = io[i][1] * 0.125_r * 0.25_r;
      par.decay[i] = fastgrowth (par.decay[i]);
      par.decay[i] = 0.4_r + par.decay[i] * 0.59_r;
    }
    run_cascade<0, 2> (_eng, xspan {l_in.data(), io.size()});
    run_cascade<3, 5> (_eng, xspan {r_in.data(), io.size()});

    auto   lf   = upar.lf_amt;
    auto   hf   = upar.hf_amt;
    sample flo  = load_float<sample> (0.8f + lf * lf * 0.1f);
    sample glo  = load_float<sample> (0.75f + lf * 0.245f);
    sample fhi1 = load_float<sample> (0.9f - hf * hf * 0.4f);
    sample fhi2 = load_float<sample> (0.93f - hf * hf * 0.4f);
    sample fhi3 = load_float<sample> (0.95f - hf * hf * 0.4f);
    sample ghi  = load_float<sample> (
      0.85f + hf * 0.1f + as_float (par.decay[0]) * 0.05f);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      m_in[i] = (l_in[i] + r_in[i]) * 0.5_r;
    }

    xspan loop {loopm.data(), io.size()};

    _eng.fetch (sl<38> {}, loop);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<6, 7> {}, loop, lfo1);
    _eng.run (sl<8> {}, loop, loop.to_const(), overwrite, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<9> {}, loop, flo, glo, fhi1, 1_r, ghi);

    _eng.run (sl<10, 11> {}, loop);
    _eng.run (sl<12> {}, loop, loop.to_const(), overwrite, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, l_in);

    _eng.run (sl<13, 14> {}, loop, lfo2);
    _eng.run (sl<15> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);

    _eng.run (sl<16, 17> {}, loop);
    _eng.run (sl<18> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<19, 20> {}, loop, lfo3);
    _eng.run (sl<21> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<22> {}, loop, flo, glo, fhi2, 1_r, ghi);

    _eng.run (sl<23, 24> {}, loop);
    _eng.run (sl<25> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, r_in);

    _eng.run (sl<26, 27> {}, loop, lfo4);
    _eng.run (sl<28> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);

    _eng.run (sl<29, 30> {}, loop);
    _eng.run (sl<31> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<32, 33> {}, loop, lfo1);
    _eng.run (sl<34> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<35> {}, loop, flo, glo, fhi3, 1_r, ghi);

    _eng.run (sl<36, 37> {}, loop);
    _eng.run (sl<38> {}, loop, loop.to_const(), add_to, l2, r1);

    auto v1 = make_array (&l1[0], &r1[0]);
    auto v2 = make_array<sample const*> (&l2[0], &r2[0]);
    ep_crossfade<sample> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return 0.1_r + v * 0.9_r;
    });
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l1[i];
      spls[1] = r1[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 8.f);
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
class broken_hall : public algorithm {
private:
  using engine
    = detail::tpaco::algo_engine<broken_hall<Dt>, Dt, max_block_size>;

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
      vec_set<f32_x2> (405.f),
      vec_set<f32_x2> (0.78f),
      vec_set<f32_x2> (2.5f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.73f - mod * 0.53f;
    auto f2 = 1.51f - mod * 0.43f;
    _lfo.set_freq (f32_x4 {f1, f1, f2, f2}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g1 = 0.75f / 4.f;
    constexpr float g2 = 0.5f / 4.f;
    constexpr float g3 = 0.25f / 4.f;
    constexpr float g4 = 0.125f / 4.f;
    constexpr float sc = 0.906f; // algo sizes are scaled down by this factor

    return make_array<stage_data> (
      make_lp (0.16), // 0
      make_hp (0.99), // 1
      make_lp (0.16), // 2
      make_hp (0.99), // 3

      make_comb (1022 * sc, 0.f), // 4
      make_crossover2(), // 5
      make_ap (781 * sc, 0.), // 6
      make_ap (436 * sc, 0.), // 7
      make_ap (312 * sc, -0.05), // 8
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        97, g1, 0,
        112, g3, 0,
        342, -g1, 48,
        377, g2, -41,
        474 - 46, g4, 89,
        472 - 45, -g4, -73,
        522, g2, 27,
        540, -g3, -28
        ), // 9
      // clang-format on
      make_comb (1053 * sc, 0.f), // 10
      make_crossover2(), // 11
      make_ap (807 * sc, 0.), // 12
      make_ap (445 * sc, 0.), // 13
      make_ap (316 * sc, 0.05), // 14
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        170, g2, 0,
        160, -g1, 0,
        396, -g1, -77,
        397, g3, 88,
        438, g3, -16,
        471, -g2, 14,
        512 - 18, g4, -89,
        513 - 18, -g4, 73
        ), // 15
      // clang-format on
      make_comb (955 * sc, 0.f), // 16
      make_crossover2(), // 17
      make_ap (734 * sc, 0.), // 18
      make_ap (420 * sc, 0.), // 19
      make_ap (282 * sc, 0.05), // 20
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        422, -g1, 0,
        415, -g2, 0,
        646, g1, -57,
        647, -g2, 58,
        673, g4, -7,
        671, -g4, 13,
        858, -g3, -89,
        836, g3, 73
        ), // 21
      // clang-format on
      make_comb (1433 * sc, 0.f), // 22
      make_crossover2(), // 23
      make_ap (1111 * sc, 0.), // 24
      make_ap (635 * sc, 0.), // 25
      make_ap (433 * sc, 0.05), // 26
      // clang-format off
      make_mod_parallel_delay (
        2,
        interpolation::thiran,
        blank,
        267, g2, 9,
        297, -g1, -10,
        776 + 39, g4, 0,
        775 + 39, g4, 0,
        842, g2, 39,
        877, g1, -33,
        1103, -g3, 77,
        1131, -g3, -98
        ), // 27
      // clang-format on
      make_ap (32, 0., (157 * 3) - 32, interpolation::zero_order_hold), // 28
      make_ap (32, 0., (243 * 3) - 32, interpolation::zero_order_hold), // 29
      make_ap (32, 0., (373 * 3) - 32, interpolation::zero_order_hold), // 30
      make_ap (32, 0., (167 * 3) - 32, interpolation::zero_order_hold), // 31
      make_ap (32, 0., (254 * 3) - 32, interpolation::zero_order_hold), // 32
      make_ap (32, 0., (383 * 3) - 32, interpolation::zero_order_hold) // 33
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
    uint                         srate)
  {
    block_arr<sample> l_in, r_in, m_in, l, r, lfo1, lfo1m, lfo2, lfo2m, k1, k2,
      k4;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo         = tick_lfo<sample> (_lfo);
      lfo1[i]          = sample {lfo[0]};
      lfo2[i]          = sample {lfo[2]};
      lfo1m[i]         = lfo1[i] * par.mod[i];
      lfo2m[i]         = lfo2[i] * par.mod[i];
      l_in[i]          = io[i][0] * 0.125_r * 0.25_r;
      r_in[i]          = io[i][1] * 0.125_r * 0.25_r;
      k1[i]            = 0.08_r + par.decay[i] * 0.1_r;
      k2[i]            = 0.08_r + par.decay[i] * 0.08_r;
      auto fg_decay    = fastgrowth (par.decay[i]);
      par.character[i] = 1_r - par.character[i];
      k4[i]            = 0.15_r + fg_decay * 0.15_r + par.character[i] * 0.2_r;
    }
    run_cascade<0, 1> (_eng, xspan {l_in.data(), io.size()});
    run_cascade<2, 3> (_eng, xspan {r_in.data(), io.size()});

    float dec   = as_float (par.decay[0]);
    auto  gains = _eng.get_gain_for_rt60 (
      sl<4, 10, 16, 22> {}, 0.3f + dec * dec * 4.75f, srate);
    sample flo  = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo  = load_float<sample> (0.85f + upar.lf_amt * 0.145f);
    sample fhi1 = load_float<sample> (0.9f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample fhi2 = load_float<sample> (0.93f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample ghi  = load_float<sample> (0.65f + upar.hf_amt * 0.35f);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      m_in[i] = (l_in[i] + r_in[i]) * 0.5_r;
    }

    block_arr<sample> combmem;
    xspan             comb {combmem.data(), io.size()};
    auto              in = xspan {m_in.data(), io.size()}.to_const();

    _eng.fetch (sl<4> {}, comb, blank, gains[0]);
    _eng.run (sl<5> {}, comb, flo, glo, fhi1, 1_r, ghi);
    _eng.run (sl<6, 7, 8> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<4> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<9> {},
      comb.to_const(),
      overwrite,
      l,
      r,
      blank,
      blank,
      lfo1,
      lfo1,
      lfo1m,
      lfo1m,
      lfo1m,
      lfo1m);

    _eng.fetch (sl<10> {}, comb, blank, -gains[1]);
    _eng.run (sl<11> {}, comb, flo, glo, fhi2, 1_r, ghi);
    _eng.run (sl<12, 13, 14> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<10> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<15> {},
      comb.to_const(),
      add_to,
      l,
      r,
      blank,
      blank,
      lfo1m,
      lfo1m,
      lfo1,
      lfo1,
      lfo1m,
      lfo1m);

    in = xspan {l_in.data(), io.size()}.to_const();
    _eng.fetch (sl<16> {}, comb, blank, gains[2]);
    _eng.run (sl<17> {}, comb, flo, glo, fhi1, 1_r, ghi);
    _eng.run (sl<18, 19, 20> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<16> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<21> {},
      comb.to_const(),
      add_to,
      l,
      r,
      blank,
      blank,
      lfo2m,
      lfo2m,
      lfo2,
      lfo2,
      lfo2m,
      lfo2m);
    // TODO: why is it distorting?
    //_eng.fetch (sl<16> {}, comb, 215); // ER L
    // span_add_with_factor (xspan {l.data(), io.size()}, comb, 0.1_r);

    in = xspan {r_in.data(), io.size()}.to_const();
    _eng.fetch (sl<22> {}, comb, blank, gains[3]);
    _eng.run (sl<23> {}, comb, flo, glo, fhi2, 1_r, ghi);
    _eng.run (sl<24, 25, 26> {}, comb, blank, k1, blank, k2);
    _eng.push (sl<22> {}, comb, comb.to_const(), in);
    _eng.run (
      sl<27> {},
      comb.to_const(),
      add_to,
      l,
      r,
      lfo2,
      lfo2,
      blank,
      blank,
      lfo2m,
      lfo2m,
      lfo2m,
      lfo2m);
    // TODO: why is it distorting?
    //_eng.fetch (sl<22> {}, comb, 298); // ER R
    // span_add_with_factor (xspan {r.data(), io.size()}, comb, -0.1_r);

    xspan lspan {l.data(), io.size()};
    xspan rspan {r.data(), io.size()};
    _eng.run (sl<28> {}, lspan, par.character, [&k4] (uint i) {
      return -k4[i];
    });
    run_cascade<29, 30> (_eng, lspan, par.character, k4);
    _eng.run (sl<31> {}, rspan, par.character, [&k4] (uint i) {
      return -k4[i];
    });
    run_cascade<32, 33> (_eng, rspan, par.character, k4);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 16.f);
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
class room : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<room<Dt>, Dt, max_block_size>;

  struct stage {
    uint  d1, d2, dd, tap_l, tap_r;
    float fac_l, fac_r, k1, k2;
  };

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
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.f, 0.f, 0.f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (450.f),
      vec_set<f32_x2> (0.78f),
      vec_set<f32_x2> (3.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.27f - mod * 0.53f;
    auto f2 = 1.51f - mod * 0.43f;
    auto f3 = 1.74f - mod * 0.37f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f3}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    float coeff1 = 0.27f;
    float coeff2 = 0.16f;
    // just to avoid copy-paste issues from the prototype
    stage b1 {}, b2 {}, b3 {}, b4 {}, b5 {}, b6 {}, b7 {}, b8 {}, b9 {}, b10 {},
      b11 {}, b12 {};
    b1.d1     = 281;
    b1.d2     = 199;
    b1.dd     = 487;
    b1.fac_l  = -0.58830350534718;
    b1.fac_r  = 0.41169649465282;
    b1.k1     = coeff1;
    b1.k2     = -coeff2;
    b1.tap_l  = 65;
    b1.tap_r  = 63;
    b2.d1     = 317;
    b2.d2     = 173;
    b2.dd     = 547;
    b2.fac_l  = -0.49055438174181;
    b2.fac_r  = 0.50944561825819;
    b2.k1     = coeff1;
    b2.k2     = coeff2;
    b2.tap_l  = 284;
    b2.tap_r  = 257;
    b3.d1     = 307;
    b3.d2     = 179;
    b3.dd     = 503;
    b3.fac_l  = -0.79811344337606;
    b3.fac_r  = 0.20188655662304;
    b3.k1     = -coeff1;
    b3.k2     = coeff2;
    b3.tap_l  = 115;
    b3.tap_r  = 97;
    b4.d1     = 293;
    b4.d2     = 139;
    b4.dd     = 443;
    b4.fac_l  = 0.17068128943692;
    b4.fac_r  = 0.82931871056308;
    b4.k1     = coeff1;
    b4.k2     = coeff2;
    b4.tap_l  = 113;
    b4.tap_r  = 91;
    b5.d1     = 359;
    b5.d2     = 251;
    b5.dd     = 587;
    b5.fac_l  = 0.11278212701734;
    b5.fac_r  = 0.88721787298266;
    b5.k1     = -coeff1;
    b5.k2     = -coeff2;
    b5.tap_l  = 124;
    b5.tap_r  = 101;
    b6.d1     = 337;
    b6.d2     = 199;
    b6.dd     = 571;
    b6.fac_l  = -0.89145277461304;
    b6.fac_r  = 0.10854722538696;
    b6.k1     = coeff1;
    b6.k2     = -coeff2;
    b6.tap_l  = 273;
    b6.tap_r  = 240;
    b7.d1     = 359;
    b7.d2     = 229;
    b7.dd     = 617;
    b7.fac_l  = -0.80255597312994;
    b7.fac_r  = 0.19744402687006;
    b7.k1     = -coeff1;
    b7.k2     = -coeff2;
    b7.tap_l  = 146;
    b7.tap_r  = 172;
    b8.d1     = 373;
    b8.d2     = 241;
    b8.dd     = 613;
    b8.fac_l  = -0.85686479435695;
    b8.fac_r  = 0.14313520564305;
    b8.k1     = -coeff1;
    b8.k2     = -coeff2;
    b8.tap_l  = 322;
    b8.tap_r  = 294;
    b9.d1     = 433;
    b9.d2     = 281;
    b9.dd     = 691;
    b9.fac_l  = -0.87407659806173;
    b9.fac_r  = 0.12592340193827;
    b9.k1     = -coeff1;
    b9.k2     = coeff2;
    b9.tap_l  = 355;
    b9.tap_r  = 345;
    b10.d1    = 449;
    b10.d2    = 277;
    b10.dd    = 733;
    b10.fac_l = -0.20988618447676;
    b10.fac_r = 0.79011381552324;
    b10.k1    = -coeff1;
    b10.k2    = coeff2;
    b10.tap_l = 143;
    b10.tap_r = 164;
    b11.d1    = 449;
    b11.d2    = 293;
    b11.dd    = 743;
    b11.fac_l = -0.30971853768214;
    b11.fac_r = 0.69028146231786;
    b11.k1    = -coeff1;
    b11.k2    = -coeff2;
    b11.tap_l = 38;
    b11.tap_r = 62;
    b12.d1    = 421;
    b12.d2    = 233;
    b12.dd    = 709;
    b12.fac_l = -0.24329721095117;
    b12.fac_r = 0.75670278904883;
    b12.k1    = -coeff1;
    b12.k2    = coeff2;
    b12.tap_l = 167;
    b12.tap_r = 151;

    // on the prototype fetch is done after pushing, on here it's done before
    b1.tap_l -= 1;
    b2.tap_l -= 1;
    b3.tap_l -= 1;
    b4.tap_l -= 1;
    b5.tap_l -= 1;
    b6.tap_l -= 1;
    b7.tap_l -= 1;
    b8.tap_l -= 1;
    b9.tap_l -= 1;
    b10.tap_l -= 1;
    b11.tap_l -= 1;
    b12.tap_l -= 1;
    b1.tap_r -= 1;
    b2.tap_r -= 1;
    b3.tap_r -= 1;
    b4.tap_r -= 1;
    b5.tap_r -= 1;
    b6.tap_r -= 1;
    b7.tap_r -= 1;
    b8.tap_r -= 1;
    b9.tap_r -= 1;
    b10.tap_r -= 1;
    b11.tap_r -= 1;
    b12.tap_r -= 1;

    return make_array<stage_data> (
      make_lp (0.15), // 0
      make_hp (0.99), // 1
      make_lp (0.15), // 2
      make_hp (0.99), // 3

      make_ap (b1.d1, b1.k1, 43), // 4
      make_ap (b1.d2, b1.k2), // 5
      make_parallel_delay (
        2, b1.dd, b1.tap_l, b1.fac_l, b1.tap_r, b1.fac_r), // 6
      make_crossover2(), // 7

      make_ap (b2.d1, b2.k1), // 8
      make_ap (b2.d2, b2.k2), // 9
      make_parallel_delay (
        2, b2.dd, b2.tap_l, b2.fac_l, b2.tap_r, b2.fac_r), // 10

      make_ap (b3.d1, b3.k1), // 11
      make_ap (b3.d2, b3.k2), // 12
      make_parallel_delay (
        2, b3.dd, b3.tap_l, b3.fac_l, b3.tap_r, b3.fac_r), // 13

      make_ap (b4.d1, b4.k1), // 14
      make_ap (b4.d2, b4.k2), // 15
      make_parallel_delay (
        2, b4.dd, b4.tap_l, b4.fac_l, b4.tap_r, b4.fac_r), // 16

      make_ap (b5.d1, b5.k1, 37), // 17
      make_ap (b5.d2, b5.k2), // 18
      make_parallel_delay (
        2, b5.dd, b5.tap_l, b5.fac_l, b5.tap_r, b5.fac_r), // 19
      make_crossover2(), // 20

      make_ap (b6.d1, b6.k1), // 21
      make_ap (b6.d2, b6.k2), // 22
      make_parallel_delay (
        2, b6.dd, b6.tap_l, b6.fac_l, b6.tap_r, b6.fac_r), // 23

      make_ap (b7.d1, b7.k1), // 24
      make_ap (b7.d2, b7.k2), // 25
      make_parallel_delay (
        2, b7.dd, b7.tap_l, b7.fac_l, b7.tap_r, b7.fac_r), // 26

      make_ap (b8.d1, b8.k1), // 27
      make_ap (b8.d2, b8.k2), // 28
      make_parallel_delay (
        2, b8.dd, b8.tap_l, b8.fac_l, b8.tap_r, b8.fac_r), // 29

      make_ap (b9.d1, b9.k1, 27), // 30
      make_ap (b9.d2, b9.k2), // 31
      make_parallel_delay (
        2, b9.dd, b9.tap_l, b9.fac_l, b9.tap_r, b9.fac_r), // 32
      make_crossover2(), // 33

      make_ap (b10.d1, b10.k1), // 34
      make_ap (b10.d2, b10.k2), // 35
      make_parallel_delay (
        2, b10.dd, b10.tap_l, b10.fac_l, b10.tap_r, b10.fac_r), // 36

      make_ap (b11.d1, b11.k1), // 37
      make_ap (b11.d2, b11.k2), // 38
      make_parallel_delay (
        2, b11.dd, b11.tap_l, b11.fac_l, b11.tap_r, b11.fac_r), // 39

      make_ap (b12.d1, b12.k1), // 40
      make_ap (b12.d2, b12.k2), // 41
      make_parallel_delay (
        2,
        b12.dd,
        b12.tap_l,
        b12.fac_l,
        b12.tap_r,
        b12.fac_r) // 42
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
    uint                         srate)
  {
    block_arr<sample> l_in, r_in, m_in, loopm, l1, r1, l2, r2, lfo1, lfo2, lfo3;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo     = tick_lfo<sample> (_lfo);
      lfo1[i]      = sample {lfo[0]} * par.mod[i];
      lfo2[i]      = sample {lfo[1]} * par.mod[i];
      lfo3[i]      = sample {lfo[2]} * par.mod[i];
      l_in[i]      = io[i][0] * 0.125_r * 0.25_r;
      r_in[i]      = io[i][1] * 0.125_r * 0.25_r;
      par.decay[i] = fastgrowth (par.decay[i]);
      par.decay[i] = 0.4_r + par.decay[i] * 0.59_r;
    }
    run_cascade<0, 1> (_eng, xspan {l_in.data(), io.size()});
    run_cascade<2, 3> (_eng, xspan {r_in.data(), io.size()});

    auto   lf   = upar.lf_amt;
    auto   hf   = upar.hf_amt;
    sample flo  = load_float<sample> (0.8f + lf * lf * 0.1f);
    sample glo  = load_float<sample> (0.75f + lf * 0.245f);
    sample fhi1 = load_float<sample> (0.9f - hf * hf * 0.4f);
    sample fhi2 = load_float<sample> (0.93f - hf * hf * 0.4f);
    sample fhi3 = load_float<sample> (0.95f - hf * hf * 0.45f);
    sample ghi  = load_float<sample> (
      0.85f + hf * 0.1f + as_float (par.decay[0]) * 0.05f);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      m_in[i] = (l_in[i] + r_in[i]) * 0.5_r;
    }

    xspan loop {loopm.data(), io.size()};

    // 709 = b12.dd
    _eng.fetch (sl<42> {}, loop);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<4, 5> {}, loop, lfo1);
    _eng.run (sl<6> {}, loop, loop.to_const(), overwrite, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<7> {}, loop, flo, glo, fhi1, 1_r, ghi);

    _eng.run (sl<8, 9> {}, loop);
    _eng.run (sl<10> {}, loop, loop.to_const(), overwrite, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, l_in);

    _eng.run (sl<11, 12> {}, loop);
    _eng.run (sl<13> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);

    _eng.run (sl<14, 15> {}, loop);
    _eng.run (sl<16> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<17, 18> {}, loop, lfo2);
    _eng.run (sl<19> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<20> {}, loop, flo, glo, fhi2, 1_r, ghi);

    _eng.run (sl<21, 22> {}, loop);
    _eng.run (sl<23> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, r_in);

    _eng.run (sl<24, 25> {}, loop);
    _eng.run (sl<26> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);

    _eng.run (sl<27, 28> {}, loop);
    _eng.run (sl<29> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add (loop, m_in);

    _eng.run (sl<30, 31> {}, loop, lfo3);
    _eng.run (sl<32> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);
    _eng.run (sl<33> {}, loop, flo, glo, fhi3, 1_r, ghi);

    _eng.run (sl<34, 35> {}, loop);
    _eng.run (sl<36> {}, loop, loop.to_const(), add_to, l2, r1);
    span_mul (loop, par.decay);

    span_add_with_factor (loop, m_in, -1_r);

    _eng.run (sl<37, 38> {}, loop);
    _eng.run (sl<39> {}, loop, loop.to_const(), add_to, l1, r2);
    span_mul (loop, par.decay);

    _eng.run (sl<40, 41> {}, loop);
    _eng.run (sl<42> {}, loop.to_const(), add_to, l2, r1);

    auto v1 = make_array (&l1[0], &r1[0]);
    auto v2 = make_array<sample const*> (&l2[0], &r2[0]);
    ep_crossfade<sample> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return 0.1_r + v * 0.9_r;
    });
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l1[i];
      spls[1] = r1[i];
    });
  }
  //----------------------------------------------------------------------------
  void post_process_block (xspan<std::array<float, 2>> io)
  {
    span_visit (io, [&] (auto& v, uint) {
      v = vec_to_array (_eq.tick_cascade (vec_from_array (v)) * 8.f);
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
class chorus_a : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<chorus_a<Dt>, Dt, max_block_size>;

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
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.0f, 0.f, 0.0f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (440.f),
      vec_set<f32_x2> (0.28f),
      vec_set<f32_x2> (2.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    mod     = mod * mod * mod;
    auto f1 = 0.14f + mod * 6.7f;
    auto f2 = 0.13f + mod * 6.9f;
    auto f3 = 0.17f + mod * 6.4f;
    auto f4 = 0.15f + mod * 6.7f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f4}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp(), // 0
      make_hp(), // 1
      make_lp(), // 2
      make_hp(), // 3
      make_ap (73, 0.5), // 4
      make_ap (77, -0.5), // 5
      make_ap (631, 0.2, 507), // 6
      make_ap (641, -0.2, 513), // 7
      make_ap (967, -0.1, 409), // 8
      make_ap (937, 0.1, 379) // 9
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
    using arr = block_arr<sample>;
    arr sig1, sig2, aux1, aux2, aux3, aux4;

    auto l = xspan {sig1.data(), io.size()};
    auto r = xspan {sig2.data(), io.size()};

    for (uint i = 0; i < io.size(); ++i) {
      auto mod = 0.2_r + par.decay[i] * 0.8_r;
      mod *= (1_r - par.mod[i] * 0.7_r);
      mod *= mod;
      auto lfo = tick_lfo<sample> (_lfo);
      aux1[i]  = sample {lfo[0]} * mod;
      aux2[i]  = sample {lfo[1]} * mod;
      aux3[i]  = sample {lfo[2]} * -mod;
      aux4[i]  = sample {lfo[3]} * -mod;
      l[i]     = io[i][0];
      r[i]     = io[i][1];
    }
    auto hi_g
      = load_float<sample> (0.1f + (1.f - upar.hf_amt * upar.hf_amt) * 0.65f);
    auto lo_g = load_float<sample> (0.8f + 0.18f * fastgrowth (upar.lf_amt));

    _eng.run (sl<0> {}, l, hi_g);
    _eng.run (sl<1> {}, l, lo_g);
    _eng.run (sl<2> {}, r, hi_g);
    _eng.run (sl<3> {}, r, lo_g);

    _eng.run (sl<4> {}, l, blank, [&] (uint i) {
      return par.character[i] * 0.5_r;
    });
    _eng.run (sl<5> {}, r, blank, [&] (uint i) {
      return par.character[i] * -0.5_r;
    });

    auto m = xspan {sig1.data(), io.size()};
    auto g = xspan {par.character.data(), io.size()};
    span_visit (m, [&] (auto& v, uint i) {
      v    = (l[i] + r[i]) * 0.5_r;
      g[i] = g[i] * 0.2_r;
    });
    _eng.run (sl<6> {}, aux1, m.to_const(), aux1, par.character);

    span_visit (g, [&] (auto& v, uint i) { v = -v; });
    _eng.run (sl<7> {}, aux2, m.to_const(), aux2, par.character);

    span_visit (g, [&] (auto& v, uint i) { v = v * 0.5_r; });
    _eng.run (sl<8> {}, aux3, m.to_const(), aux3, par.character);

    span_visit (g, [&] (auto& v, uint i) { v = -v; });
    _eng.run (sl<9> {}, aux4, m.to_const(), aux4, par.character);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = aux1[i] + aux3[i];
      spls[1] = aux4[i] - aux2[i];
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
class chorus_b : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<chorus_b<Dt>, Dt, max_block_size>;

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
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.0f, 0.f, 0.0f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (440.f),
      vec_set<f32_x2> (0.28f),
      vec_set<f32_x2> (2.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto m1 = mod * mod;
    auto m2 = m1 + mod;
    auto f1 = 0.147f + m1 * 0.723f;
    auto f2 = 0.131f + m1 * 0.913f;
    auto f3 = 0.473f + m2 * 3.763f;
    auto f4 = 0.457f + m2 * 3.427f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f4}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp(), // 0
      make_hp(), // 1
      make_ap (2331, 0.2, 407), // 2
      make_ap (2341, -0.2, 413), // 3
      make_ap (1967, -0.15, 509), // 4
      make_ap (1937, 0.15, 579), // 5
      make_ap (173, 0.5), // 6
      make_ap (177, 0.5) // 7
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
    using arr = block_arr<sample>;
    arr sig1, sig2, aux1, aux2, aux3, aux4;

    auto m = xspan {sig1.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      auto mod = 0.2_r + par.decay[i] * 0.8_r;
      mod *= (1_r - par.mod[i] * 0.8_r);
      mod *= mod;
      auto lfo = tick_lfo<sample> (_lfo);
      aux2[i]  = sample {lfo[0] - lfo[2]} * mod * 0.5_r;
      aux4[i]  = sample {lfo[1] + lfo[3]} * mod * 0.5_r;
      m[i]     = (io[i][0] + io[i][1]) * 0.5_r;
    }

    auto hi_g
      = load_float<sample> (0.1f + (1.f - upar.hf_amt * upar.hf_amt) * 0.65f);
    auto lo_g = load_float<sample> (0.8f + 0.18f * fastgrowth (upar.lf_amt));
    _eng.run (sl<0> {}, m, hi_g);
    _eng.run (sl<1> {}, m, lo_g);

    _eng.run (sl<2> {}, aux1, m.to_const(), aux2, [&] (uint i) {
      return par.character[i] * 0.2_r;
    });
    span_visit (xspan {aux2.data(), io.size()}, [&] (auto& v, uint i) {
      v = -v;
    });
    _eng.run (sl<3> {}, aux2, m.to_const(), aux2, [&] (uint i) {
      return par.character[i] * -0.2_r;
    });
    _eng.run (sl<4> {}, aux3, m.to_const(), aux4, [&] (uint i) {
      return par.character[i] * -0.15_r;
    });
    span_visit (xspan {aux4.data(), io.size()}, [&] (auto& v, uint i) {
      v = -v;
    });
    _eng.run (sl<5> {}, aux4, m.to_const(), aux4, [&] (uint i) {
      return par.character[i] * 0.15_r;
    });

    auto l = xspan {sig1.data(), io.size()};
    auto r = xspan {sig2.data(), io.size()};
    span_visit (l, [&] (auto&, uint i) {
      l[i] = aux1[i] + aux3[i];
      r[i] = aux4[i] - aux2[i];
    });

    _eng.run (sl<6> {}, l, blank, [&] (uint i) {
      return par.character[i] * 0.5_r;
    });
    _eng.run (sl<7> {}, r, blank, [&] (uint i) {
      return par.character[i] * -0.5_r;
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
}}} // namespace artv::detail::tpaco
