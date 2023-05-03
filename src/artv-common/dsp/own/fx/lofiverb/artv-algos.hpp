#pragma once

#include <type_traits>

#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"
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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
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
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
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
      late_in[i] = (sample) ((io[i][0] + io[i][1]) * 0.5_r);
      auto mod   = (sample) (0.25_r + (1_r - par.mod[i]) * 0.75_r);
      // ER + late lfo
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = sample {lfo[0]};
      lfo2[i] = (sample) (sample {lfo[1]} * (0.5_r + par.character[i] * 0.5_r));
      lfo3[i] = (sample) (sample {lfo[2]} * mod);
      lfo4[i] = (sample) (sample {lfo[3]} * mod);

      // decay fixup
      auto decay   = (sample) (_eng.one - par.decay[i]);
      decay        = (sample) (_eng.one - decay * decay);
      par.decay[i] = (sample) (0.6_r + decay * 0.38_r);
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

    _eng.fetch_block (sl<9> {}, er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = (sample) (late_in[i] * 0.5_r + er2[i]);
    });
    er2.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<5> {}, er1, xspan {lfo2});
    _eng.run (sl<6> {}, er1, flo, glo, fhi, _eng.one, ghi);
    xspan_memcpy (er1b, er1);
    _eng.run (sl<7> {}, er1b, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    _eng.run (sl<8> {}, er1b, xspan {lfo1});
    _eng.push (sl<9> {}, er1b.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto e_l = (sample) (-er1[i] * 0.66_r - er2[i] * 0.34_r);
      auto e_r = (sample) (-er1[i] * 0.66_r + er2[i] * 0.34_r);
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
    _eng.fetch_block (sl<18> {}, l, 1); // feedback, fetching block + 1 samples
    _eng.fetch_block (sl<28> {}, r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (sample) (loopsig * (_eng.one - er_amt) + er_sig * er_amt);
      g[i]
        = (sample) (0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r));
    }
    r.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<10> {}, late, xspan {lfo3}, g);
    _eng.run (sl<11> {}, late);
    _eng.run (sl<12> {}, late, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<13> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    _eng.run (sl<14> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<15> {}, late);
    _eng.run (sl<16> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    _eng.run (sl<17> {}, late);
    _eng.push (sl<18> {}, late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = (sample) (loopsig * (_eng.one - er_amt) + er_sig * er_amt);
    }
    l.cut_head (1); // drop feedback sample from previous block
    _eng.run (sl<19> {}, late, xspan {lfo4}, g);
    _eng.run (sl<20> {}, late);
    _eng.run (sl<21> {}, late, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<22> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    _eng.run (sl<23> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<24> {}, late);
    _eng.run (sl<25> {}, late, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    _eng.run (sl<26> {}, late);
    _eng.run (sl<27> {}, late);
    _eng.push (sl<28> {}, late.to_const()); // feedback point

    // Mixdown
    auto v1 = make_array (&l[0], &er1[0]);
    auto v2 = make_array<sample const*> (&r[0], &er2[0]);
    ep_crossfade<sample> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return fastgrowth ((sample) (v * 0.1_r), _eng.one);
    });

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
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
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
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
    _eng.fetch_block (sl<50> {}, a, 1); // feedback, fetching block + 1 samples
    _eng.fetch_block (sl<64> {}, d, 1); // feedback, fetching block + 1 samples

    xspan_memcpy (b, a.advanced (1)); // L out on b
    a.cut_tail (1); // fb samples on a.

    xspan_memcpy (c, d.advanced (1)); // R out on c
    d.cut_tail (1); // fb samples on d.

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // invert character first
      par.character[i] = (sample) (_eng.one - par.character[i]);
      tmp1[i]          = (sample) (0.4_r + 0.3_r * par.character[i]);
      tmp2[i]          = (sample) (-0.6_r * par.character[i]);
      tmp3[i]          = (sample) (0.6_r * par.character[i]);
      // decay fixup
      auto d       = _eng.one - par.decay[i];
      par.decay[i] = (sample) (0.985_r - d * d * 0.21_r);
      // lfo
      auto lfo = tick_lfo<sample> (lfo_obj);
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
      return fastgrowth ((sample) (v * 0.5_r), _eng.one);
    });

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      b[i] = (sample) (b[i] + cho1[i]);
      c[i] = (sample) (c[i] + cho2[i]);
    }

    _eng.run (sl<10> {}, b, blank, tmp3);
    _eng.run (sl<11> {}, c, blank, tmp3); // tmp3 free

    xspan i1 {tmp1.data(), io.size()};
    xspan i2 {tmp2.data(), io.size()};
    // output write + preparations
    span_visit (io, [&] (auto& spls, uint i) {
      i1[i]   = (sample) ((spls[0] + spls[1]) * 0.25_r); // m
      i2[i]   = (sample) ((spls[0] - spls[1]) * 0.25_r); // s
      spls[0] = b[i];
      spls[1] = c[i];
      b[i]    = (sample) (0.8_r + 0.1_r * par.character[i]); // character 1 on b
      c[i]    = (sample) (0.8_r * par.character[i]); // character 2 on c
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
    _eng.fetch_block (sl<55> {}, b, 1); // feedback
    _eng.fetch_block (sl<60> {}, c, 1); // feedback

    // quantized decay
    _eng.run (sl<22> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<23> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<24> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<25> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

    // sum inputs to feedback
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // a single hadamard iteration can duplicate the range on one
      // of the channels, so halving once more (TODO needed?).
      b[i] = (sample) (b[i] + tmp1[i] * 0.5_r);
      d[i] = (sample) (d[i] + tmp2[i] * 0.5_r);
    }
    hadamard4 (make_array (a.data(), b.data(), c.data(), d.data()), a.size());

    // channel a block 1
    _eng.run (sl<26, 27, 28> {}, a);
    _eng.run (sl<29> {}, a);

    // channel b block 1
    _eng.run (sl<30> {}, b);

    // channel c block 1
    _eng.run (sl<31, 32, 33> {}, c);
    _eng.run (sl<34> {}, c);

    // channel d block 1
    _eng.run (sl<35, 36, 37> {}, d);
    _eng.run (sl<38> {}, d);

    // quantized decay
    _eng.run (sl<39> {}, a, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<40> {}, b, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<41> {}, c, [&] (auto v, uint i) { return v * par.decay[i]; });
    _eng.run (sl<42> {}, d, [&] (auto v, uint i) { return v * par.decay[i]; });

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
    sample ghi = load_float<sample> (0.6f + upar.hf_amt * 0.18f);

    _eng.run (sl<43> {}, a, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<44> {}, b, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<45> {}, c, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<46> {}, d, flo, glo, fhi, _eng.one, ghi);

    // channel a block 2
    _eng.run (sl<47> {}, a);
    _eng.run (sl<48> {}, a);
    _eng.run (sl<49> {}, a);
    _eng.push (sl<50> {}, a.to_const());

    // channel b block 2
    _eng.run (sl<51, 52, 53> {}, b);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      sample m = fastgrowth (par.mod[i], _eng.one);
      tmp1[i]  = (sample) (0.4_r + 0.4_r * lfo2[i] * m);
      tmp2[i]  = (sample) (0.6_r + 0.4_r * -lfo2[i] * m); // for block c
    }
    _eng.run (sl<54> {}, b, tmp1);
    _eng.push (sl<55> {}, b.to_const());

    // channel c block 2
    _eng.run (sl<56, 57, 58> {}, c);
    _eng.run (sl<59> {}, c, tmp2);
    _eng.push (sl<60> {}, c.to_const());

    // channel d block 2
    _eng.run (sl<61, 62, 63> {}, d);
    _eng.push (sl<64> {}, d.to_const());
  }

private:
  engine _eng;
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
template <delay::data_type Dt>
class room : public algorithm {
private:
  using engine = detail::lofiverb::algo_engine<room<Dt>, Dt, max_block_size>;
  static constexpr std::array<u8, 32> ffwd_table {get_room_ffwd_table()};

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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
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
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
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
      l[i] = (sample) (spl[0] * 0.25_r);
      r[i] = (sample) (spl[1] * 0.25_r);
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
      spl      = (sample) ((l[i] + r[i]) * 0.25_r);
      auto inv = _eng.one - par.decay[i];
      l[i]     = (sample) (l[i] + t1[i] * inv);
      r[i]     = (sample) (r[i] - t2[i] * inv);
      ltdec[i] = (sample) (par.decay[i] * par.decay[i]);
      t1[i]    = (sample) (0.39_r * ltdec[i]);
      t2[i]    = (sample) (-0.28_r * ltdec[i]);
      t3[i]    = (sample) (0.19_r * ltdec[i]);
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
      _eng.one,
      load_float<sample> (0.35f + upar.hf_amt * 0.14f),
      par.character,
      t3);

    span_visit (xspan {ltdec.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (sample) (0.3_r + v * 0.3_r);
      t2[i] = (sample) (-0.25_r - v * 0.3_r);
      t3[i] = (sample) (0.25_r + v * 0.3_r);
    });
    _eng.run (sl<24> {}, m, blank, t1);
    _eng.run (sl<25> {}, m, blank, t2);
    _eng.run (sl<26> {}, m, blank, t3);

    span_visit (xspan {par.mod.data(), io.size()}, [&] (auto& v, uint i) {
      t1[i] = (sample) (sample {tick_lfo<sample> (lfo_obj)[0]} * v);
      t2[i] = (sample) (0.365_r + v * 0.2_r);
    });
    _eng.run (sl<27> {}, xspan {t3.data(), io.size()}, m.to_const(), t1, t2);
    span_visit (io, [&] (auto& spl, uint i) {
      spl[0] = (sample) (l[i] + par.decay[i] * t3[i]);
      spl[1] = (sample) (r[i] - par.decay[i] * m[i]);
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f1 = 1.73f - mod * 0.63f;
    auto f2 = 1.53f - mod * 0.43f;
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
      make_lp (0.26f), // 0
      make_hp (0.975f), // 1
      make_lp (0.26f), // 2
      make_hp (0.975f), // 3

      make_block_delay (1447), // feedback point -> h1 4
      make_quantizer(), // 5

      make_ap (376, 0, 27), // 6
      make_delay (1007), // 7 -> a1

      make_crossover2(), // 8
      make_ap (363), // 9
      make_delay (1107), // 10 -> b1
      make_quantizer(), // 11

      make_ap (414, 0, 23), // 12
      make_delay (1207), // 13 -> c1

      make_crossover2(), // 14
      make_ap (477), // 15
      make_delay (1307), // 16 -> d1
      make_quantizer(), // 17

      make_ap (420, 0, 21), // 18
      make_delay (1407), // 19 -> e1

      make_crossover2(), // 20
      make_ap (252), // 21
      make_delay (1507), // 22 -> f1
      make_quantizer(), // 23

      make_ap (413, 0, 22), // 24
      make_delay (1607), // 25 -> g1

      make_crossover2(), // 26
      make_ap (813) // 27
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<sample>;
    using arr_fb = fb_block_arr<sample>;

    arr loop_mem;
    arr l, r, r_cp, m, s, k1, k2, tmp_mem, loop2_mem;
    arr lfo1, lfo2;

    xspan loop {loop_mem.data(), io.size()};
    xspan loop2 {loop2_mem.data(), io.size()};
    xspan tmp {tmp_mem.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      l[i]    = (sample) (spl[0] * 0.25_r);
      r[i]    = (sample) (spl[1] * 0.25_r);
      r_cp[i] = r[i];
    });
    _eng.run (sl<0> {}, xspan {l.data(), io.size()});
    _eng.run (sl<1> {}, xspan {l.data(), io.size()});
    _eng.run (sl<2> {}, xspan {r.data(), io.size()});
    _eng.run (sl<3> {}, xspan {r.data(), io.size()});
    _eng.fetch_block (sl<4> {}, loop, 1); // feedback signal

    constexpr auto sqrt8_recip = 0.3535534_r; // 1/sqrt(8), equal power mixing

    _eng.run (sl<5> {}, loop, [&] (auto v, uint i) {
      par.decay[i] = (sample) (0.05_r + par.decay[i] * 0.92_r);
      return v * par.decay[i];
    });
    span_visit (loop, [&] (auto& v, uint i) {
      m[i]     = (sample) ((l[i] + r[i]) * 0.5_r);
      s[i]     = (sample) ((l[i] - r[i]) * 0.5_r);
      k1[i]    = (sample) (0.4_r + par.character[i] * 0.2_r);
      k2[i]    = (sample) (0.4_r + par.character[i] * 0.15_r);
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      v        = (sample) (v + l[i] * 0.75_r + m[i] * 0.25_r);
    });

    auto flo = load_float<sample> (0.90 + upar.lf_amt * upar.lf_amt * 0.05);
    auto glo = load_float<sample> (0.75 + upar.lf_amt * 0.23);
    auto fhi = load_float<sample> (0.9 - upar.hf_amt * upar.hf_amt * 0.4);
    auto ghi = load_float<sample> (0.7 + upar.hf_amt * 0.23);
    // A
    _eng.run (sl<6> {}, loop, lfo1, k1);
    xspan_memdump (l.data(), loop);
    _eng.fetch (sl<7> {}, r, 400 - 32);
    span_mul (xspan {l.data(), io.size()}, sqrt8_recip);
    span_mul (xspan {r.data(), io.size()}, sqrt8_recip * 0.4_r);
    _eng.run (sl<7> {}, loop);
    // B
    _eng.run (sl<8> {}, loop, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<9> {}, loop, blank, k1);
    _eng.fetch (sl<10> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] + v * 0.4_r * sqrt8_recip);
      r[i] = (sample) (r[i] + loop[i] * sqrt8_recip);
    });
    _eng.run (sl<10> {}, loop);
    _eng.run (sl<11> {}, loop, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    // C
    span_visit (loop, [&] (auto& v, uint i) { v = (sample) (s[i] + v); });
    _eng.run (sl<12> {}, loop, lfo2, [&] (uint i) { return -k2[i]; });
    _eng.fetch (sl<13> {}, tmp, 1001 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] - loop[i] * sqrt8_recip);
      r[i] = (sample) (r[i] - v * 0.42_r * sqrt8_recip);
    });
    _eng.run (sl<13> {}, loop);
    // D
    _eng.run (sl<14> {}, loop, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<15> {}, loop, blank, k2);
    _eng.fetch (sl<16> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] + v * 0.42_r * sqrt8_recip);
      r[i] = (sample) (r[i] - loop[i] * 0.49_r * sqrt8_recip);
    });
    _eng.run (sl<16> {}, loop);
    _eng.run (sl<17> {}, loop, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    // E
    span_visit (loop, [&] (auto& v, uint i) {
      v = (sample) (v + r_cp[i] * 0.75_r + m[i] * 0.25_r);
    });
    _eng.run (
      sl<18> {}, loop, [&] (uint i) { return -lfo1[i]; }, k1);
    _eng.fetch (sl<19> {}, tmp, 801 - 37 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] + loop[i] * 0.49_r * sqrt8_recip);
      r[i] = (sample) (r[i] - v * sqrt8_recip);
    });
    _eng.run (sl<19> {}, loop);
    // F
    _eng.run (sl<20> {}, loop, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<21> {}, loop, blank, k1);
    _eng.fetch (sl<22> {}, tmp, 777 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] + v * sqrt8_recip);
      r[i] = (sample) (r[i] + loop[i] * sqrt8_recip);
    });
    _eng.run (sl<22> {}, loop);
    _eng.run (sl<23> {}, loop, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    // G
    span_visit (loop, [&] (auto& v, uint i) { v = (sample) (v + s[i]); });
    _eng.run (
      sl<24> {}, loop, [&] (uint i) { return -lfo2[i]; }, k2);
    _eng.fetch (sl<25> {}, tmp, 1001 - 27 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] - loop[i] * sqrt8_recip);
      r[i] = (sample) (r[i] - v * sqrt8_recip);
    });
    _eng.run (sl<25> {}, loop);
    // H
    _eng.run (sl<26> {}, loop, flo, glo, fhi, _eng.one, ghi);
    _eng.run (sl<27> {}, loop, blank, k2);
    _eng.fetch (sl<4> {}, tmp, 1001 - 32);
    span_visit (tmp, [&] (auto v, uint i) {
      l[i] = (sample) (l[i] + loop[i] * 0.2_r * sqrt8_recip);
      r[i] = (sample) (r[i] + v * 0.2_r * sqrt8_recip);
    });
    _eng.push (sl<4> {}, loop.to_const());

    span_visit (io, [&] (auto& spl, uint i) {
      spl[0] = l[i];
      spl[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
}}} // namespace artv::detail::lofiverb
