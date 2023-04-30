#pragma once

#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"

namespace artv { namespace detail { namespace lofiverb {

/* reminder_ sample rates that result in a low number of interpolator/decimator
  tables.
   srate t44k t48k
  7200 49    20
  8400 21    40
  9000 49    16
  10500 21   32
  10800 49   40
  13500 49   32
  14400 49   10
  16800 21   20
  18000 49   08
  20160 35   50
  21000 21   16
  21600 49   20
  22500 49   32
  25200 07   40
  27000 49   16
  28800 49   05
  31500 07   32
  32400 49   40
  33600 21   10
  36000 49   04
  39200 09   60
  39600 49   40
  40320 35   25
  40500 49   32
  42000 21   08
  43200 49   10
  45000 50   16
  46800 52   40
  47040 16   50
  49000 10   49
  49500 55   33
  50400 08   21
  52500 25   35
  54000 60   09
  57600 64   06
  58500 65   39
  58800 04   49
  60480 48   63
  61200 68   51
  63000 10   21
  64800 72   27
  67200 32   07
  67500 75   45
  68400 76   57
*/
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class debug : public algorithm {
private:
  using engine = detail::lofiverb::algo_engine<debug<Dt>, Dt, max_block_size>;

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
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------

}}} // namespace artv::detail::lofiverb
