#pragma once

#include "artv-common/dsp/own/fx/lofiverb/algorithm-engine.hpp"
#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/primes_table.hpp"
#include "artv-common/misc/xspan.hpp"
#include <type_traits>

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
    stage b1 {}, b2 {}, b3 {}, b4 {}, b5 {}, b6 {}, b7 {}, b8 {};
    b1.d1    = 571;
    b1.d2    = 359;
    b1.dd    = 857;
    b1.fac_l = -0.40780916679832;
    b1.fac_r = 0.59219083320168;
    b1.k1    = -0.33;
    b1.k2    = 0.14;
    b1.tap_l = 329;
    b1.tap_r = 343;
    b2.d1    = 467;
    b2.d2    = 257;
    b2.dd    = 707;
    b2.fac_l = -0.49429981212465;
    b2.fac_r = 0.50570018787535;
    b2.k1    = -0.32;
    b2.k2    = 0.14;
    b2.tap_l = 316;
    b2.tap_r = 286;
    b3.d1    = 577;
    b3.d2    = 383;
    b3.dd    = 941;
    b3.fac_l = 0.25340987584866;
    b3.fac_r = -0.74659012415134;
    b3.k1    = -0.33;
    b3.k2    = 0.14;
    b3.tap_l = 492;
    b3.tap_r = 456;
    b4.d1    = 449;
    b4.d2    = 283;
    b4.dd    = 757;
    b4.fac_l = -0.77724914480402;
    b4.fac_r = 0.22275085519598;
    b4.k1    = 0.32;
    b4.k2    = 0.14;
    b4.tap_l = 287;
    b4.tap_r = 295;
    b5.d1    = 619;
    b5.d2    = 353;
    b5.dd    = 991 + 100;
    b5.fac_l = -0.4534214228516;
    b5.fac_r = 0.5465785771484;
    b5.k1    = -0.33;
    b5.k2    = 0.13;
    b5.tap_l = 353;
    b5.tap_r = 388;
    b6.d1    = 443;
    b6.d2    = 241;
    b6.dd    = 691 + 100;
    b6.fac_l = 0.38615093528902;
    b6.fac_r = -0.61384906471098;
    b6.k1    = 0.32;
    b6.k2    = -0.14;
    b6.tap_l = 285;
    b6.tap_r = 264;
    b7.d1    = 599;
    b7.d2    = 373;
    b7.dd    = 919;
    b7.fac_l = -0.93732491320402;
    b7.fac_r = 0.06267508679598;
    b7.k1    = -0.33;
    b7.k2    = 0.14;
    b7.tap_l = 424;
    b7.tap_r = 459;
    b8.d1    = 467;
    b8.d2    = 269;
    b8.dd    = 709 + 84;
    b8.fac_l = -0.44557774496395;
    b8.fac_r = 0.554422;
    b8.k1    = -0.32;
    b8.k2    = 0.14;
    b8.tap_l = 302;
    b8.tap_r = 265;

    // on the prototype fetch is done after pushing, on here it's done before
    b1.tap_l -= 1;
    b2.tap_l -= 1;
    b3.tap_l -= 1;
    b4.tap_l -= 1;
    b5.tap_l -= 1;
    b6.tap_l -= 1;
    b7.tap_l -= 1;
    b8.tap_l -= 1;
    b1.tap_r -= 1;
    b2.tap_r -= 1;
    b3.tap_r -= 1;
    b4.tap_r -= 1;
    b5.tap_r -= 1;
    b6.tap_r -= 1;
    b7.tap_r -= 1;
    b8.tap_r -= 1;

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
        2, b8.dd, b8.tap_l, b8.fac_l, b8.tap_r, b8.fac_r) // 31
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
    sample ghi  = load_float<sample> (
      0.85f + hf * 0.1f + as_float (par.decay[0]) * 0.05f);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      m_in[i] = (l_in[i] + r_in[i]) * 0.5_r;
    }

    xspan loop {loopm.data(), io.size()};

    _eng.fetch (sl<31> {}, loop);
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
}}} // namespace artv::detail::lofiverb
