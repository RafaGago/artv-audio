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
