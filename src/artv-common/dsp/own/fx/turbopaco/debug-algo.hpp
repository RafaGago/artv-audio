#pragma once

#include "artv-common/dsp/own/fx/turbopaco/algorithm-engine.hpp"
#include "artv-common/dsp/own/fx/turbopaco/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/primes_table.hpp"
#include "artv-common/misc/xspan.hpp"
#include <type_traits>

namespace artv { namespace detail { namespace tpaco {

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
// TODO: don't use 23400 as base sample rate anymore. 25600 and 21600 require
// less sync interpolation tables.
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class debug : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<debug<Dt>, Dt, max_block_size>;

  struct stage {
    uint  d1, d2, d3, dd1, dd2;
    float k1, k2, k3;
  };

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (
      8400, 10500, 13500, 16800, 21000, 23400, 33600, 42000, 50400);
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
      f32_x2 {440.f, 460.f},
      vec_set<f32_x2> (0.48f),
      vec_set<f32_x2> (1.5f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 2.97f - mod * 1.43f;
    auto f2 = 0.21f + mod * 0.53f;
    auto f3 = 3.34f - mod * 1.17f;
    auto f4 = 0.84f + mod * 1.01f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f4}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    // just to avoid copy-paste issues from the prototype
    stage b1 {}, b2 {}, b3 {}, b4 {};
    b1.d1  = 634;
    b1.d2  = 787;
    b1.d3  = 474;
    b1.dd1 = 797;
    b1.dd2 = 1051;
    b1.k1  = -0.5;
    b1.k2  = -0.5;
    b1.k3  = -0.3;
    b2.d1  = 945;
    b2.d2  = 1191;
    b2.d3  = 689;
    b2.dd1 = 1237;
    b2.dd2 = 1571;
    b2.k1  = 0.475;
    b2.k2  = 0.475;
    b2.k3  = -0.285;
    b3.d1  = 1568;
    b3.d2  = 1998;
    b3.d3  = 1166;
    b3.dd1 = 2011;
    b3.dd2 = 2617;
    b3.k1  = -0.45;
    b3.k2  = 0.45;
    b3.k3  = 0.27;

    return make_array<stage_data> (
      make_lp (0.15), // 0
      make_hp (0.985), // 1
      make_lp (0.15), // 2
      make_hp (0.985), // 3

      make_ap (173, 0.5), // 4
      make_ap (173, 0.5), // 5

      make_crossover2(), // 6
      make_ap (b1.d1, b1.k1, 27), // 7
      make_delay (b1.dd1), // 8
      make_ap (b1.d2, b1.k2), // 9
      make_ap (b1.d3, b1.k3), // 10
      make_delay (b1.dd2), // 11

      make_crossover2(), // 12
      make_ap (b2.d1, b2.k1, 18), // 13
      make_delay (b2.dd1), // 14
      make_ap (b2.d2, b2.k2), // 15
      make_ap (b2.d3, b2.k3), // 16
      make_delay (b2.dd2), // 17

      make_crossover2(), // 18
      make_ap (b3.d1, b3.k1, 21), // 19
      make_delay (b3.dd1), // 20
      make_ap (b3.d2, b3.k2), // 21
      make_ap (b3.d3, b3.k3), // 22
      make_delay (b3.dd2), // 23

      make_block_delay (52 + 143), // 24

      make_ap (171, 0.7), // 25
      make_ap (175, 0.7) // 26
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
    block_arr<sample> mm, lm, rm;
    block_arr<sample> lfo1, lfo2, tankmem;

    xspan m {mm.data(), io.size()};
    xspan l {lm.data(), io.size()};
    xspan r {rm.data(), io.size()};
    xspan tank {tankmem.data(), io.size()};

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto lfo         = tick_lfo<sample> (_lfo);
      auto mod         = fastgrowth (par.mod[i]);
      lfo1[i]          = sample {lfo[0] + lfo[1]} * mod * 0.5_r;
      lfo2[i]          = sample {lfo[2] + lfo[3]} * mod * 0.5_r;
      l[i]             = io[i][0] * 0.125_r;
      r[i]             = io[i][1] * 0.125_r;
      par.character[i] = fastgrowth (par.character[i]) * 0.5_r;
      par.decay[i]     = fastgrowth (par.decay[i]);
      par.decay[i]     = 0.5_r + par.decay[i] * 0.49_r;
    }

    run_cascade<0, 1> (_eng, l);
    run_cascade<2, 3> (_eng, r);

    _eng.run (sl<4> {}, l, blank, par.character);
    _eng.run (sl<5> {}, r, blank, par.character);
    span_visit (m, [&] (auto& v, uint i) { v = (l[i] + r[i]) * 0.5_r; });

    auto   lo   = upar.lf_amt;
    auto   hi   = upar.hf_amt;
    sample flo  = load_float<sample> (0.94f + lo * lo * 0.06f);
    sample glo  = load_float<sample> (0.85f + lo * 0.14f);
    sample fhi1 = load_float<sample> (0.48f - hi * hi * 0.3f);
    sample fhi2 = load_float<sample> (0.54f - hi * hi * 0.3f);
    sample fhi3 = load_float<sample> (0.52f - hi * hi * 0.34f);
    sample ghi  = load_float<sample> (0.55f + hi * hi * 0.45f);

    _eng.fetch_block (sl<24> {}, tank, 1);

    // block1
    _eng.run (sl<6> {}, tank, flo, glo, fhi1, 1_r, ghi);
    span_add (tank, m);
    auto g = m; // reuse m buffer, which is now free
    span_visit (g, [&] (auto& v, uint i) {
      v = -0.5_r * (1_r - (1_r - par.decay[i]) * 0.2_r);
    });
    _eng.run (sl<7> {}, tank, lfo1, g);
    span_mul (tank, par.decay);
    _eng.run (sl<8> {}, tank);
    _eng.run (sl<9, 10> {}, tank, blank, g);
    span_mul (tank, par.decay);
    _eng.run (sl<11> {}, tank);
    // block2
    _eng.run (sl<12> {}, tank, flo, glo, fhi2, 1_r, ghi);
    span_add (tank, l);
    span_mul (xspan {lfo1.data(), io.size()}, -1_r);
    span_visit (g, [&] (auto& v, uint i) {
      v = v * -0.95_r; // from max -0.5 to max 0.475)
    });
    _eng.run (sl<13> {}, tank, lfo1, g);
    span_mul (tank, par.decay);
    _eng.fetch (sl<14> {}, l, 748); // b2.tap_a
    _eng.run (sl<14> {}, tank);
    _eng.run (sl<15, 16> {}, tank, blank, g);
    span_mul (tank, par.decay);
    _eng.run (sl<17> {}, tank);
    // block3
    _eng.run (sl<18> {}, tank, flo, glo, fhi3, 1_r, ghi);
    span_add (tank, r);
    span_visit (g, [&] (auto& v, uint i) {
      v = v * -0.9473684210526316_r; // from max 0.475 to max -0.45)
    });
    _eng.run (sl<19> {}, tank, lfo2, g);
    span_mul (tank, par.decay);
    _eng.fetch (sl<20> {}, r, 1213); // b3.tap_a
    _eng.run (sl<20> {}, tank);
    span_mul (g, -1_r); // from -0.45 to 0.45)
    _eng.run (sl<21, 22> {}, tank, blank, g);
    span_mul (tank, par.decay);
    _eng.run (sl<23> {}, tank);
    _eng.push (sl<24> {}, tank.to_const());

    span_visit (xspan {par.character.data(), io.size()}, [] (auto& v, uint i) {
      // character is from 0 to 0.5 now. . 0.5 * 0.6 = 0.3
      v = 0.45_r + v * 0.6_r;
    });
    _eng.run (sl<25> {}, l, blank, par.character);
    _eng.run (sl<26> {}, r, blank, par.character);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
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
}}} // namespace artv::detail::tpaco
