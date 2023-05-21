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
      auto lfo      = tick_lfo<sample> (_lfo);
      lfo1[i]       = sample {lfo[0]};
      lfo2[i]       = sample {lfo[2]};
      lfo1m[i]      = lfo1[i] * par.mod[i];
      lfo2m[i]      = lfo2[i] * par.mod[i];
      l_in[i]       = io[i][0] * 0.125_r * 0.25_r;
      r_in[i]       = io[i][1] * 0.125_r * 0.25_r;
      k1[i]         = 0.08_r + par.decay[i] * 0.1_r;
      k2[i]         = 0.08_r + par.decay[i] * 0.08_r;
      auto fg_decay = fastgrowth (par.decay[i]);
      k4[i]         = 0.15_r + fg_decay * 0.15_r + par.character[i] * 0.2_r;
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
}}} // namespace artv::detail::lofiverb
