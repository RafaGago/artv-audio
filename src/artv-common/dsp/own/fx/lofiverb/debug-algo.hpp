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
      l_in[i]       = io[i][0] * 0.5_r;
      r_in[i]       = io[i][1] * 0.5_r;
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

    // TODO: ER
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
}}} // namespace artv::detail::lofiverb
