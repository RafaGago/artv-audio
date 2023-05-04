#pragma once

#include "artv-common/dsp/own/fx/lofiverb/algorithm-engine.hpp"
#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
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
    _eq.reset_states_cascade();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f1 = 1.73f - mod * 0.63f;
    auto f2 = 1.53f - mod * 0.43f;
    _lfo.set_freq (f32_x4 {f1, f1, f1, f1}, t_spl);
  }
  //----------------------------------------------------------------------------
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
      auto lfo = tick_lfo<sample> (_lfo);
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
