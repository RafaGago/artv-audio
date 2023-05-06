#pragma once

#include "artv-common/dsp/own/fx/lofiverb/algorithm-engine.hpp"
#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
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
      par.decay[i] = 0.05_r + par.decay[i] * 0.92_r;
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
}}} // namespace artv::detail::lofiverb
