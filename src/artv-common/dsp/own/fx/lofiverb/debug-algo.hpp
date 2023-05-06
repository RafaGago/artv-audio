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
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_lp (0.15), // 0
      make_hp (0.985), // 1

      make_comb (3821 - 1, 0.f, 54), // 2
      make_crossover2(), // 3
      make_parallel_delay (
        1, 429, -g, 1000, g, 1472, -g, 2088, g, 2765, -g, 3311, g), // 4

      make_comb (4036 + 1, 0.f, 53), // 5
      make_crossover2(), // 6
      make_parallel_delay (
        1,
        616,
        -g,
        1225,
        g,
        1691,
        -g,
        2434,
        g,
        3122,
        -g,
        3631,
        g), // 7

      make_comb (4059, 0.f, 44.f), // 8
      make_crossover2(), // 9
      make_parallel_delay (
        1,
        657,
        -g,
        1359,
        g,
        2184,
        -g,
        2744,
        g,
        3411,
        -g,
        3934,
        g), // 10

      make_ap (282, -0.7), // 11
      make_ap (343, -0.7), // 12
      // L
      make_delay (311, 311), // 13
      make_ap (233, -0.7), // 14
      make_ap (273, -0.7), // 15
      make_ap (534, -0.7), // 16
      // R
      make_delay (277, 400), // 17
      make_ap (194, -0.7), // 18
      make_ap (426, -0.7), // 19
      make_ap (566, -0.7)); // 20
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
    block_arr<sample> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan             in {in_mem.data(), io.size()};
    xspan             l {l_mem.data(), io.size()};
    xspan             r {r_mem.data(), io.size()};

    span_visit (in, [&] (auto& spl, uint i) {
      auto lfo = tick_lfo<sample> (_lfo);
      lfo1[i]  = sample {lfo[0]} * par.mod[i];
      lfo2[i]  = sample {lfo[1]} * par.mod[i];
      lfo3[i]  = sample {lfo[2]} * par.mod[i];
      spl      = (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<0> {}, in);
    _eng.run (sl<1> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = _eng.get_gain_for_rt60 (sl<2, 5, 8> {}, 0.25f + dec2 * 10.f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.88f + upar.lf_amt * 0.1f);
    sample fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample ghi = load_float<sample> (0.4f + dec2 * 0.4f + upar.hf_amt * 0.15f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<2> {}, comb_fb, lfo1, -gains[0]);
    _eng.run (sl<3> {}, comb_fb, flo, glo, fhi, 1_r, ghi);
    _eng.push (sl<2> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<4> {}, comb_fb.to_const(), overwrite, tank);

    _eng.fetch (sl<5> {}, comb_fb, lfo2, -gains[1]);
    _eng.run (sl<6> {}, comb_fb, flo, glo, fhi, 1_r, ghi);
    _eng.push (sl<5> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_add (comb_fb, in);
    span_visit (comb_fb, [&] (auto& s, uint i) { s += in[i]; });
    _eng.run (sl<7> {}, comb_fb.to_const(), add_to, tank);

    _eng.fetch (sl<8> {}, comb_fb, lfo3, -gains[2]);
    _eng.run (sl<9> {}, comb_fb, flo, glo, fhi, 1_r, ghi);
    _eng.push (sl<8> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_add (comb_fb, in);
    _eng.run (sl<10> {}, comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (c - 1_r * 0.5_r) * 2_r;
      c        = (1_r - c * c) * 0.4_r;
      eramt[i] = c;
    });

    xspan stank {tank.data(), io.size()};
    _eng.run (sl<11> {}, stank);
    _eng.run (sl<12> {}, l, stank.to_const());
    xspan_memdump (r.data(), l);

    xspan er {tmp2.data(), io.size()};
    _eng.run (sl<13> {}, er, in.to_const(), par.character);
    crossfade (l, er, eramt);
    _eng.run (sl<14> {}, l);
    _eng.run (sl<15> {}, l);
    _eng.run (sl<16> {}, l);

    _eng.run (sl<17> {}, er, in.to_const(), par.character);
    crossfade (r, er, eramt);
    _eng.run (sl<18> {}, r);
    _eng.run (sl<19> {}, r);
    _eng.run (sl<20> {}, r);

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
