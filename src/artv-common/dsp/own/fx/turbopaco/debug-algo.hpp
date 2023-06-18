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

//------------------------------------------------------------------------------
template <delay::data_type Dt>
class debug : public algorithm {
private:
  using engine = detail::tpaco::algo_engine<debug<Dt>, Dt, max_block_size>;

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
  void reset (xspan<u8> mem, float t_spl)
  {
    _eng.reset_memory (mem);
    _lfo.reset();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (1000.f),
      vec_set<f32_x2> (0.42f),
      vec_set<f32_x2> (3.85f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto f_er   = 0.3f + mod * 0.3f;
    auto f_late = 0.1f + mod * 1.2f;
    _lfo.set_freq (
      f32_x4 {f_er, f_er * 1.1234f, f_late * 1.1234f, f_late}, t_spl);
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
      make_ap (1787, 0.5, 261), // 7
      make_block_delay (max_block_size), // 8 to allow block processing
      // loop1
      make_ap (977, 0.5 /*overridden*/, 51), // 9
      make_delay (2819), // 10
      make_crossover2(), // 11
      make_ap (863, -0.5 /*overridden*/), // 12
      make_delay (1021), // 13
      make_ap (1453, 0.618), // 14
      make_block_delay (
        787), // 15 delay (allows block processing) (> blocksz + 1)
      // loop2
      make_ap (947, 0.5 /*overridden*/, 67), // 16
      make_delay (3191), // 17
      make_crossover2(), // 18
      make_ap (887, -0.5 /*overridden*/), // 19
      make_delay (1049), // 20
      make_ap (1367, 0.618), // 21
      make_hp (0.98), // 22
      make_block_delay (
        647)); // 23 delay (allows block processing) (> blocksz + 1)
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

    arr late_in_arr;
    arr lfo1;
    arr lfo2;
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (io[i][0] + io[i][1]) * 0.5_r;
      auto mod   = 0.25_r + (1_r - par.mod[i]) * 0.75_r;
      // ER + late lfo
      auto lfo = tick_lfo<sample> (_lfo);
      lfo1[i]  = sample {lfo[0]};
      lfo2[i]  = sample {lfo[1]} * (0.5_r + par.character[i] * 0.5_r);
      lfo3[i]  = sample {lfo[2]} * mod;
      lfo4[i]  = sample {lfo[3]} * mod;
      // decay fixup
      auto decay   = 1_r - par.decay[i];
      decay        = 1_r - decay * decay;
      par.decay[i] = 0.6_r + decay * 0.38_r;
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

    _eng.fetch_block (sl<8> {}, er2, 1); // feedback, fetching block + 1 samples

    span_visit (er1, [&] (auto& v, uint i) {
      v = late_in[i] * 0.5_r + er2[i];
    });
    er2.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<5> {}, er1, xspan {lfo2});
    _eng.run (sl<6> {}, er1, flo, glo, fhi, 1_r, ghi);
    xspan_memcpy (er1b, er1);
    span_mul (er1b, par.decay);
    _eng.run (sl<7> {}, er1b, xspan {lfo1});
    _eng.push (sl<8> {}, er1b.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto e_l = -er1[i] * 0.66_r - er2[i] * 0.34_r;
      auto e_r = -er1[i] * 0.66_r + er2[i] * 0.34_r;
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
    _eng.fetch_block (sl<15> {}, l, 1); // feedback, fetching block + 1 samples
    _eng.fetch_block (sl<23> {}, r, 1); // feedback, fetching block + 1 samples

    for (uint i = 0; i < io.size(); ++i) {
      auto loopsig = late_in[i] * 0.5_r + r[i];
      auto er_sig  = (er1[i] + er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = loopsig * (1_r - er_amt) + er_sig * er_amt;
      g[i]         = 0.618_r + par.character[i] * ((0.707_r - 0.618_r) * 2_r);
    }
    r.cut_head (1); // drop feedback sample from previous block

    _eng.run (sl<9> {}, late, xspan {lfo3}, g);
    _eng.run (sl<10> {}, late);
    _eng.run (sl<11> {}, late, flo, glo, fhi, 1_r, ghi);
    span_mul (late, par.decay);
    _eng.run (sl<12> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<13> {}, late);
    span_mul (late, par.decay);
    _eng.run (sl<14> {}, late);
    _eng.push (sl<15> {}, late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      auto loopsig = late_in[i] * 0.5_r + l[i];
      auto er_sig  = (er1[i] - er2[i]) * 0.25_r;
      auto er_amt  = par.character[i] * 0.5_r;
      late[i]      = loopsig * (1_r - er_amt) + er_sig * er_amt;
    }
    l.cut_head (1); // drop feedback sample from previous block
    _eng.run (sl<16> {}, late, xspan {lfo4}, g);
    _eng.run (sl<17> {}, late);
    _eng.run (sl<18> {}, late, flo, glo, fhi, 1_r, ghi);
    span_mul (late, par.decay);
    _eng.run (sl<19> {}, late, blank, [g] (uint i) { return -g[i]; });
    _eng.run (sl<20> {}, late);
    span_mul (late, par.decay);
    _eng.run (sl<21> {}, late);
    _eng.run (sl<22> {}, late);
    _eng.push (sl<23> {}, late.to_const()); // feedback point

    // Mixdown
    auto v1 = make_array (&l[0], &er1[0]);
    auto v2 = make_array<sample const*> (&r[0], &er2[0]);
    ep_crossfade<sample> (v1, v2, &par.character[0], io.size(), [=] (auto v) {
      return fastgrowth (v * 0.06_r);
    });

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
}}} // namespace artv::detail::tpaco
