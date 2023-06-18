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
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.0f, 0.f, 0.0f});
    _eq.reset_states_cascade();
    _eq.reset_coeffs (
      vec_set<f32_x2> (440.f),
      vec_set<f32_x2> (0.28f),
      vec_set<f32_x2> (2.f),
      t_spl,
      bell_tag {});
  }
  //----------------------------------------------------------------------------
  void mod_changed (float mod, float t_spl)
  {
    auto m1 = mod * mod;
    auto m2 = m1 + mod;
    auto f1 = 0.147f + m1 * 0.723f;
    auto f2 = 0.131f + m1 * 0.913f;
    auto f3 = 0.473f + m2 * 3.763f;
    auto f4 = 0.457f + m2 * 3.427f;
    _lfo.set_freq (f32_x4 {f1, f2, f3, f4}, t_spl);
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      make_lp(), // 0
      make_hp(), // 1
      make_ap (2331, 0.2, 407), // 2
      make_ap (2341, -0.2, 413), // 3
      make_ap (1967, -0.15, 509), // 4
      make_ap (1937, 0.15, 579), // 5
      make_ap (173, 0.5), // 6
      make_ap (177, 0.5) // 7
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
    using arr = block_arr<sample>;
    arr sig1, sig2, aux1, aux2, aux3, aux4;

    auto m = xspan {sig1.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      auto mod = 0.2_r + par.decay[i] * 0.8_r;
      mod *= (1_r - par.mod[i] * 0.8_r);
      mod *= mod;
      auto lfo = tick_lfo<sample> (_lfo);
      aux2[i]  = sample {lfo[0] - lfo[2]} * mod * 0.5_r;
      aux4[i]  = sample {lfo[1] + lfo[3]} * mod * 0.5_r;
      m[i]     = (io[i][0] + io[i][1]) * 0.5_r;
    }

    auto hi_g
      = load_float<sample> (0.1f + (1.f - upar.hf_amt * upar.hf_amt) * 0.65f);
    auto lo_g = load_float<sample> (0.8f + 0.18f * fastgrowth (upar.lf_amt));
    _eng.run (sl<0> {}, m, hi_g);
    _eng.run (sl<1> {}, m, lo_g);

    _eng.run (sl<2> {}, aux1, m.to_const(), aux2, [&] (uint i) {
      return par.character[i] * 0.2_r;
    });
    span_visit (xspan {aux2.data(), io.size()}, [&] (auto& v, uint i) {
      v = -v;
    });
    _eng.run (sl<3> {}, aux2, m.to_const(), aux2, [&] (uint i) {
      return par.character[i] * -0.2_r;
    });
    _eng.run (sl<4> {}, aux3, m.to_const(), aux4, [&] (uint i) {
      return par.character[i] * -0.15_r;
    });
    span_visit (xspan {aux4.data(), io.size()}, [&] (auto& v, uint i) {
      v = -v;
    });
    _eng.run (sl<5> {}, aux4, m.to_const(), aux4, [&] (uint i) {
      return par.character[i] * 0.15_r;
    });

    auto l = xspan {sig1.data(), io.size()};
    auto r = xspan {sig2.data(), io.size()};
    span_visit (l, [&] (auto&, uint i) {
      l[i] = aux1[i] + aux3[i];
      r[i] = aux4[i] - aux2[i];
    });

    _eng.run (sl<6> {}, l, blank, [&] (uint i) {
      return par.character[i] * 0.5_r;
    });
    _eng.run (sl<7> {}, r, blank, [&] (uint i) {
      return par.character[i] * -0.5_r;
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
}}} // namespace artv::detail::lofiverb
