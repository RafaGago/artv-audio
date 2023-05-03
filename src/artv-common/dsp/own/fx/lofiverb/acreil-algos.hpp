#pragma once

#include "artv-common/dsp/own/fx/lofiverb/algorithm.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
template <delay::data_type Dt>
class midifex49 : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<midifex49<Dt>, Dt, max_block_size>;

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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.2f + mod * 0.2f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // diffusors
      make_ap (321, 0.5), // 0 PreAP
      make_ap (431, 0.5), // 1 PreAP
      make_ap (968, 0.5), // 2 PreAP
      make_ap (1620, 0.5), // 3 PreAP

      make_delay (21), // 4 L1
      make_block_delay (1010, 0), // 5 R1

      make_delay (1624), // 6 Loop
      make_ap (1992, 0.5, 17), // 7 Loop

      make_block_delay (1891, 0), // 8 L2
      make_block_delay (890, 0), // 9 R2

      make_delay (2110), // 10 Loop
      make_quantizer(), // 11
      make_ap (2371, 0.5), // 12 Loop nested allpass 1
      make_ap (1378, 0.2), // 13 Loop nested allpass 2

      make_block_delay (2003, 0), // 14 L3
      make_block_delay (671, 0), // 15 R3

      make_delay (2157 - max_block_size), // 16 Loop
      make_quantizer(), // 17
      make_crossover2(), // 18
      make_ap (2712, 0.5, 22), // 19 Loop nested allpass 1
      make_ap (1783, 0.2), // 20 Loop nested allpass 2

      make_block_delay (max_block_size) // 22 delay block (== blocksz + 1)
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr = block_arr<sample>;

    arr   tmp_arr;
    xspan tmp {tmp_arr.data(), io.size()};
    arr   sig_arr;
    xspan sig {sig_arr.data(), io.size()};
    arr   l_arr;
    xspan l {l_arr.data(), io.size()};
    arr   r_arr;
    xspan r {r_arr.data(), io.size()};
    arr   lfo1_arr;
    xspan lfo1 {lfo1_arr.data(), io.size()};
    arr   lfo2_arr;
    xspan lfo2 {lfo2_arr.data(), io.size()};

    span_visit (io, [&] (auto& spl, uint i) {
      // Midside signal
      sig[i]
        = (sample) (((spl[0] + spl[1]) * 0.25_r)); // gain = 0.5 + feedback = 1
      // LFO
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      // decay fixup
      auto decay   = (sample) (_eng.one - par.decay[i]);
      decay        = (sample) (_eng.one - decay * decay);
      par.decay[i] = (sample) (0.1_r + decay * 0.8375_r);
    });

    _eng.run (sl<0> {}, sig);
    _eng.run (sl<1> {}, sig);
    _eng.run (sl<2> {}, sig);
    _eng.run (sl<3> {}, sig);

    // feedback handling, fetching the block with a negative offset of 1
    _eng.fetch_block (sl<21> {}, tmp, 1);
    span_visit (sig, [&] (auto& v, uint i) { v = (sample) (v + tmp[i]); });

    // 1st output point for L and R signal
    xspan_memcpy (l, sig); // delay LT a block -> might overlap, requires copy
    _eng.run (sl<4> {}, l);
    _eng.fetch_block (sl<5> {}, r, 0); // delay GT a block will never overlap
    _eng.push (sl<5> {}, sig.to_const());

    // continuing the loop
    _eng.run (sl<6> {}, sig);
    _eng.run (sl<7> {}, sig, lfo1, blank);

    // 2nd output point for L and R signal
    _eng.fetch_block (sl<8> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (sample) ((v + tmp[i]) * (2_r / 3_r));
    });
    _eng.push (sl<8> {}, sig.to_const());
    _eng.fetch_block (sl<9> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (sample) ((v + tmp[i]) * (2_r / 3_r));
    });
    _eng.push (sl<9> {}, sig.to_const());

    // continuing the loop
    _eng.run (sl<10> {}, sig);
    _eng.run (sl<11> {}, sig, [&] (auto v, uint i) {
      return v * par.decay[i];
    });
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (sample) (par.character[i] * 0.14_r);
    });
    _eng.run (sl<12, 13> {}, sig, blank, blank, blank, tmp);

    // 3rd output point for L and R signal
    _eng.fetch_block (sl<14> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (l, [&] (auto& v, uint i) {
      v = (sample) (v + tmp[i] * (1_r / 3_r));
    });
    _eng.push (sl<14> {}, sig.to_const());
    _eng.fetch_block (sl<15> {}, tmp, 0); // delay GT a block will never overlap
    span_visit (r, [&] (auto& v, uint i) {
      v = (sample) (v + tmp[i] * (1_r / 3_r));
    });
    _eng.push (sl<15> {}, sig.to_const());

    // continuing the loop
    _eng.run (sl<16> {}, sig);
    _eng.run (sl<17> {}, sig, [&] (auto v, uint i) {
      return v * par.decay[i];
    });

    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.35f + upar.lf_amt * 0.6f);
    sample fhi = load_float<sample> (0.8f - upar.hf_amt * upar.hf_amt * 0.1f);
    sample ghi = load_float<sample> (0.05f + upar.hf_amt * 0.6f);

    _eng.run (sl<18> {}, sig, flo, glo, fhi, _eng.one, ghi);
    span_visit (tmp, [&] (auto& v, uint i) {
      v = (sample) (par.character[i] * 0.2_r);
    });
    _eng.run (sl<19, 20> {}, sig, lfo2, blank, blank, tmp);
    // push to delay feedback
    _eng.push (sl<21> {}, sig.to_const());

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }

private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class midifex50 : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<midifex50<Dt>, Dt, max_block_size>;

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
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.3f + mod * 0.1f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    return make_array<stage_data> (
      // outputs first
      make_quantizer(), // 0
      make_ap (147, 0.5), // 1 L diff
      make_ap (43, 0.5), // 2 L diff
      make_ap (55, 0.5), // 3 L diff
      make_ap (249, 0.5), // 4 R diff
      make_ap (48, 0.5), // 5 R diff
      make_ap (21, 0.5), // 6 R diff
      // loop
      make_ap (13, 0.5), // 7
      make_ap (83, 0.5), // 8
      make_ap (116, 0.5), // 9
      make_ap (239, 0.5), // 10
      make_ap (339, 0.5), // 11
      make_ap (481, 0.5), // 12
      make_ap (555, 0.5), // 13
      make_ap (823, 0.5), // 14
      make_ap (999, 0.5), // 15
      make_ap (1100, 0.5), // 16
      make_ap (1347, 0.5), // 17
      make_ap (1563, 0.5), // 18
      make_ap (1841 - 32, 0.5), // 19
      make_block_delay (32), // 20 C1
      make_ap (2001 - 32, 0.5, 67), // 21
      make_block_delay (32), // 22 C2
      make_ap (2083 - 32, 0.5, 127), // 23
      make_crossover2(), // 24
      make_block_delay (max_block_size) // 25 (FB point) (== blocksz + 1)
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint)
  {
    using arr    = block_arr<sample>;
    using arr_fb = fb_block_arr<sample>;

    arr_fb tmp1;
    xspan  loop {tmp1.data(), io.size() + 1};
    arr    l_arr;
    xspan  l {l_arr.data(), io.size()};
    arr    r_arr;
    xspan  r {r_arr.data(), io.size()};

    // xspan_memcpy (loop, io);
    _eng.fetch_block (sl<25> {}, loop, 1);
    xspan_memcpy (r, loop.advanced (1));
    xspan_memcpy (l, loop.advanced (1));

    loop.cut_tail (1); // feedback samples
    _eng.run (sl<0> {}, loop, [&] (sample fb_spl, uint i) {
      // decay fixup
      auto decay = (sample) (_eng.one - par.decay[i]);
      decay      = (sample) (_eng.one - decay * decay);
      decay      = (sample) - (0.3_r + decay * 0.45_r);
      return (fb_spl * decay) + ((io[i][0] + io[i][1]) * 0.25_r); // gain = 1
    }); // feedback + input summing with quantizer

    _eng.run (sl<1> {}, l);
    _eng.run (sl<2> {}, l);
    _eng.run (sl<3> {}, l);

    _eng.run (sl<4> {}, r);
    _eng.run (sl<5> {}, r);
    _eng.run (sl<6> {}, r);

    _eng.run (sl<7> {}, loop);
    _eng.run (sl<8> {}, loop);
    _eng.run (sl<9> {}, loop);
    _eng.run (sl<10> {}, loop);
    _eng.run (sl<11> {}, loop);
    _eng.run (sl<12> {}, loop);
    _eng.run (sl<13> {}, loop);
    _eng.run (sl<14> {}, loop);
    _eng.run (sl<15> {}, loop);
    _eng.run (sl<16> {}, loop);
    _eng.run (sl<17> {}, loop);
    _eng.run (sl<18> {}, loop);
    _eng.run (sl<19> {}, loop);

    arr tmp2;
    arr tmp3;

    xspan c1 {tmp2.data(), loop.size()};
    xspan lfo2 {tmp3.data(), loop.size()};
    _eng.fetch_block (sl<20> {}, c1, 0);
    _eng.push (sl<20> {}, loop.to_const());
    auto lfo1 = loop; // data inserted already (lfo1 -> tmp1)
    loop      = c1; // avoid a copy. loop -> tmp2
    span_visit (l, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = _eng.one - c;
      v      = (sample) (k * v + loop[i] * c); // L done
      // unrelated but done here to skip one iteration
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
    });
    _eng.run (sl<21> {}, loop, xspan {lfo1}, blank); // tmp1 free
    xspan c2 {tmp1.data(), loop.size()};
    _eng.fetch_block (sl<22> {}, c2, 0);
    _eng.push (sl<22> {}, loop.to_const());
    loop = c2; // avoid copy. tmp2 (loop) free.
    span_visit (r, [&] (auto& v, uint i) {
      auto c = par.character[i] * 0.5_r;
      auto k = _eng.one - c;
      v      = (sample) (k * v + loop[i] * c); // R done
    });
    // outputs done, dump now that they have been recently touched
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
    _eng.run (sl<23> {}, loop, xspan {lfo2}, blank); // tmp3 free

    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.35f + upar.lf_amt * 0.6f);
    sample fhi = load_float<sample> (0.62f - upar.hf_amt * upar.hf_amt * 0.2f);
    sample ghi = load_float<sample> (0.35f + upar.hf_amt * 0.3f);

    _eng.run (sl<24> {}, loop, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<25> {}, loop.to_const());
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class dre2000a : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<dre2000a<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.75f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.12), // 1
      make_hp (0.99), // 2
      make_parallel_delay (1, 42, -g, 586, g, 1099, -g), // 3 (L)
      make_parallel_delay (1, 105, g, 490, -g, 1290, g), // 4 (R)

      make_comb (5719 - 1, 0.f, 43), // 5
      make_crossover2(), // 6
      make_parallel_delay (
        1, 640, -g, 1494, g, 2199, -g, 3122, g, 4135, -g, 4952, g), // 7

      make_comb (5779 - 1, 0.f, 44), // 8
      make_crossover2(), // 9
      make_parallel_delay (
        1,
        902,
        -g,
        1830,
        g,
        2528,
        -g,
        3641,
        g,
        4670,
        -g,
        5432,
        g), // 10

      make_comb (5905 - 1, 0.f, 42), // 11
      make_crossover2(), // 12
      make_parallel_delay (
        1,
        979,
        -g,
        2031,
        g,
        3267,
        -g,
        4104,
        g,
        5103,
        -g,
        5886,
        g), // 13
      // L
      make_ap (224, -0.7), // 14
      make_delay (311, 311), // 15
      make_ap (227, -0.7), // 16
      make_ap (1343, -0.7), // 17
      // R
      make_ap (224, -0.7), // 18
      make_delay (277, 400), // 19
      make_ap (91, -0.7), // 20
      make_ap (1182, -0.7)); // 21
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<sample> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan             in {in_mem.data(), io.size()};
    xspan             l {l_mem.data(), io.size()};
    xspan             r {r_mem.data(), io.size()};

    _eng.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      lfo3[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<1> {}, in);
    _eng.run (sl<2> {}, l, in.to_const());
    xspan_memdump (r.data(), l);
    _eng.run (sl<3> {}, l, overwrite);
    _eng.run (sl<4> {}, r, overwrite);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = _eng.get_gain_for_rt60 (sl<5, 8, 11> {}, 0.9f + dec2 * 19.1f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.8f + upar.lf_amt * 0.18f);
    sample fhi = load_float<sample> (0.7f - upar.hf_amt * upar.hf_amt * 0.05f);
    sample ghi = load_float<sample> (dec2 * 0.25f + upar.hf_amt * 0.45f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<5> {}, comb_fb, lfo1, -gains[0]);
    _eng.run (sl<6> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<5> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<7> {}, comb_fb.to_const(), overwrite, tank);

    _eng.fetch (sl<8> {}, comb_fb, lfo2, -gains[1]);
    _eng.run (sl<9> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<8> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<10> {}, comb_fb.to_const(), add_to, tank);

    _eng.fetch (sl<11> {}, comb_fb, lfo3, -gains[2]);
    _eng.run (sl<12> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<11> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<13> {}, comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      l[i]     = (sample) (l[i] + tank[i]);
      r[i]     = (sample) (r[i] + tank[i]);
      c        = (sample) ((c - _eng.one * 0.5_r) * 2_r);
      c        = (sample) ((_eng.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    auto stank = xspan {tank.data(), io.size()}.to_const();
    _eng.run (sl<14> {}, l, stank);
    xspan er {tmp2.data(), io.size()};
    _eng.run (sl<15> {}, er, in.to_const(), par.character);
    crossfade (l, er, eramt, _eng.one);
    _eng.run (sl<16> {}, l);
    _eng.run (sl<17> {}, l);

    _eng.run (sl<18> {}, r, stank);
    _eng.run (sl<19> {}, er, in.to_const(), par.character);
    crossfade (r, er, eramt, _eng.one);
    _eng.run (sl<20> {}, r);
    _eng.run (sl<21> {}, r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class dre2000b : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<dre2000b<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.15f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.3333f, 0.6666f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.15), // 1
      make_hp (0.985), // 2

      make_comb (3821 - 1, 0.f, 54), // 3
      make_crossover2(), // 4
      make_parallel_delay (
        1, 429, -g, 1000, g, 1472, -g, 2088, g, 2765, -g, 3311, g), // 5

      make_comb (4036 - 1, 0.f, 53), // 6
      make_crossover2(), // 7
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
        g), // 8

      make_comb (4059 - 1, 0.f, 44.f), // 9
      make_crossover2(), // 10
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
        g), // 11

      make_ap (282, -0.7), // 12
      make_ap (343, -0.7), // 13
      // L
      make_delay (311, 311), // 14
      make_ap (233, -0.7), // 15
      make_ap (273, -0.7), // 16
      make_ap (534, -0.7), // 17
      // R
      make_delay (277, 400), // 18
      make_ap (194, -0.7), // 19
      make_ap (426, -0.7), // 20
      make_ap (566, -0.7)); // 21
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<sample> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, tmp1, tmp2, tank;
    xspan             in {in_mem.data(), io.size()};
    xspan             l {l_mem.data(), io.size()};
    xspan             r {r_mem.data(), io.size()};

    _eng.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      lfo3[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<1> {}, in);
    _eng.run (sl<2> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = _eng.get_gain_for_rt60 (sl<3, 6, 9> {}, 0.25f + dec2 * 10.f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.88f + upar.lf_amt * 0.1f);
    sample fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.4f);
    sample ghi = load_float<sample> (0.4f + dec2 * 0.4f + upar.hf_amt * 0.15f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<3> {}, comb_fb, lfo1, -gains[0]);
    _eng.run (sl<4> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<3> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<5> {}, comb_fb.to_const(), overwrite, tank);

    _eng.fetch (sl<6> {}, comb_fb, lfo2, -gains[1]);
    _eng.run (sl<7> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (sample) (s + in[i]); });
    _eng.run (sl<8> {}, comb_fb.to_const(), add_to, tank);

    _eng.fetch (sl<9> {}, comb_fb, lfo3, -gains[2]);
    _eng.run (sl<10> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    span_visit (comb_fb, [&] (auto& s, uint i) { s = (sample) (s + in[i]); });
    _eng.run (sl<11> {}, comb_fb.to_const(), add_to, tank);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (sample) ((c - _eng.one * 0.5_r) * 2_r);
      c        = (sample) ((_eng.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    xspan stank {tank.data(), io.size()};
    _eng.run (sl<12> {}, stank);
    _eng.run (sl<13> {}, l, stank.to_const());
    xspan_memdump (r.data(), l);

    xspan er {tmp2.data(), io.size()};
    _eng.run (sl<14> {}, er, in.to_const(), par.character);
    crossfade (l, er, eramt, _eng.one);
    _eng.run (sl<15> {}, l);
    _eng.run (sl<16> {}, l);
    _eng.run (sl<17> {}, l);

    _eng.run (sl<18> {}, er, in.to_const(), par.character);
    crossfade (r, er, eramt, _eng.one);
    _eng.run (sl<19> {}, r);
    _eng.run (sl<20> {}, r);
    _eng.run (sl<21> {}, r);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class dre2000c : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<dre2000c<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.75f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //------------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.18), // 1
      make_hp (0.985), // 2
      make_ap (115, -0.7f), // 3
      make_ap (160, -0.7f), // 4
      make_ap (231, -0.7f), // 5

      make_comb (3794 - 1, 0.f, 23), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 20 + 13, -g, 959, g, 1817, -g, 2855, g), // 8 R/L

      make_comb (3838 - 1, 0.f, 24), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 339, -g, 1309, g, 2271, -g, 3221, g), // 11 R/L

      make_comb (3861 - 1, 0.f, 25), // 12
      make_crossover2(), // 13
      make_parallel_delay (2, 130, -g, 1104, g, 2065, -g, 3391, g), // 14 R/L

      make_comb (3894 - 1, 0.f, 26), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 499, -g, 1445, g, 2071, -g, 2885, g), // 17 R/L
      // L
      make_ap (140, -0.7), // 18
      make_ap (160, -0.7), // 19
      make_delay (212, 237), // 20
      // R
      make_ap (160, -0.7), // 21
      make_ap (227, -0.7), // 22
      make_delay (213, 126) // 23
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<sample> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1;
    xspan             in {in_mem.data(), io.size()};
    xspan             l {l_mem.data(), io.size()};
    xspan             r {r_mem.data(), io.size()};

    _eng.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      lfo3[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      lfo4[i]  = (sample) (sample {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<1> {}, in);
    _eng.run (sl<2> {}, in);
    _eng.run (sl<3> {}, in);
    _eng.run (sl<4> {}, in);
    _eng.run (sl<5> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains = _eng.get_gain_for_rt60 (
      sl<6, 9, 12, 15> {}, 0.25f + dec2 * 20.f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.88f + upar.lf_amt * 0.1f);
    sample fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    sample ghi
      = load_float<sample> (0.3f + dec2 * 0.25f + upar.hf_amt * 0.4499f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<6> {}, comb_fb, lfo1, gains[0]);
    _eng.run (sl<7> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<8> {}, comb_fb.to_const(), overwrite, r, l);

    _eng.fetch (sl<9> {}, comb_fb, lfo2, gains[1]);
    _eng.run (sl<10> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<11> {}, comb_fb.to_const(), add_to, r, l);

    _eng.fetch (sl<12> {}, comb_fb, lfo3, gains[2]);
    _eng.run (sl<13> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<12> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<14> {}, comb_fb.to_const(), add_to, r, l);

    _eng.fetch (sl<15> {}, comb_fb, lfo4, gains[3]);
    _eng.run (sl<16> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<15> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<17> {}, comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (sample) ((c - _eng.one * 0.5_r) * 2_r);
      c        = (sample) ((_eng.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    _eng.run (sl<18> {}, l);
    _eng.run (sl<19> {}, l);
    _eng.run (sl<20> {}, lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt, _eng.one);

    _eng.run (sl<21> {}, r);
    _eng.run (sl<22> {}, r);
    _eng.run (sl<23> {}, in, par.character);
    crossfade (r, in, eramt, _eng.one);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class dre2000d : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<dre2000d<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (16800, 21000, 25200, 32400, 40320, 57600, 57600);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.35f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.185), // 1
      make_hp (0.99), // 2
      make_ap (115, -0.7f), // 3
      make_ap (214, -0.7f), // 4
      make_ap (365, -0.7f), // 5

      make_comb (5652 - 1, 0.f, 57), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 745, -g, 2164, g, 2922, -g, 5005, g), // 8 R/L

      make_comb (5679 - 1, 0.f), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 625, -g, 2080, g, 3523, -g, 4948, g), // 11 R/L

      make_comb (5689 - 1, 0.f), // 12
      make_crossover2(), // 13
      make_parallel_delay (
        2, 11 + 21, -g, 1435, g, 2769, -g, 4279, g), // 14 R/L

      make_comb (5711 - 1, 0.f, 57), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 59, -g, 1557, g, 3539, -g, 4304, g), // 17 R/L
      // L
      make_ap (214, -0.7, 30), // 18
      make_ap (469, -0.7), // 19
      make_delay (234, 337), // 20
      // R
      make_ap (245, -0.7, 29), // 21
      make_ap (426, -0.7), // 22
      make_delay (434, 247) // 23
    );
  }
  //----------------------------------------------------------------------------

  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<sample> in_mem, l_mem, r_mem, lfo1, lfo2, lfo3, lfo4, tmp1, tmp2;
    xspan             in {in_mem.data(), io.size()};
    xspan             l {l_mem.data(), io.size()};
    xspan             r {r_mem.data(), io.size()};

    _eng.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      lfo3[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      lfo4[i]  = (sample) (sample {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<1> {}, in);
    _eng.run (sl<2> {}, in);
    _eng.run (sl<3> {}, in);
    _eng.run (sl<4> {}, in);
    _eng.run (sl<5> {}, in);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains
      = _eng.get_gain_for_rt60 (sl<6, 9, 12, 15> {}, 0.25f + dec2 * 5.f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.80f + upar.lf_amt * 0.199f);
    sample fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    sample ghi = load_float<sample> (0.4f + dec2 * 0.1f + upar.hf_amt * 0.49f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<6> {}, comb_fb, lfo1, gains[0]);
    _eng.run (sl<7> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<6> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<8> {}, comb_fb.to_const(), overwrite, r, l);

    _eng.fetch (sl<9> {}, comb_fb, blank, gains[1]);
    _eng.run (sl<10> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<9> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<11> {}, comb_fb.to_const(), add_to, r, l);

    _eng.fetch (sl<12> {}, comb_fb, blank, gains[2]);
    _eng.run (sl<13> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<12> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<14> {}, comb_fb.to_const(), add_to, r, l);

    _eng.fetch (sl<15> {}, comb_fb, lfo4, gains[3]);
    _eng.run (sl<16> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<15> {}, comb_fb, comb_fb.to_const(), in.to_const());
    _eng.run (sl<17> {}, comb_fb.to_const(), add_to, r, l);
    xspan eramt {tmp1.data(), io.size()};
    span_visit (xspan {par.character.data(), io.size()}, [&] (auto c, uint i) {
      c        = (sample) ((c - _eng.one * 0.5_r) * 2_r);
      c        = (sample) ((_eng.one - c * c) * 0.4_r);
      eramt[i] = c;
    });

    _eng.run (sl<18> {}, l, lfo2);
    _eng.run (sl<19> {}, l);
    _eng.run (sl<20> {}, lfo1, in.to_const(), par.character);
    crossfade (l, lfo1, eramt, _eng.one);

    _eng.run (sl<21> {}, r, lfo3);
    _eng.run (sl<22> {}, r);
    _eng.run (sl<23> {}, in, par.character);
    crossfade (r, in, eramt, _eng.one);

    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------
template <delay::data_type Dt>
class rev5_l_hall : public algorithm {
private:
  using engine
    = detail::lofiverb::algo_engine<rev5_l_hall<Dt>, Dt, max_block_size>;

public:
  //----------------------------------------------------------------------------
  static constexpr auto get_sample_rates()
  {
    return make_array (25200, 33600, 40320, 44100, 45000, 52500, 57600);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    return engine::get_required_bytes();
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<u8> mem) { _eng.reset_memory (mem); }
  //----------------------------------------------------------------------------
  static void reset_lfo_freq (lfo<4>& lfo, float mod, float t_spl)
  {
    auto f = 0.25f + mod * 0.35f;
    lfo.set_freq (f32_x4 {f, f, f, f}, t_spl);
  }
  //----------------------------------------------------------------------------
  static void reset_lfo_phase (lfo<4>& lfo)
  {
    lfo.set_phase (
      phase<4> {phase_tag::normalized {}, 0.f, 0.25f, 0.5f, 0.75f});
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_spec()
  {
    constexpr float g = 0.125f;
    return make_array<stage_data> (
      make_quantizer(), // 0
      make_lp (0.18), // 1
      make_hp (0.987), // 2

      make_quantizer(), // 3
      make_ap (1257, 0.75f), // 4
      make_ap (978, 0.75f), // 5

      // rev main
      make_comb (2076, 0., 43), // 6
      make_crossover2(), // 7
      make_parallel_delay (2, 949, g, 837, g), // 8
      make_comb (2894, 0., 43), // 9
      make_crossover2(), // 10
      make_parallel_delay (2, 1330, g, 1412, g), // 11
      make_comb (3295, 0.), // 12
      make_crossover2(), // 13
      make_parallel_delay (2, 2393, g, 2704, g), // 14
      make_comb (3919, 0.), // 15
      make_crossover2(), // 16
      make_parallel_delay (2, 3263, g, 3011, g), // 17
      make_comb (4570, 0.), // 18
      make_crossover2(), // 19
      make_parallel_delay (2, 3667, g, 3667, g), // 20
      // sub reverb cond
      make_variable_delay (533, 533 + 2411), // 21
      make_ap (1212, 0.75f), // 22
      make_ap (960, 0.75f), // 23
      // sub reverb
      make_comb (3169, 0., 44), // 24
      make_crossover2(), // 25
      make_parallel_delay (2, 1757, g, 1644, g), // 26
      make_comb (3753, 0., 44), // 27
      make_crossover2(), // 28
      make_parallel_delay (2, 1900, g, 1981, g), // 29
      make_comb (4280, 0.), // 30
      make_crossover2(), // 31
      make_parallel_delay (2, 2838, g, 3148, g), // 32
      make_comb (4491, 0.), // 33
      make_crossover2(), // 34
      make_parallel_delay (2, 3798, g, 3545, g), // 35
      make_comb (5091, 0.), // 36
      make_crossover2(), // 37
      make_parallel_delay (2, 4298, g, 4298, g), // 38
      make_ap (682, 0.75f), // 39
      make_ap (830, 0.75f), // 40
      make_ap (695, 0.75f), // 41
      make_ap (844, 0.75f), // 42
      make_parallel_delay (
        2,
        18 + 14,
        0.95f,
        18 + 14,
        0.95f,
        895,
        0.95f,
        695,
        0.95f,
        1238,
        0.75f,
        1487,
        0.75f,
        1870,
        0.75f,
        1936,
        0.75f,
        2724,
        0.5f,
        2635,
        0.5f,
        3365,
        0.25f,
        3609,
        0.25f) // 43 L/R
    );
  }
  //----------------------------------------------------------------------------
  using value_type = typename engine::value_type;
  using sample     = value_type;
  //----------------------------------------------------------------------------
  void process_block (
    lfo<4>&                      lfo_obj,
    xspan<std::array<sample, 2>> io,
    smoothed_parameters<sample>& par,
    unsmoothed_parameters const& upar,
    uint                         srate)
  {
    block_arr<sample> in_mem, main_mem, sub_mem, l_mem, r_mem, lfo1, lfo2, lfo3,
      lfo4, tmp1, tmp2;
    xspan main {main_mem.data(), io.size()};
    xspan sub {sub_mem.data(), io.size()};
    xspan in {in_mem.data(), io.size()};
    xspan l {l_mem.data(), io.size()};
    xspan r {r_mem.data(), io.size()};

    _eng.run (sl<0> {}, in, [&] (auto spl, uint i) {
      auto lfo = tick_lfo<sample> (lfo_obj);
      lfo1[i]  = (sample) (sample {lfo[0]} * par.mod[i]);
      lfo2[i]  = (sample) (sample {lfo[1]} * par.mod[i]);
      lfo3[i]  = (sample) (sample {lfo[2]} * par.mod[i]);
      lfo4[i]  = (sample) (sample {lfo[3]} * par.mod[i]);
      return (io[i][0] + io[i][1]) * 0.25_r;
    });
    _eng.run (sl<1> {}, in);
    _eng.run (sl<2> {}, in);
    xspan_memdump (main.data(), in);
    xspan_memdump (sub.data(), in);
    _eng.run (sl<3> {}, main, [&] (auto v, uint i) {
      auto c  = par.character[i];
      tmp1[i] = (sample) (0.5_r + 0.2_r * c); // cg3
      tmp2[i] = (sample) (0.35_r + 0.7_r * clamp (c, 0._r, 0.4999_r)); // cg2
      return v * 0.25_r;
    });
    _eng.run (sl<4> {}, main, blank, tmp1);
    _eng.run (sl<5> {}, main, blank, tmp2);

    float dec2 = as_float (par.decay[0]);
    dec2 *= dec2;
    auto gains = _eng.get_gain_for_rt60 (
      sl<6, 9, 12, 15, 18> {}, 1.f + dec2 * 10.f, srate);
    sample flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    sample glo = load_float<sample> (0.9f + upar.lf_amt * 0.09f);
    sample fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.1f);
    sample ghi = load_float<sample> (0.75f + upar.hf_amt * 0.2f);

    xspan comb_fb {tmp1.data(), io.size()};
    _eng.fetch (sl<6> {}, comb_fb, lfo1, gains[0]);
    _eng.run (
      sl<7> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.85_r));
    _eng.push (sl<6> {}, comb_fb, comb_fb.to_const(), main.to_const());
    _eng.run (sl<8> {}, comb_fb.to_const(), overwrite, l, r);

    _eng.fetch (sl<9> {}, comb_fb, lfo3, gains[1]);
    _eng.run (
      sl<10> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.9_r));
    _eng.push (sl<9> {}, comb_fb, comb_fb.to_const(), main.to_const());
    _eng.run (sl<11> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<12> {}, comb_fb, blank, gains[2]);
    _eng.run (
      sl<13> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.95_r));
    _eng.push (sl<12> {}, comb_fb, comb_fb.to_const(), main.to_const());
    _eng.run (sl<14> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<15> {}, comb_fb, blank, gains[3]);
    _eng.run (
      sl<16> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.97_r));
    _eng.push (sl<15> {}, comb_fb, comb_fb.to_const(), main.to_const());
    _eng.run (sl<17> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<18> {}, comb_fb, blank, gains[4]);
    _eng.run (sl<19> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<18> {}, comb_fb, comb_fb.to_const(), main.to_const());
    _eng.run (sl<20> {}, comb_fb.to_const(), add_to, l, r);

    _eng.run (sl<21> {}, sub, par.character);
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto c  = par.character[i];
      tmp1[i] = (sample) (0.5_r + 0.2_r * c); // cg3
      // tmp2 still not overwritten...
      // tmp2[i] = (sample) (0.35_r + 0.999_r * clamp (c, 0._r, 0.333_r)); //
      // cg2
    }
    _eng.run (sl<22> {}, main, tmp2);
    _eng.run (sl<23> {}, main, tmp1);
    gains = _eng.get_gain_for_rt60 (
      sl<24, 27, 30, 33, 36> {}, 1.f + dec2 * 6.f, srate);
    flo = load_float<sample> (0.9f + upar.lf_amt * upar.lf_amt * 0.05f);
    glo = load_float<sample> (0.9f + upar.lf_amt * 0.04f);
    fhi = load_float<sample> (0.82f - upar.hf_amt * upar.hf_amt * 0.2f);
    ghi = load_float<sample> (0.7f + upar.hf_amt * 0.25f);

    // sub  reverb
    _eng.fetch (sl<24> {}, comb_fb, lfo2, gains[0]);
    _eng.run (
      sl<25> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.9_r));
    _eng.push (sl<24> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    _eng.run (sl<26> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<27> {}, comb_fb, lfo4, gains[1]);
    _eng.run (
      sl<28> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.93_r));
    _eng.push (sl<27> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    _eng.run (sl<29> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<30> {}, comb_fb, blank, gains[2]);
    _eng.run (
      sl<31> {}, comb_fb, flo, glo, fhi, _eng.one, (sample) (ghi * 0.97_r));
    _eng.push (sl<30> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    _eng.run (sl<32> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<33> {}, comb_fb, blank, gains[3]);
    _eng.run (sl<34> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<33> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    _eng.run (sl<35> {}, comb_fb.to_const(), add_to, l, r);

    _eng.fetch (sl<36> {}, comb_fb, blank, gains[4]);
    _eng.run (sl<37> {}, comb_fb, flo, glo, fhi, _eng.one, ghi);
    _eng.push (sl<36> {}, comb_fb, comb_fb.to_const(), sub.to_const());
    _eng.run (sl<38> {}, comb_fb.to_const(), add_to, l, r);

    _eng.run (sl<39> {}, l, blank, [&] (uint i) {
      return (sample) (0.65_r + 0.1_r * lfo1[i]);
    });
    _eng.run (sl<40> {}, l, blank, [&] (uint i) {
      auto c  = par.character[i];
      tmp1[i] = (sample) (0.325_r + 0.999_r * clamp (c, 0._r, 0.333_r)); // cg1
      return (sample) (tmp1[i] + 0.05_r * lfo2[i]);
    });
    _eng.run (sl<41> {}, r, blank, [&] (uint i) {
      return (sample) (0.65_r + 0.1_r * lfo3[i]);
    });
    _eng.run (sl<42> {}, r, blank, [&] (uint i) {
      return (sample) (tmp1[i] + 0.05_r * lfo4[i]);
    });

    _eng.run (sl<43> {}, in.to_const(), add_to, l, r); // ER
    span_visit (io, [&] (auto& spls, uint i) {
      spls[0] = l[i];
      spls[1] = r[i];
    });
  }
  //----------------------------------------------------------------------------
private:
  engine _eng;
};
//------------------------------------------------------------------------------

}}} // namespace artv::detail::lofiverb
