#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
namespace detail { namespace lofiverb {
//------------------------------------------------------------------------------
using q0_15          = fixpt_p<1, 0, 15>;
using q1_14          = fixpt_p<1, 1, 14>;
using fixpt_spls     = fixpt_p<0, 14, 0>;
using fixpt_spls_mod = fixpt_p<0, 9, 0>;
//------------------------------------------------------------------------------
struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // in samples
  q0_15          g;
};
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_dd (
  u16   spls,
  float g   = 0.,
  u16   mod = 0)
{
  assert (spls < 8 * 1024);
  return delay_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    q0_15::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr auto get_rev1_delays()
{
  return make_array<delay_data> (
    // diffusors
    make_dd (147, -0.707), // AP
    make_dd (183, 0.707), // AP
    make_dd (389, -0.6), // AP
    make_dd (401, 0.6), // AP
    // er
    make_dd (1367, 0.35, 71 + 70), // AP + chorus
    make_dd (0), // damp filter
    make_dd (1787, 0., 261), // Delay + chorus
    make_dd (33), // delay to allow block processing
    // loop1
    make_dd (977, 0.5 /*overridden*/, 51), // AP with variable gain + chorus
    make_dd (2819), // Delay
    make_dd (0), // damp filter
    make_dd (863, -0.5 /*overridden*/), // AP with variable gain
    make_dd (1021), // Delay
    make_dd (1453, 0.618), // Delay + chorus
    make_dd (787), // delay (allows block processing)
    // loop2
    make_dd (947, 0.5 /*overridden*/, 67), // AP with variable gain + chorus
    make_dd (3191), // Delay
    make_dd (0), // damp filter
    make_dd (887, -0.5 /*overridden*/), // AP with variable gain
    make_dd (1049), // Delay
    make_dd (1367, 0.618), // Delay + chorus
    make_dd (647)); // delay (allows block processing)
}

struct rev1_spec {
  static constexpr auto values {get_rev1_delays()};
};
//------------------------------------------------------------------------------
template <class Spec>
class reverb_tool {
public:
  //----------------------------------------------------------------------------
  static constexpr uint get_required_size()
  {
    uint ret = 0;
    for (uint i = 0; i < decltype (_stage) {}.size(); ++i) {
      ret += get_buffer_size (i);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<s16> mem)
  {
    for (uint i = 0; i < _stage.size(); ++i) {
      _stage[i].z = mem.cut_head (get_buffer_size (i)).data();
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (xspan<T> io, q0_15 g = q0_15 {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) == 1);

    for (uint i = 0; i < io.size(); ++i) {
      auto y1 = T::from (_stage[Idx].z[0]);
      auto gv = (g.value() == 0) ? dd.g : g;
      auto v  = (T) ((gv.max() - gv) * io[i]);
      v       = (T) (v + y1 * gv);
      io[i]   = v;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (xspan<T> io, q0_15 g = q0_15 {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) == 1);

    for (uint i = 0; i < io.size(); ++i) {
      auto y1 = T::from (_stage[Idx].z[0]);
      auto gv = (g.value() == 0) ? dd.g : g;
      auto v  = (T) ((gv.max() - gv) * io[i]);
      v       = (T) (v + y1 * gv);
      io[i]   = (T) (io[i] - v);
    }
  }
  //----------------------------------------------------------------------------
  // Pure delays with "size > blocksize" are placed before feedback loops to
  // enable block processing, so it is possible to fetch the future
  // feedbacks/outputs (minus the first) at once and then to push them all at
  // once. This is exactly what fetch/push accomplish.
  template <uint Idx, class T>
  void fetch (xspan<T> dst)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (sz >= dst.size());

    for (uint i = 0; i < dst.size(); ++i) {
      dst[i] = get<Idx, T> (dd.spls.to_int() - i);
    }
  }
  // see comment on fetch
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (sz >= src.size());

    for (uint i = 0; i < src.size(); ++i) {
      _stage[Idx].z[_stage[Idx].pos] = src[i].value();
      advance_pos<Idx>();
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run (xspan<T> io, xspan<q0_15> g = xspan<q0_15> {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    static_assert (dd.mod.value() == 0);

    for (uint i = 0; i < io.size(); ++i) {
      // no interpolation
      T push {io[i]};
      T out = run_plain_delay<Idx, T>();

      if constexpr (dd.g.value() != 0) {
        auto [outv, pushv] = run_allpass<Idx, T> (push, out, g ? g[i] : dd.g);
        out                = outv;
        push               = pushv;
      }
      io[i]                          = out;
      _stage[Idx].z[_stage[Idx].pos] = push.value();
      advance_pos<Idx>();
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  auto run_mod (xspan<T> io, xspan<q0_15> lfo, xspan<q0_15> g = xspan<q0_15> {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    static_assert (dd.mod.value() != 0);

    for (uint i = 0; i < io.size(); ++i) {
      // interpolation
      T push {io[i]};
      T out = run_thiran<Idx, T> (lfo[i]);

      if constexpr (dd.g.value() != 0) {
        auto [outv, pushv] = run_allpass<Idx, T> (push, out, g ? g[i] : dd.g);
        out                = outv;
        push               = pushv;
      }
      io[i]                          = out;
      _stage[Idx].z[_stage[Idx].pos] = push.value();
      advance_pos<Idx>();
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <uint Idx>
  void advance_pos()
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    ++_stage[Idx].pos;
    if (_stage[Idx].pos == sz) {
      _stage[Idx].pos = 0;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T run_plain_delay()
  {
    constexpr delay_data dd = Spec::values[Idx];
    return get<Idx, T> (dd.spls.to_int());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T run_thiran (q0_15 lfo)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    auto fixpt_spls  = dd.spls.add_sign() + (lfo * dd.mod.add_sign());
    auto n_spls      = fixpt_spls.to_int();
    auto n_spls_frac = fixpt_spls.fractional().to_lossless();
    auto d           = n_spls_frac + num {0.418f};
    auto one         = fixpt_p<1, 1, 0>::from_int (1);
    auto co_thiran   = (one - d) / (one + d); // 0.4104 to -1
    auto a           = co_thiran.cast<fixpt_p<1, 0, 23>>();

    auto z0 = get<Idx, T> (n_spls - 1);
    auto z1 = get<Idx, T> (n_spls);
    auto y1 = T::from (_stage[Idx].z[y1_pos]);

    auto ret              = (T) (z0 * a + z1 - a * y1);
    auto ret_dbg          = ret.to_float();
    _stage[Idx].z[y1_pos] = ret.value();
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T get (uint delay_spls)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    assert (delay_spls < sz);

    uint z = _stage[Idx].pos - delay_spls;
    z += (z >= sz) ? sz : 0;
    return T::from (_stage[Idx].z[z]);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  auto run_allpass (T in, T yn, q0_15 g)
  {
    auto u = in + yn * g;
    auto x = yn - u * g;
    return std::make_tuple (x.cast (T {}), u.cast (T {}));
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_size (uint i)
  {
    uint ret = Spec::values[i].spls.to_int() + 1; // spls + 1 = size
    ret += Spec::values[i].mod.to_int(); // add mod spls to the size
    return ret;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_buffer_size (uint i)
  {
    uint ret = get_delay_size (i);
    ret += (uint) Spec::values[i].mod.value() != 0; // thiran state
    return ret;
  }
  //----------------------------------------------------------------------------
  struct stage {
    s16* z {};
    uint pos {};
  };
  std::array<stage, Spec::values.size()> _stage;
};
}} // namespace detail::lofiverb
//------------------------------------------------------------------------------
class lofiverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct mode_tag {};
  void set (mode_tag, int v)
  {
    if (v == _param.mode) {
      return;
    }
    switch (v) {
    case mode::rev1: {
      _param.mode_headroom_gain = 2.f + 0.8f + 0.855f + 0.443f;
      auto& rev                 = _modes.emplace<rev1_type>();
      rev.reset_memory (_mem_reverb);
    } break;
    case mode::dummy: {
    } break;
    default:
      return;
    }
    _param.mode          = v;
    _n_processed_samples = 0; // trigger the control block on first sample
    xspan_memset (_mem_reverb, 0);
    xspan_memset (xspan {_feedback}, 0);
  }
  struct mode {
    enum { rev1, dummy };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (0, make_cstr_array ("Reverb 1", "Dummy!"), 48);
  }
  //----------------------------------------------------------------------------
  struct decay_tag {};
  void set (decay_tag, float v) { _param_smooth.target().decay = v * 0.01f; }

  static constexpr auto get_parameter (decay_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct mod_tag {};
  void set (mod_tag, float v)
  {
    v *= 0.01f;
    if (v == _param_smooth.target().mod) {
      return;
    }
    _param_smooth.target().mod = v;
    update_mod();
  }

  static constexpr auto get_parameter (mod_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct character_tag {};
  void set (character_tag, float v)
  {
    _param_smooth.target().character = v * 0.01f;
  }

  static constexpr auto get_parameter (character_tag)
  {
    return float_param ("%", 0., 100., 50., 0.001);
  }
  //----------------------------------------------------------------------------
  struct damp_tag {};
  void set (damp_tag, float v) { _param.damp = v * 0.01f; }

  static constexpr auto get_parameter (damp_tag)
  {
    return float_param ("%", 0., 100., 30., 0.001);
  }
  //----------------------------------------------------------------------------
  struct tilt_tag {};
  void set (tilt_tag, float v)
  {
    if (v == _param.tilt_db) {
      return;
    }
    _param.tilt_db = v;
    _tilt.reset_coeffs (
      vec_set<2> (330.f), vec_set<2> (0.5f), vec_set<2> ((float) -v), t_spl);
  }

  static constexpr auto get_parameter (tilt_tag)
  {
    return float_param ("dB", -14., 14., 0., 0.1);
  }
  //----------------------------------------------------------------------------
  static constexpr uint max_predelay_qb = 4;
  struct predelay_tag {};
  void set (predelay_tag, float v) { _param.predelay = v; }

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("quarters", 0., max_predelay_qb, 0., 0.001);
  }
  //----------------------------------------------------------------------------
  struct er_tag {};
  void set (er_tag, float v)
  {
    v *= 0.01f;
    v *= v;
    _param_smooth.target().er = v;
  }

  static constexpr auto get_parameter (er_tag)
  {
    return float_param ("%", 0., 100., 25., 0.001);
  }
  //----------------------------------------------------------------------------
  struct stereo_tag {};
  void set (stereo_tag, float v) { _param_smooth.target().stereo = v * 0.01f; }

  static constexpr auto get_parameter (stereo_tag)
  {
    return float_param ("%", -100., 100., 100., 0.01);
  }
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    auto beat_hz         = pc.get_play_state().bpm * (1.f / 60.f);
    _1_4beat_spls        = (0.25f / beat_hz) * srate;
    _n_processed_samples = 0;

    _resampler.reset (
      srate,
      pc.get_sample_rate(),
      15000,
      15000,
      32,
      16,
      210,
      true,
      max_block_size,
      6 * 1024);

    _tilt.reset_states();

    // resize memory
    uint predelay_spls = std::ceil (_1_4beat_spls * max_predelay_qb) * 2;
    uint rev_spls      = 0;
    mp11::mp_for_each<decltype (_modes)> ([&] (auto mode) {
      rev_spls = std::max (mode.get_required_size(), rev_spls);
    });
    _mem.clear();
    _mem.resize (predelay_spls + rev_spls);
    _mem_reverb = xspan {_mem};
    _predelay.reset (_mem_reverb.cut_head (predelay_spls), 2);

    // set defaults
    mp11::mp_for_each<parameters> ([&] (auto param) {
      set (param, get_parameter (param).min);
      if constexpr (!is_choice<decltype (get_parameter (param))>) {
        set (param, get_parameter (param).max); // max might be not yet impl.
      }
      else {
        set (param, get_parameter (param).min + 1);
      }
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
#if MOD_DBG_DENORMALS
    feenableexcept (FE_INVALID);
#endif
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    _resampler.process (outs, ins, samples, [=] (auto io) {
      process_block (io);
    });
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    mode_tag,
    character_tag,
    damp_tag,
    decay_tag,
    predelay_tag,
    tilt_tag,
    er_tag,
    mod_tag,
    stereo_tag>;
  //----------------------------------------------------------------------------
private:
  using q0_15 = detail::lofiverb::q0_15;
  using q1_14 = detail::lofiverb::q1_14;
  //----------------------------------------------------------------------------
  static constexpr uint  max_block_size = 32;
  static constexpr uint  n_channels     = 2;
  static constexpr uint  srate          = 33600;
  static constexpr float t_spl          = (float) (1. / srate);
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  void process_block (xspan<std::array<float, 2>> io)
  {
    // TODO: test __fp16 (half-precision float storage).
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    using fixptype = q0_15;
    array2d<fixptype, 2, max_block_size> fixedp;
    loop_parameters                      pars;

    // tilt! + param smoothing + int conversion
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = vec_to_array (_tilt.tick (vec_from_array (io[i])));
      _param_smooth.tick();
      pars.stereo[i] = _param_smooth.get().stereo;
      pars.er[i].load_float (_param_smooth.get().er);
      pars.decay[i].load_float (_param_smooth.get().decay);
      pars.character[i].load_float (_param_smooth.get().character);
      pars.mod[i].load_float (_param_smooth.get().mod);
    }

    // predelay + int conversion
    if (_param.predelay != 0) {
      uint  predelay_spls = _1_4beat_spls * _param.predelay;
      float gain          = 1.f / _param.mode_headroom_gain;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        fixedp[i][0].load (_predelay.get (predelay_spls, 0));
        fixedp[i][1].load (_predelay.get (predelay_spls, 1));
        auto l = std::clamp (io[i][0], -0.98f, 0.98f);
        auto r = std::clamp (io[i][1], -0.98f, 0.98f);
        // the maximum gain the loop has is computed theoretically and applied
        // before input, so fixed point scaling is not needed: no need for
        // integer bits.
        l *= gain;
        r *= gain;
        auto push = make_array (
          fixptype::from_float (l).value(), fixptype::from_float (r).value());
        _predelay.push (xspan {push});
      }
    }
    else {
      float gain = 1.f / _param.mode_headroom_gain;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        auto l = std::clamp (io[i][0], -0.98f, 0.98f);
        auto r = std::clamp (io[i][1], -0.98f, 0.98f);
        l *= gain;
        r *= gain;
        fixedp[i][0].load_float (l);
        fixedp[i][1].load_float (r);
      }
    }

    // main loop
    process_rev1 (xspan {fixedp.data(), io.size()}, pars);

    // float conversion
    auto gain = _param.mode_headroom_gain;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = fixedp[i][0].to_float();
      auto r = fixedp[i][1].to_float();
      l      = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r * gain;
        io[i][1] = l * gain;
      }
      else {
        io[i][0] = l * gain;
        io[i][1] = r * gain;
      }
    }
  }
  //----------------------------------------------------------------------------
  struct loop_parameters;
  //----------------------------------------------------------------------------
  void process_rev1 (xspan<std::array<q0_15, 2>> io, loop_parameters& par)
  {
    auto& rev = std::get<rev1_type> (_modes);

    using fixptarr    = std::array<q0_15, max_block_size>;
    using fixptarr_fb = std::array<q0_15, max_block_size + 1>;

    fixptarr late_in_arr;
    fixptarr lfo1;
    fixptarr lfo2;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = ((io[i][0] + io[i][1]) >> 1).resize<-1>();
      // ER lfo
      auto lfo = _lfo_er.tick_sine_fixpt().spec_cast<q0_15>().value();
      lfo1[i].load (lfo[0]);
      lfo2[i].load (lfo[1]);
      auto er_amt       = (q1_14) (num {0.5} + (par.er[i] >> 2));
      auto lfo2final    = lfo2[i] * er_amt;
      auto lfo2finalflt = lfo2final.to_float();

      lfo2[i] = (q0_15) (lfo2final);
      // decay fixup
      auto decay   = (q1_14) (num {1} - par.decay[i]);
      decay        = (q1_14) (num {1} - decay * decay);
      par.decay[i] = (q1_14) (num {0.6f} + decay * num {0.39f});
    }

    // diffusion -----------------------------
    rev.run<0> (late_in);
    rev.run<1> (late_in);
    rev.run<2> (late_in);
    rev.run<3> (late_in);

    // ER -----------------------------
    fixptarr    early1_arr;
    fixptarr    early1b_arr;
    fixptarr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    er2[0] = _feedback[0];
    rev.fetch<7> (er2.advanced (1));
    _feedback[0] = er2.back();

    for (uint i = 0; i < io.size(); ++i) {
      // apply ER feedback
      er1[i] = (q0_15) (er2[i] * num {0.2});
      er1[i] = (q0_15) ((late_in[i] + er1[i]) * par.decay[i]);
    }
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run_mod<4> (er1, lfo1);
    rev.run_lp<5> (er1, q0_15::from_float (0.0001f + 0.17f * _param.damp));
    xspan_memcpy (er1b, er1);
    rev.run_mod<6> (er1b, lfo2);
    rev.push<7> (er1b.to_const()); // feedback point

    // Late -----------------------------
    fixptarr    late_arr;
    fixptarr_fb l_arr;
    fixptarr_fb r_arr;
    fixptarr    g_arr;

    auto late = xspan {late_arr.data(), io.size()};
    auto g    = xspan {g_arr.data(), io.size()};
    auto l    = xspan {l_arr.data(), io.size() + 1}; // +1: Feedback on head
    auto r    = xspan {r_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    l[0] = _feedback[1];
    r[0] = _feedback[2];
    rev.fetch<14> (l.advanced (1));
    rev.fetch<21> (r.advanced (1));
    _feedback[1] = l.back();
    _feedback[2] = r.back();

    for (uint i = 0; i < io.size(); ++i) {
      // clang-format off
      auto mod = (q1_14)
        (num {0.25} + (num {1} - par.mod[i]) * num {0.75});
      // clang-format on
      auto lfo = _lfo.tick_sine_fixpt().spec_cast<q0_15>().value();
      lfo1[i].load (lfo[0]);
      lfo1[i] = (q0_15) (lfo1[i] * mod);
      lfo2[i].load (lfo[1]);
      lfo2[i] = (q0_15) (lfo1[i] * mod);
      // prepare input with feedback
      late[i] = (q0_15) (late_in[i] + r[i] * par.decay[i]);
      // add ER blend to late in
      auto fact  = (q0_15) (par.er[i] * num {0.4});
      auto er_in = (q0_15) ((er1[i] + er2[i]) * fact);
      late[i]    = (q0_15) (late[i] + er_in);

      // g character
      // clang-format off
      g[i] = (q0_15)
        (num {0.618} + par.character[i] * num {(0.707 - 0.618) * 2.});
      // clang-format on
    }
    r.cut_head (1); // drop feedback sample from previous block

    auto late_damp_1 = q1_14::from_float (_param.damp);
    late_damp_1      = (q1_14) (num {0.9} - (num {0.9} * late_damp_1));
    late_damp_1      = (q1_14) (num {1} - late_damp_1 * late_damp_1);
    late_damp_1      = (q1_14) (late_damp_1 * num {0.4});
    auto late_damp   = (q0_15) late_damp_1;

    // hp skipped for now.
    rev.run_mod<8> (late, lfo1, g);
    rev.run<9> (late);
    rev.run_lp<10> (late, late_damp);
    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run<11> (late, g);
    rev.run<12> (late);
    rev.run<13> (late);
    rev.push<14> (late.to_const()); // feedback point

    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      late[i] = (q0_15) (late_in[i] + l[i] * par.decay[i]);

      auto fact  = (q0_15) (par.er[i] * num {0.4});
      auto er_in = (q0_15) ((er1[i] - er2[i]) * fact);
      late[i]    = (q0_15) (late[i] + er_in);
    }
    l.cut_head (1); // drop feedback sample from previous block

    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run_mod<15> (late, lfo2, g);
    rev.run<16> (late);
    rev.run_lp<17> (late, late_damp);
    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run<18> (late, g);
    rev.run<19> (late);
    rev.run<20> (late);
    rev.push<21> (late.to_const()); // feedback point

    // Mixdown
    for (uint i = 0; i < io.size(); ++i) {
      io[i][0] = (q0_15) (-er1[i] * num {0.825f});
      io[i][0] = (q0_15) (io[i][0] - er2[i] * num {0.423});
      io[i][0] = (q0_15) (io[i][0] * par.er[i]);
      io[i][0] = (q0_15) (io[i][0] + l[i]);

      io[i][1] = (q0_15) (-er1[i] * num {0.855f});
      io[i][1] = (q0_15) (io[i][1] + er2[i] * num {0.443});
      io[i][1] = (q0_15) (io[i][1] * par.er[i]);
      io[i][1] = (q0_15) (io[i][1] + r[i]);
    }
  }
  //----------------------------------------------------------------------------
  void update_mod()
  {
    auto mod = _param_smooth.target().mod;
    // TODO: Two Lfos are not required just to change sign (180 separation).
    _lfo.set_freq (vec_set<2> (0.1f + mod * 1.2f), t_spl);
    _lfo.set_phase (_lfo.get_phase().set_at_uniform_spacing());
    _lfo_er.set_freq (vec_set<2> (0.3f + mod * 0.3f), t_spl);
    _lfo_er.set_phase (_lfo_er.get_phase().set_at_uniform_spacing());
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   mode;
    float mode_headroom_gain;
    float tilt_db;
    float predelay;
    float damp;
  };
  //----------------------------------------------------------------------------
  struct smoothed_parameters {
    float mod;
    float stereo;
    float er;
    float decay;
    float character;
    // dry, wet, ducker/gate
  };
  //----------------------------------------------------------------------------
  struct loop_parameters {
    std::array<q1_14, max_block_size> er;
    std::array<q1_14, max_block_size> decay;
    std::array<q1_14, max_block_size> character;
    std::array<q1_14, max_block_size> mod;
    std::array<float, max_block_size> stereo;
  };
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>            _resampler {};
  part_class_array<tilt_eq, f32_x2>    _tilt {};
  static_delay_line<s16, false, false> _predelay;

  lfo<2> _lfo;
  lfo<2> _lfo_er;

  std::array<q0_15, 8> _feedback;

  using rev1_type = detail::lofiverb::reverb_tool<detail::lofiverb::rev1_spec>;
  std::variant<rev1_type> _modes;

  uint  _n_processed_samples;
  float _1_4beat_spls;

  std::vector<s16, overaligned_allocator<s16, 16>> _mem;
  xspan<s16>                                       _mem_reverb;
};
//------------------------------------------------------------------------------
} // namespace artv
