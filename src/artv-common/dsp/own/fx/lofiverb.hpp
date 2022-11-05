#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "artv-common/dsp/own/classes/block_resampler.hpp"
#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/ducker.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/value_smoother.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/float.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {
namespace detail { namespace lofiverb {
//------------------------------------------------------------------------------
using q0_15          = fixpt_dt<1, 0, 15>;
using q0_15r         = fixpt_dr<1, 0, 15>;
using q1_14          = fixpt_dt<1, 1, 14>;
using fixpt_spls     = fixpt_dt<0, 14, 0>;
using fixpt_spls_mod = fixpt_dt<0, 9, 0>;
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
    make_dd (0, 0.98), // fixed highpass filter
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
  // Pure delays with "size > blocksize" are placed before feedback loops to
  // enable block processing, so it is possible to fetch the future
  // feedbacks/outputs (minus the first) at once and then to push them all at
  // once. This is exactly what fetch/push accomplish.
  //
  // dst has to contain one element more at the head, where the previous output
  // will be placed. Once the feedback is applied to a current input, this head
  // feedback sample (last output) can be dropped.
  template <uint Idx, class T>
  void fetch_block_plus_one (xspan<T> dst)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (dst && dst.size());
    assert (sz >= dst.size());

    constexpr auto out_prev = 1;

    uint block1 = _stage[Idx].pos - dd.spls.to_int() - out_prev;
    block1 += (block1 >= sz) ? sz : 0;
    uint end = _stage[Idx].pos - dd.spls.to_int() + dst.size() - out_prev - 1;
    end += (end >= sz) ? sz : 0;

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = dst.size();
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = dst.size() - block1sz;
    }
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      memcpy (dst.data(), &_stage[Idx].z[block1], block1sz * sizeof dst[0]);
      memcpy (
        dst.data() + block1sz, &_stage[Idx].z[0], block2sz * sizeof dst[0]);
    }
    else {
      for (uint i = 0; i < block1sz; ++i) {
        dst[i] = float16::decode (_stage[Idx].z[block1 + i]);
      }
      for (uint i = 0; i < block2sz; ++i) {
        dst[i + block1sz] = float16::decode (_stage[Idx].z[i]);
      }
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
    assert (src && src.size());
    assert (sz >= src.size());

    uint block1 = _stage[Idx].pos;
    uint end    = block1 + src.size() - 1;
    end -= (end >= sz) ? sz : 0;
    _stage[Idx].pos = end;
    advance_pos<Idx>();

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = src.size();
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = src.size() - block1sz;
    }
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      memcpy (&_stage[Idx].z[block1], src.data(), block1sz * sizeof src[0]);
      memcpy (
        &_stage[Idx].z[0], src.data() + block1sz, block2sz * sizeof src[0]);
    }
    else {
      for (uint i = 0; i < block1sz; ++i) {
        _stage[Idx].z[block1 + i] = float16::encode (src[i]);
      }
      for (uint i = 0; i < block2sz; ++i) {
        _stage[Idx].z[i] = float16::encode (src[i + block1sz]);
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  void run_lp (xspan<T> io, q0_15 g = q0_15 {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);

    auto y1 = T::from (_stage[Idx].z[0]);

    for (uint i = 0; i < io.size(); ++i) {
      auto gv = (g.value() == 0) ? dd.g : g;
      auto v  = (T) ((gv.max() - gv) * io[i]);
      v       = (T) (v + y1 * gv);
      y1      = v;
      io[i]   = v;
    }
    _stage[Idx].z[0] = y1.value();
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_lp (xspan<float> io, float g = 0.f)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);

    auto y1 = float16::decode (_stage[Idx].z[0]);
    for (uint i = 0; i < io.size(); ++i) {
      auto gv = (g == 0.f) ? dd.g.to_float() : g;
      auto v  = (1.f - gv) * io[i];
      v       = v + y1 * gv;
      y1      = v;
      io[i]   = v;
    }
    _stage[Idx].z[0] = float16::encode (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  void run_hp (xspan<T> io, q0_15 g = q0_15 {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);

    auto y1 = T::from (_stage[Idx].z[0]);

    for (uint i = 0; i < io.size(); ++i) {
      auto gv = (g.value() == 0) ? dd.g : g;
      auto v  = (T) ((gv.max() - gv) * io[i]);
      v       = (T) (v + y1 * gv);
      y1      = v;
      io[i]   = (T) (io[i] - v);
    }
    _stage[Idx].z[0] = y1.value();
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_hp (xspan<float> io, float g = 0.f)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);

    auto y1 = float16::decode (_stage[Idx].z[0]);
    for (uint i = 0; i < io.size(); ++i) {
      auto gv = (g == 0.f) ? dd.g.to_float() : g;
      auto v  = (1.f - gv) * io[i];
      v       = v + y1 * gv;
      y1      = v;
      io[i]   = io[i] - v;
    }
    _stage[Idx].z[0] = float16::encode (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
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
  template <uint Idx>
  void run (xspan<float> io, xspan<float> g = xspan<float> {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    static_assert (dd.mod.value() == 0);

    for (uint i = 0; i < io.size(); ++i) {
      // no interpolation
      float push {io[i]};
      float out = run_plain_delay<Idx, float>();

      if constexpr (dd.g.value() != 0) {
        auto [outv, pushv]
          = run_allpass<Idx> (push, out, g ? g[i] : dd.g.to_float());
        out  = outv;
        push = pushv;
      }
      io[i]                          = out;
      _stage[Idx].z[_stage[Idx].pos] = float16::encode (push);
      advance_pos<Idx>();
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
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
  template <uint Idx>
  auto run_mod (
    xspan<float> io,
    xspan<float> lfo,
    xspan<float> g = xspan<float> {})
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);

    static_assert (dd.mod.value() != 0);

    for (uint i = 0; i < io.size(); ++i) {
      // interpolation
      float push {io[i]};
      float out = run_thiran<Idx> (lfo[i]);

      if constexpr (dd.g.value() != 0) {
        auto [outv, pushv]
          = run_allpass<Idx> (push, out, g ? g[i] : dd.g.to_float());
        out  = outv;
        push = pushv;
      }
      io[i]                          = out;
      _stage[Idx].z[_stage[Idx].pos] = float16::encode (push);
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
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  T run_thiran (q0_15 lfo)
  {
    // y1 doesn't need to be converted once per sample, can be done per
    // block (at the expense of one extra buffer). Small blocks would then be
    // more noisy!
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    auto fixpt_spls  = dd.spls.add_sign() + (lfo * dd.mod.add_sign());
    auto n_spls      = fixpt_spls.to_int();
    auto n_spls_frac = fixpt_spls.fractional().to_dynamic();
    auto d           = n_spls_frac + num {0.418f};
    auto one         = fixpt_dt<1, 1, 0>::from_int (1);
    auto co_thiran   = (one - d) / (one + d); // 0.4104 to -1
    auto a           = co_thiran.cast<fixpt<1, 0, 23>>();

    auto z0 = get<Idx, T> (n_spls - 1);
    auto z1 = get<Idx, T> (n_spls);
    auto y1 = T::from (_stage[Idx].z[y1_pos]);

    auto ret              = (T) (z0 * a + z1 - a * y1);
    _stage[Idx].z[y1_pos] = ret.value();
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  float run_thiran (float lfo)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    // y1 doesn't need to be converted once per sample, can be done per
    // block (at the expense of one extra buffer). Small blocks would then be
    // more noisy!
    float fixpt_spls  = dd.spls.to_float() + (lfo * dd.mod.to_float());
    uint  n_spls      = (uint) fixpt_spls;
    float n_spls_frac = fixpt_spls - n_spls;
    float d           = n_spls_frac + 0.418f;
    float a           = (1.f - d) / (1.f + d); // 0.4104 to -1

    float z0 = get<Idx, float> (n_spls - 1);
    float z1 = get<Idx, float> (n_spls);
    float y1 = float16::decode (_stage[Idx].z[y1_pos]);

    float ret             = z0 * a + z1 - a * y1;
    _stage[Idx].z[y1_pos] = float16::encode (ret);
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
    if constexpr (std::is_floating_point_v<T>) {
      return float16::decode (_stage[Idx].z[z]);
    }
    else {
      return T::from (_stage[Idx].z[z]);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  auto run_allpass (T in, T yn, q0_15 g)
  {
    auto u = in + yn * g;
    auto x = yn - u * g;
    return std::make_tuple (x.cast (T {}), u.cast (T {}));
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  auto run_allpass (float in, float yn, float g)
  {
    auto u = in + yn * g;
    auto x = yn - u * g;
    return std::make_tuple (x, u);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_size (uint i)
  {
    uint ret = Spec::values[i].spls.to_int() + 1; // spls + 1 = size
    // pure delays have 1 extra sample to be able to return the previous output
    // (delaying one cycle the time it is overwritten).
    auto mod_val = Spec::values[i].mod.value();
    auto g_val   = Spec::values[i].g.value();
    ret += (uint) (mod_val == 0 && g_val == 0);
    ret += Spec::values[i].mod.to_int(); // add mod spls to the size
    return ret;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_buffer_size (uint i)
  {
    uint ret = get_delay_size (i);
    ret += (uint) (Spec::values[i].mod.value() != 0); // thiran state
    return ret;
  }
  //----------------------------------------------------------------------------
  using float16 = f16pack<5, -1, f16pack_dftz>;
  struct stage {
    s16* z {};
    uint pos {};
  };
  std::array<stage, Spec::values.size()> _stage;
};
}} // namespace detail::lofiverb
//------------------------------------------------------------------------------
// A reverb using 16-bit fixed-point arithmetic on the main loop. One design
// criteria has been for it to be extremely CPU friendly.
// Notice: __fp16 (half-precision floating point storage) could achieve the
// same memory savings with more dynamic range. This is still kept as a 16-bit
// fixed-point reverb.
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
    case mode::rev1_flt:
    case mode::rev1: {
      _param.mode_headroom_gain = 2.f + 0.8f + 0.855f + 0.443f;
      auto& rev                 = _modes.emplace<rev1_type>();
      rev.reset_memory (_mem_reverb);
    } break;
    default:
      return;
    }
    _param.mode          = v;
    _n_processed_samples = 0; // trigger the control block on first sample
    xspan_memset (_mem_reverb, 0);
  }
  struct mode {
    enum { rev1_flt, rev1 };
  };

  static constexpr auto get_parameter (mode_tag)
  {
    return choice_param (0, make_cstr_array ("Reverb1", "Reverb1 16-bit"), 48);
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
    v *= 0.01f;
    if (v == _param.tilt) {
      return;
    }
    _param.tilt = v;
    _tilt.reset_coeffs (
      vec_set<2> (330.f),
      vec_set<2> (0.5f),
      vec_set<2> ((float) v * -14.f),
      t_spl);
  }

  static constexpr auto get_parameter (tilt_tag)
  {
    return float_param ("%", -100, 100., 0., 0.1);
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
  struct ducking_threshold_tag {};
  void set (ducking_threshold_tag, float v)
  {
    if (v == _param.ducking_threshold) {
      return;
    }
    _param.ducking_threshold = v;
    _ducker.set_threshold (vec_set<2> (v));
  }

  static constexpr auto get_parameter (ducking_threshold_tag)
  {
    return float_param ("dB", -40.f, 12.f, 12.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct ducking_speed_tag {};
  void set (ducking_speed_tag, float v)
  {
    if (v == _param.ducking_speed) {
      return;
    }
    _param.ducking_speed = v;
    v *= 0.01f;
    _ducker.set_speed (vec_set<2> (v * v), t_spl);
  }

  static constexpr auto get_parameter (ducking_speed_tag)
  {
    return float_param ("%", 0.f, 100.f, 10.f, 0.01f);
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
    _ducker.reset();

    // resize memory
    uint predelay_spls = std::ceil (_1_4beat_spls * max_predelay_qb) * 2;
    uint rev_spls      = 0;
    mp11::mp_for_each<decltype (_modes)> ([&] (auto mode) {
      rev_spls = std::max (mode.get_required_size(), rev_spls);
    });
    uint predelay_flt = predelay_spls * (sizeof (float) / sizeof _mem[0]);

    _mem.clear();
    _mem.resize (predelay_spls + predelay_flt + rev_spls);

    _mem_reverb = xspan {_mem};
    _predelay.reset (_mem_reverb.cut_head (predelay_spls), 2);
    _predelay_flt.reset (_mem_reverb.cut_head (predelay_flt).cast<float>(), 2);

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
      if ((_param.mode % 2) == 0) {
        process_float (io);
      }
      else {
        process_fixed_point (io);
      }
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
    stereo_tag,
    ducking_threshold_tag,
    ducking_speed_tag>;
  //----------------------------------------------------------------------------
private:
  // q0_15 truncates when dropping fractional bits (leaky).
  // q0_15r rounds to the nearest (never reaches full zero).
  //
  // Using q0_15 as default, with q0_15r at some points compensate the
  // truncating leakage.
  using q0_15  = detail::lofiverb::q0_15;
  using q0_15r = detail::lofiverb::q0_15r;
  using q1_14  = detail::lofiverb::q1_14;
  //----------------------------------------------------------------------------
  static constexpr uint  max_block_size = 32;
  static constexpr uint  n_channels     = 2;
  static constexpr uint  srate          = 33600;
  static constexpr float t_spl          = (float) (1. / srate);
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters;
  struct smoothed_parameters;
  //----------------------------------------------------------------------------
  void process_fixed_point (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    using fixptype = q0_15;
    array2d<fixptype, 2, max_block_size> wet;
    std::array<f32_x2, max_block_size>   ducker_gain;
    loop_parameters                      pars;

    // tilt! + param smoothing + int conversion
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wet     = {io[i][0], io[i][1]};
      ducker_gain[i] = _ducker.tick (wet);
      io[i]          = vec_to_array (_tilt.tick (vec_from_array (io[i])));
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
        wet[i][0].load (_predelay.get (predelay_spls, 0));
        wet[i][1].load (_predelay.get (predelay_spls, 1));
        auto l = std::clamp (io[i][0], -0.98f, 0.98f) * gain;
        auto r = std::clamp (io[i][1], -0.98f, 0.98f) * gain;
        // the maximum gain the loop has is computed theoretically and applied
        // before input, so fixed point scaling is not needed: no need for
        // integer bits.
#if 0 // will process afterwards
        auto            wn    = _whitenoise();
        constexpr float scale = 1.f / 0xffff;
        l += wn[0] + wn[1]; // TPDF
        r += wn[0] - wn[1];
#endif
        auto push = make_array (
          fixptype::from_float (l).value(), fixptype::from_float (r).value());
        _predelay.push (xspan {push});
      }
    }
    else {
      float gain = 1.f / _param.mode_headroom_gain;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        auto l = std::clamp (io[i][0], -0.98f, 0.98f) * gain;
        auto r = std::clamp (io[i][1], -0.98f, 0.98f) * gain;
        wet[i][0].load_float (l);
        wet[i][1].load_float (r);
      }
    }

    // main loop
    process_rev1 (xspan {wet.data(), io.size()}, pars);

    // float conversion
    auto gain = _param.mode_headroom_gain;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0].to_float() * gain * ducker_gain[i][0];
      auto r = wet[i][1].to_float() * gain * ducker_gain[i][1];
      l      = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
  //----------------------------------------------------------------------------
  void process_float (xspan<std::array<float, 2>> io)
  {
    assert (io.size() <= max_block_size);

    // clip + convert to u16
    array2d<float, 2, max_block_size>  wet;
    std::array<f32_x2, max_block_size> ducker_gain;
    loop_parameters_flt                pars;

    // tilt! + param smoothing + int conversion
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wet     = {io[i][0], io[i][1]};
      ducker_gain[i] = _ducker.tick (wet);
      io[i]          = vec_to_array (_tilt.tick (vec_from_array (io[i])));
      _param_smooth.tick();
      pars.stereo[i]    = _param_smooth.get().stereo;
      pars.er[i]        = _param_smooth.get().er;
      pars.decay[i]     = _param_smooth.get().decay;
      pars.character[i] = _param_smooth.get().character;
      pars.mod[i]       = _param_smooth.get().mod;
    }

    // predelay + int conversion
    if (_param.predelay != 0) {
      uint  predelay_spls = _1_4beat_spls * _param.predelay;
      float gain          = 1.f / _param.mode_headroom_gain;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0] = _predelay_flt.get (predelay_spls, 0);
        wet[i][1] = _predelay_flt.get (predelay_spls, 1);
        auto l    = std::clamp (io[i][0], -0.999f, 0.999f) * gain;
        auto r    = std::clamp (io[i][1], -0.999f, 0.999f) * gain;
        auto push = make_array (l, r);
        _predelay_flt.push (xspan {push});
      }
    }
    else {
      float gain = 1.f / _param.mode_headroom_gain;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        auto l    = std::clamp (io[i][0], -0.999f, 0.999f) * gain;
        auto r    = std::clamp (io[i][1], -0.999f, 0.999f) * gain;
        wet[i][0] = l;
        wet[i][1] = r;
      }
    }

    // main loop
    process_rev1_flt (xspan {wet.data(), io.size()}, pars);

    auto gain = _param.mode_headroom_gain;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0] * ducker_gain[i][0] * gain;
      auto r = wet[i][1] * ducker_gain[i][1] * gain;
      l      = r * (1 - abs (pars.stereo[i])) + l * abs (pars.stereo[i]);
      if (pars.stereo[i] < 0) {
        io[i][0] = r;
        io[i][1] = l;
      }
      else {
        io[i][0] = l;
        io[i][1] = r;
      }
    }
  }
  //----------------------------------------------------------------------------
  struct loop_parameters;
  //----------------------------------------------------------------------------
  void process_rev1 (xspan<std::array<q0_15, 2>> io, loop_parameters& par)
  {
    auto& rev = std::get<rev1_type> (_modes);

    using arr    = std::array<q0_15, max_block_size>;
    using arr_fb = std::array<q0_15, max_block_size + 1>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = ((io[i][0] + io[i][1]) >> 1).resize<-1>();
      // ER lfo
      auto lfo = _lfo_er.tick_sine_fixpt().spec_cast<q0_15>().value();
      lfo1[i].load (lfo[0]);
      lfo2[i].load (lfo[1]);
      auto er_amt    = (q1_14) (num {0.5} + (par.er[i] >> 1));
      auto lfo2final = lfo2[i] * er_amt;
      lfo2[i]        = (q0_15) (lfo2final);
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
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    rev.fetch_block_plus_one<7> (er2);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // apply ER feedback
      er1[i] = (q0_15r) (er2[i] * num {0.2});
      er1[i] = (q0_15r) ((late_in[i] + er1[i]) * par.decay[i]);
    }
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run_mod<4> (er1, lfo2);
    rev.run_lp<5> (er1, q0_15::from_float (0.0001f + 0.17f * _param.damp));
    xspan_memcpy (er1b, er1);
    rev.run_mod<6> (er1b, lfo1);
    rev.push<7> (er1b.to_const()); // feedback point

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
    rev.fetch_block_plus_one<14> (l);
    rev.fetch_block_plus_one<22> (r);

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
      late[i] = (q0_15r) (late_in[i] + r[i] * par.decay[i]);
      // add ER blend to late in
      auto fact  = (q0_15) (par.er[i] * num {0.4});
      auto er_in = (q0_15r) ((er1[i] + er2[i]) * fact);
      late[i]    = (q0_15r) (late[i] - er_in);

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

    rev.run_mod<8> (late.cast<q0_15r>(), lfo1, g);
    rev.run<9> (late);
    rev.run_lp<10> (late, late_damp);
    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run<11> (late, g);
    rev.run<12> (late);
    rev.run<13> (late);
    rev.push<14> (late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      // prepare input with feedback
      late[i] = (q0_15r) (late_in[i] + l[i] * par.decay[i]);

      auto fact  = (q0_15) (par.er[i] * num {0.4});
      auto er_in = (q0_15r) ((er1[i] - er2[i]) * fact);
      late[i]    = (q0_15r) (late[i] + er_in);
    }
    l.cut_head (1); // drop feedback sample from previous block

    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run_mod<15> (late.cast<q0_15r>(), lfo2, g);
    rev.run<16> (late);
    rev.run_lp<17> (late, late_damp);
    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run<18> (late, g);
    rev.run<19> (late);
    rev.run<20> (late);
    rev.run_hp<21> (late);
    rev.push<22> (late.to_const()); // feedback point

    // Mixdown
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
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
  struct loop_parameters_flt;
  void process_rev1_flt (
    xspan<std::array<float, 2>> io,
    loop_parameters_flt&        par)
  {
    auto& rev = std::get<rev1_type> (_modes);

    using arr    = std::array<float, max_block_size>;
    using arr_fb = std::array<float, max_block_size + 1>;

    arr late_in_arr;
    arr lfo1;
    arr lfo2;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (io[i][0] + io[i][1]) * 0.5f;
      // ER lfo
      auto lfo = _lfo_er.tick_sine();
      lfo1[i]  = lfo[0];
      lfo2[i]  = lfo[1] * (0.5f + par.er[i] * 0.5f);
      // decay fixup
      auto decay   = 1.f - par.decay[i];
      decay        = 1.f - decay * decay;
      par.decay[i] = 0.6f + decay * 0.39f;
    }

    // diffusion -----------------------------
    rev.run<0> (late_in);
    rev.run<1> (late_in);
    rev.run<2> (late_in);
    rev.run<3> (late_in);

    // ER -----------------------------
    arr    early1_arr;
    arr    early1b_arr;
    arr_fb early2_arr;

    auto er1  = xspan {early1_arr.data(), io.size()};
    auto er1b = xspan {early1b_arr.data(), io.size()};
    auto er2 = xspan {early2_arr.data(), io.size() + 1}; // +1: Feedback on head

    // feedback handling
    rev.fetch_block_plus_one<7> (er2);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      er1[i] = (late_in[i] + er2[i] * 0.2f) * par.decay[i]; // apply ER feedback
    }
    er2.cut_head (1); // drop feedback sample from previous block

    rev.run_mod<4> (er1, lfo1);
    rev.run_lp<5> (er1, 0.0001f + 0.17f * _param.damp);
    xspan_memcpy (er1b, er1);
    rev.run_mod<6> (er1b, lfo2);
    rev.push<7> (er1b.to_const()); // feedback point

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
    rev.fetch_block_plus_one<14> (l);
    rev.fetch_block_plus_one<22> (r);

    for (uint i = 0; i < io.size(); ++i) {
      auto mod = 0.25f + (1.f - par.mod[i]) * 0.75f;
      auto lfo = _lfo.tick_sine();
      lfo1[i]  = lfo[0] * mod;
      lfo2[i]  = lfo[1] * mod;
      late[i]  = late_in[i] + r[i] * par.decay[i];
      late[i] -= (er1[i] + er2[i]) * par.er[i] * 0.4f;
      g[i] = 0.618f + par.character[i] * (0.707f - 0.618f) * 2.f;
    }
    r.cut_head (1); // drop feedback sample from previous block

    float late_damp = (0.9f - _param.damp * 0.9f);
    late_damp       = 1.f - late_damp * late_damp;
    late_damp *= 0.4f;

    rev.run_mod<8> (late, lfo1, g);
    rev.run<9> (late);
    rev.run_lp<10> (late, late_damp);
    std::for_each (g.begin(), g.end(), [] (auto& v) { v = -v; }); // g negate
    rev.run<11> (late, g);
    rev.run<12> (late);
    rev.run<13> (late);
    rev.push<14> (late.to_const()); // feedback point

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      late[i] = late_in[i] + l[i] * par.decay[i];
      late[i] += (er1[i] + er2[i]) * par.er[i] * 0.4f;
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
    rev.run_hp<21> (late);
    rev.push<22> (late.to_const()); // feedback point

    // Mixdown
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      io[i][0] = l[i] + (-er1[i] * 0.825f - er2[i] * 0.423f) * par.er[i];
      io[i][1] = r[i] + (-er1[i] * 0.855f + er2[i] * 0.443f) * par.er[i];
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
    float tilt;
    float predelay;
    float damp;
    float ducking_threshold;
    float ducking_speed;
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
  struct loop_parameters_flt {
    std::array<float, max_block_size> er;
    std::array<float, max_block_size> decay;
    std::array<float, max_block_size> character;
    std::array<float, max_block_size> mod;
    std::array<float, max_block_size> stereo;
  };
  //----------------------------------------------------------------------------
  unsmoothed_parameters                      _param;
  value_smoother<float, smoothed_parameters> _param_smooth;

  block_resampler<float, 2>              _resampler {};
  part_class_array<tilt_eq, f32_x2>      _tilt {};
  static_delay_line<s16, false, false>   _predelay;
  static_delay_line<float, false, false> _predelay_flt;

  lfo<2> _lfo;
  lfo<2> _lfo_er;

  using rev1_type = detail::lofiverb::reverb_tool<detail::lofiverb::rev1_spec>;
  std::variant<rev1_type> _modes;
  ducker<f32_x2>          _ducker;

  uint  _n_processed_samples;
  float _1_4beat_spls;

  std::vector<s16, overaligned_allocator<s16, 16>> _mem;
  xspan<s16>                                       _mem_reverb;
#if 0
  white_noise_generator                            _whitenoise;
#endif
};
//------------------------------------------------------------------------------
} // namespace artv
