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
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  q0_15          g; // If 0 the delay has no allpass
};
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::delay_data make_dd (
  u16   spls,
  float g   = 0.,
  u16   mod = 0)
{
  assert (spls < fixpt_spls::max_int());
  assert (mod < fixpt_spls_mod::max_int());
  return delay_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    q0_15::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr auto get_rev1_spec()
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
  static constexpr auto values {get_rev1_spec()};
};
//------------------------------------------------------------------------------
template <class Spec, uint Max_block_size>
class engine {
public:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
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
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (dst);
    decode_read (dst, get_read_buffers<Idx> (dst.size(), 1));
  }
  // see comment on fetch_block_plus_one
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (dd.g.value() == 0, "Not possible on allpasses");
    static_assert (dd.mod.value() == 0, "Not possible on Modulated delays");
    assert (src);
    encode_write (prepare_block_insertion<Idx> (src.size()), src);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_lp (xspan<T> io, U g)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);
    T y1;
    decode_read (y1, _stage[Idx].z[0]);
    for (uint i = 0; i < io.size(); ++i) {
      T v;
      if constexpr (is_fixpt_v<T>) {
        auto gv = (g.value() == 0) ? dd.g : g;
        v       = (T) ((gv.max() - gv) * io[i]);
        v       = (T) (v + y1 * gv);
      }
      else {
        float gv = (g == 0.f) ? dd.g.to_float() : g;
        v        = (1.f - gv) * io[i];
        v        = v + y1 * gv;
      }
      y1    = v;
      io[i] = v;
    }
    encode_write (_stage[Idx].z[0], y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (xspan<T> io)
  {
    run_lp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_hp (xspan<T> io, U g)
  {
    constexpr delay_data dd = Spec::values[Idx];
    static_assert (get_delay_size (Idx) >= 1);

    assert (io);
    T y1;
    decode_read (y1, _stage[Idx].z[0]);
    for (uint i = 0; i < io.size(); ++i) {
      T v;
      if constexpr (is_fixpt_v<T>) {
        auto gv = (g.value() == 0) ? dd.g : g;
        v       = (T) ((gv.max() - gv) * io[i]);
        v       = (T) (v + y1 * gv);
      }
      else {
        float gv = (g == 0.f) ? dd.g.to_float() : g;
        v        = (1.f - gv) * io[i];
        v        = v + y1 * gv;
      }
      y1    = v;
      io[i] = io[i] - v;
    }
    encode_write (_stage[Idx].z[0], y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (xspan<T> io)
  {
    run_hp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U = T>
  void run (xspan<T> io, xspan<U> g = xspan<U> {})
  {
    run_impl<Idx> (io, g, xspan<U> {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  auto run_mod (xspan<T> io, xspan<U> lfo, xspan<U> g = xspan<U> {})
  {
    run_impl<Idx> (io, g, lfo);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U, std::size_t N>
  auto run_mod (xspan<T> io, std::array<U, N>& lfo, xspan<U> g = xspan<U> {})
  {
    run_impl<Idx> (io, g, xspan {lfo});
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  void decode_read (T& dst, s16 src)
  {
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      static_assert (sizeof dst == sizeof src);
      memcpy (&dst, &src, sizeof src);
    }
    else {
      static_assert (std::is_same_v<T, float>);
      dst = float16::decode (src);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void decode_read (xspan<T> dst, std::array<xspan<s16>, 2> src)
  {
    assert (dst.size() == (src[0].size() + src[1].size()));
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      static_assert (sizeof (T) == sizeof src[0][0]);
      xspan_memdump (dst.data(), src[0]);
      xspan_memdump (dst.data() + src[0].size(), src[1]);
    }
    else {
      for (uint i = 0; i < src[0].size(); ++i) {
        dst[i] = float16::decode (src[0][i]);
      }
      dst.cut_head (src[0].size());
      for (uint i = 0; i < src[1].size(); ++i) {
        dst[i] = float16::decode (src[1][i]);
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (s16& dst, T src)
  {
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      static_assert (sizeof dst == sizeof src);
      memcpy (&dst, &src, sizeof src);
    }
    else {
      static_assert (std::is_same_v<T, float>);
      dst = float16::encode (src);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (std::array<xspan<s16>, 2> dst, xspan<T const> src)
  {
    assert (src.size() == (dst[0].size() + dst[1].size()));
    if constexpr (is_fixpt_v<T>) {
      // 16 bit fixed point
      static_assert (sizeof (T) == sizeof dst[0][0]);
      xspan_memdump (dst[0].data(), src.reduced (dst[1].size()));
      xspan_memdump (dst[1].data(), src.advanced (dst[0].size()));
    }
    else {
      for (uint i = 0; i < dst[0].size(); ++i) {
        dst[0][i] = float16::encode (src[i]);
      }
      src.cut_head (dst[0].size());
      for (uint i = 0; i < dst[1].size(); ++i) {
        dst[1][i] = float16::encode (src[i]);
      }
    }
  }
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
  T read_next()
  {
    constexpr delay_data dd = Spec::values[Idx];
    return get<Idx, T> (dd.spls.to_int());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (T v)
  {
    encode_write (_stage[Idx].z[_stage[Idx].pos], v);
    advance_pos();
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  void run_thiran (xspan<T> dst, xspan<q0_15> lfo)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    T y1;
    decode_read (y1, _stage[Idx].z[y1_pos]);
    for (uint i = 0; i < dst.size(); ++i) {
      auto fixpt_spls
        = (dd.spls.add_sign() + (lfo[i] * dd.mod.add_sign())) - num {i};
      auto n_spls      = fixpt_spls.to_int();
      auto n_spls_frac = fixpt_spls.fractional().to_dynamic();
      auto d           = n_spls_frac + num {0.418f};
      auto one         = fixpt_dt<1, 1, 0>::from_int (1);

      auto co_thiran = (one - d) / (one + d); // 0.4104 to -1
      auto a         = co_thiran.cast<fixpt<1, 0, 23>>();

      auto z0 = get<Idx, T> (n_spls - 1);
      auto z1 = get<Idx, T> (n_spls);
      y1      = (T) (z0 * a + z1 - a * y1);
      dst[i]  = y1;
    }
    encode_write (_stage[Idx].z[y1_pos], y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_thiran (xspan<float> dst, xspan<float> lfo)
  {
    constexpr delay_data dd     = Spec::values[Idx];
    constexpr auto       sz     = get_delay_size (Idx);
    constexpr auto       y1_pos = sz;

    float y1;
    decode_read (y1, _stage[Idx].z[y1_pos]);
    for (uint i = 0; i < dst.size(); ++i) {
      float fixpt_spls  = dd.spls.to_float() + (lfo[i] * dd.mod.to_float());
      uint  n_spls      = (uint) fixpt_spls;
      float n_spls_frac = fixpt_spls - n_spls;
      n_spls -= i;

      float d = n_spls_frac + 0.418f;
      float a = (1.f - d) / (1.f + d); // 0.4104 to -1

      float z0 = get<Idx, float> (n_spls - 1);
      float z1 = get<Idx, float> (n_spls);
      y1       = z0 * a + z1 - a * y1;
      dst[i]   = y1;
    }
    encode_write (_stage[Idx].z[y1_pos], y1);
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
    T ret;
    decode_read (ret, _stage[Idx].z[z]);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  auto run_allpass (T in, T yn, U g)
  {
    auto u = in + yn * g;
    auto x = yn - u * g;
    return std::make_tuple (static_cast<T> (x), static_cast<T> (u));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_impl (xspan<T> io, xspan<U> g, xspan<U> lfo)
  {
    constexpr delay_data dd             = Spec::values[Idx];
    constexpr auto       sz             = get_delay_size (Idx);
    constexpr auto       minsz          = sz - dd.mod.to_int();
    constexpr auto       has_allpass    = dd.g.value() != 0;
    constexpr auto       has_modulation = dd.mod.value() != 0;

    if constexpr (minsz >= max_block_size) {
      // no overlap, can run in blocks
      std::array<T, max_block_size> iocp_mem;
      xspan                         iocp {iocp_mem.data(), io.size()};
      xspan_memcpy (iocp, io);

      if constexpr (has_modulation) {
        run_thiran<Idx> (io, lfo);
      }
      else {
        decode_read (io, get_read_buffers<Idx> (io.size(), 0));
      }
      if constexpr (has_allpass) {
        for (uint i = 0; i < io.size(); ++i) {
          auto [out, push]
            = run_allpass<Idx, T> (iocp[i], io[i], g ? g[i] : (U) dd.g);
          io[i]   = out;
          iocp[i] = push;
        }
      }
      encode_write (prepare_block_insertion<Idx> (io.size()), iocp.to_const());
    }
    else {
      // overlap, needs single sample iteration
      for (uint i = 0; i < io.size(); ++i) {
        T qv;
        if constexpr (has_modulation) {
          qv = run_thiran<Idx> (xspan {&io[i], 1}, xspan {&lfo[i], 1});
        }
        else {
          T qv = read_next<Idx>();
        }
        if constexpr (has_allpass) {
          auto [out, push]
            = run_allpass<Idx, T> (io[i], qv, g ? g[i] : (U) dd.g);
          push<Idx> (push);
          io[i] = out;
        }
        else {
          push<Idx> (io[i]);
          io[i] = qv;
        }
      }
    }
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
  template <uint Idx>
  std::array<xspan<s16>, 2> get_read_buffers (uint blocksize, uint neg_offset)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= blocksize);

    uint block1 = _stage[Idx].pos - dd.spls.to_int() - neg_offset;
    block1 += (block1 >= sz) ? sz : 0;
    uint end = _stage[Idx].pos - dd.spls.to_int() + blocksize - neg_offset - 1;
    end += (end >= sz) ? sz : 0;

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = blocksize - block1sz;
    }
    return {
      xspan {&_stage[Idx].z[block1], block1sz},
      xspan {&_stage[Idx].z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
  // moves the write pointer and returns the locations to write
  template <uint Idx>
  std::array<xspan<s16>, 2> prepare_block_insertion (uint blocksize)
  {
    constexpr delay_data dd = Spec::values[Idx];
    constexpr auto       sz = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= blocksize);

    uint block1 = _stage[Idx].pos;
    uint end    = block1 + blocksize - 1;
    end -= (end >= sz) ? sz : 0;
    _stage[Idx].pos = end;
    advance_pos<Idx>();

    uint block1sz, block2sz;
    if (block1 < end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = sz - block1;
      block2sz = blocksize - block1sz;
    }
    return {
      xspan {&_stage[Idx].z[block1], block1sz},
      xspan {&_stage[Idx].z[0], block2sz}};
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
    // reminder. "normalization_gain" is the gain such that if the reverb is fed
    // values in the [-1,1] range will make it output values at most at the
    // [-1,1] range. This is important for fixed-point and truncated float
    // scaling.
    switch (v) {
    case mode::rev1_flt:
    case mode::rev1: {
      _param.normalization_gain = 1.f / (2.f + 0.8f + 0.855f + 0.443f);
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
  struct freq_balace_tag {};
  void set (freq_balace_tag, float v)
  {
    v *= 0.01f;
    if (v == _param.tilt) {
      return;
    }
    _param.tilt = v;
    auto db     = vec_set<2> ((float) v * -14.f);
    _filt.reset_coeffs<0> (vec_set<2> (350.f), vec_set<2> (0.5f), db, t_spl);
    _filt.reset_coeffs<1> (vec_set<2> (2200.f), db * -0.5f);
  }

  static constexpr auto get_parameter (freq_balace_tag)
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
    return float_param ("dB", -60.f, 0.0f, 0.f, 0.01f);
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

    _filt.reset_states_cascade();
    _ducker.reset();
    _lfo.reset();
    _lfo.set_phase (phase<4> {phase_tag::normalized {}, 0.f, 0.5f, 0.f, 0.5f});

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
        set (param, get_parameter (param).max);
      }
      else {
        // max might be not yet impl.
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
    freq_balace_tag,
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

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wetv    = _filt.tick_cascade (vec_from_array (io[i]));
      wetv           = vec_clamp (wetv, -0.98f, 0.98f);
      ducker_gain[i] = _ducker.tick (wetv);
      wetv *= _param.normalization_gain;
      io[i] = vec_to_array (wetv);

      _param_smooth.tick();
      pars.stereo[i] = _param_smooth.get().stereo;
      pars.er[i].load_float (_param_smooth.get().er);
      pars.decay[i].load_float (_param_smooth.get().decay);
      pars.character[i].load_float (_param_smooth.get().character);
      pars.mod[i].load_float (_param_smooth.get().mod);
    }
    // predelay + int conversion
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        // TODO: block fetch?
        wet[i][0].load (_predelay.get (predelay_spls, 0));
        wet[i][1].load (_predelay.get (predelay_spls, 1));
        auto push = make_array (
          fixptype::from_float (io[i][0]).value(),
          fixptype::from_float (io[i][1]).value());
        _predelay.push (xspan {push});
      }
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        wet[i][0].load_float (io[i][0]);
        wet[i][1].load_float (io[i][0]);
      }
    }
    // main loop
    process_rev1 (xspan {wet.data(), io.size()}, pars);

    // float conversion
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = wet[i][0].to_float() * ducker_gain[i][0];
      auto r = wet[i][1].to_float() * ducker_gain[i][1];
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
    std::array<f32_x2, max_block_size> ducker_gain;
    loop_parameters_flt                pars;

    // tilt + clamp + ducker measuring + param smoothing
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      f32_x2 wet     = _filt.tick_cascade (vec_from_array (io[i]));
      wet            = vec_clamp (wet, -0.98f, 0.98f);
      ducker_gain[i] = _ducker.tick (wet);
      wet *= _param.normalization_gain;
      io[i] = vec_to_array (wet);

      _param_smooth.tick();
      pars.stereo[i]    = _param_smooth.get().stereo;
      pars.er[i]        = _param_smooth.get().er;
      pars.decay[i]     = _param_smooth.get().decay;
      pars.character[i] = _param_smooth.get().character;
      pars.mod[i]       = _param_smooth.get().mod;
    }
    // predelay
    if (_param.predelay != 0) {
      uint predelay_spls = _1_4beat_spls * _param.predelay;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        // TODO: block fetch?
        std::array<float, 2> spl;
        spl[0] = _predelay_flt.get (predelay_spls, 0);
        spl[1] = _predelay_flt.get (predelay_spls, 1);
        _predelay_flt.push (xspan {io[i]});
        io[i] = spl;
      }
    }
    // main loop
    process_rev1_flt (io, pars);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto l = io[i][0] * ducker_gain[i][0];
      auto r = io[i][1] * ducker_gain[i][1];
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
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = ((io[i][0] + io[i][1]) >> 1).resize<-1>();

      auto mod = (q1_14) (num {0.25} + (num {1} - par.mod[i]) * num {0.75});
      // ER + late lfo
      auto lfo = _lfo.tick_sine_fixpt().spec_cast<q0_15>().value();
      lfo1[i].load (lfo[0]);
      lfo2[i].load (lfo[1]);
      auto er_amt    = (q1_14) (num {0.5} + (par.er[i] >> 1));
      auto lfo2final = lfo2[i] * er_amt;
      lfo2[i]        = (q0_15) (lfo2final);
      lfo3[i].load (lfo[2]);
      lfo3[i] = (q0_15) (lfo3[i] * mod);
      lfo4[i].load (lfo[3]);
      lfo4[i] = (q0_15) (lfo4[i] * mod);

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

    float late_dampf = (0.9f - _param.damp * 0.9f);
    late_dampf       = 1.f - late_dampf * late_dampf;
    late_dampf *= 0.4f;
    auto late_damp = q0_15::from_float (late_dampf);

    rev.run_mod<8> (late.cast<q0_15r>(), lfo3, g);
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
    rev.run_mod<15> (late.cast<q0_15r>(), lfo4, g);
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
    arr lfo3;
    arr lfo4;

    auto late_in = xspan {late_in_arr.data(), io.size()};
    for (uint i = 0; i < io.size(); ++i) {
      // to MS
      late_in[i] = (io[i][0] + io[i][1]) * 0.5f;
      // ER + late lfo
      auto mod = 0.25f + (1.f - par.mod[i]) * 0.75f;
      auto lfo = _lfo.tick_sine();
      lfo1[i]  = lfo[0];
      lfo2[i]  = lfo[1] * (0.5f + par.er[i] * 0.5f);
      lfo3[i]  = lfo[2] * mod;
      lfo4[i]  = lfo[3] * mod;
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
      late[i] = late_in[i] + r[i] * par.decay[i];
      late[i] -= (er1[i] + er2[i]) * par.er[i] * 0.4f;
      g[i] = 0.618f + par.character[i] * (0.707f - 0.618f) * 2.f;
    }
    r.cut_head (1); // drop feedback sample from previous block

    float late_damp = (0.9f - _param.damp * 0.9f);
    late_damp       = 1.f - late_damp * late_damp;
    late_damp *= 0.4f;

    rev.run_mod<8> (late, lfo3, g);
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
    rev.run_mod<15> (late, lfo4, g);
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
    // Reminder, phase are at 0, 180, 0, 180 on reset.
    auto mod    = _param_smooth.target().mod;
    auto f_er   = 0.3f + mod * 0.3f;
    auto f_late = 0.1f + mod * 1.2f;
    _lfo.set_freq (f32_x4 {f_er, f_er, f_late, f_late}, t_spl);
  }
  //----------------------------------------------------------------------------
  struct unsmoothed_parameters {
    u32   mode;
    float normalization_gain;
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

  block_resampler<float, 2>                                       _resampler {};
  part_classes<mp_list<tilt_eq, onepole_naive_highshelf>, f32_x2> _filt {};
  static_delay_line<s16, true, false>                             _predelay;
  static_delay_line<float, true, false>                           _predelay_flt;

  lfo<4> _lfo;

  using rev1_type
    = detail::lofiverb::engine<detail::lofiverb::rev1_spec, max_block_size>;
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