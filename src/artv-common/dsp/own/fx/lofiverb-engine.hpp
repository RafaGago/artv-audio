#pragma once

// internals of lofiverb.hpp

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <variant>
#include <vector>

#include "artv-common/dsp/own/classes/noise.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/float.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
// truncating fixed point type for computation
using fixpt_tt = fixpt_d<1, 0, 15, 0>;
// rounding fixed point type for computation
using fixpt_tr = fixpt_d<1, 0, 15, fixpt_rounding>;

using fixpt_t = fixpt_tt;

using fixpt_acum_t = fixpt_s<1, 0, 31>;

// fixed point type for storage
using fixpt_sto      = fixpt_s<1, 0, 15, 0>;
using fixpt_spls     = fixpt_m<0, 14, 0, 0>;
using fixpt_spls_mod = fixpt_m<0, 9, 0, 0>;
//------------------------------------------------------------------------------
struct allpass_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
};

struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
};

struct filter_data {
  fixpt_t g;
  bool    is_lowpass;
};

struct quantizer {
  fixpt_t g;
};

using stage_data
  = std::variant<allpass_data, delay_data, filter_data, quantizer>;
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_ap (
  u16   spls,
  float g,
  u16   mod = 0)
{
  return allpass_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_delay (u16 spls, u16 mod = 0)
{
  return delay_data {
    fixpt_spls::from_int (spls), fixpt_spls_mod::from_int (mod)};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_hp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), false};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_lp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), true};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_quantizer()
{
  return quantizer {};
}
//------------------------------------------------------------------------------

template <class Spec_array>
class spec_access {
public:
  //----------------------------------------------------------------------------
  static constexpr bool is_allpass (uint i)
  {
    return std::holds_alternative<allpass_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_delay (uint i)
  {
    return std::holds_alternative<delay_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_filter (uint i)
  {
    return std::holds_alternative<filter_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_quantizer (uint i)
  {
    return std::holds_alternative<quantizer> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_lowpass_filter (uint i)
  {
    if (is_filter (i)) {
      return std::get<filter_data> (values[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_t get_gain (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (values[i]).g;
    }
    else if (is_filter (i)) {
      return std::get<filter_data> (values[i]).g;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls get_delay_spls (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (values[i]).spls;
    }
    else if (is_delay (i)) {
      return std::get<delay_data> (values[i]).spls;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls_mod get_delay_mod_spls (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (values[i]).mod;
    }
    else if (is_delay (i)) {
      return std::get<delay_data> (values[i]).mod;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool has_modulated_delay (uint i)
  {
    return get_delay_mod_spls (i).value() != 0;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_max_delay_spls (uint i)
  {
    return get_delay_spls (i).to_int() + get_delay_mod_spls (i).to_int();
  }
  //----------------------------------------------------------------------------
  static constexpr std::size_t size() { return values.size(); }
  //----------------------------------------------------------------------------
private:
  static constexpr auto values {Spec_array::values};
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// A class to abstract 16-bit storage, queue access and common DSP operations
// when building reverbs based on allpass loops. Both on fixed and floating
// point.
template <class Spec_array, uint Max_block_size>
class engine {
public:
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
  //----------------------------------------------------------------------------
  static constexpr uint get_required_size()
  {
    uint ret = 0;
    for (uint i = 0; i < decltype (_stage) {}.size(); ++i) {
      ret += get_stage_n_elems (i);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<s16> mem)
  {
    for (uint i = 0; i < _stage.size(); ++i) {
      _stage[i].z = mem.cut_head (get_stage_n_elems (i)).data();
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
  void fetch_block (xspan<T> dst, uint negative_offset = 0)
  {
    static_assert (spec::is_delay (Idx), "Only possible on pure delays");
    static_assert (
      !spec::has_modulated_delay (Idx), "Not possible on Modulated delays");
    assert (dst);
    decode_read (dst, get_read_buffers<Idx> (dst.size(), negative_offset));
  }
  // see comment on fetch_block_plus_one
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    static_assert (spec::is_delay (Idx), "Only possible on pure delays");
    static_assert (
      !spec::has_modulated_delay (Idx), "Not possible on Modulated delays");
    assert (src);
    encode_write (prepare_block_insertion<Idx> (src.size()), src);
  }
  //----------------------------------------------------------------------------
  template <uint Start_idx, uint N, class T>
  void run_multichannel (std::array<T*, N> io, uint block_size)
  {
    mp11::mp_for_each<mp11::mp_iota_c<N>> ([&] (auto idx) {
      run<Start_idx + idx.value> (xspan {io[idx.value], block_size});
    });
  }
  //----------------------------------------------------------------------------
  template <class T, uint N>
  static void hadamard_no_norm (std::array<T*, N> io, uint block_size)
  {
    static_assert (is_pow2 (N));
    if constexpr (N > 1) {
      constexpr int half_n = N / 2;
      // all this superfluous array copying stuff should be removed by the
      // optimizer.
      hadamard (array_slice<0, half_n> (io), block_size);
      hadamard (array_slice<half_n, half_n> (io), block_size);

      for (uint n = 0; n < half_n; ++n) {
        uint n1 = n;
        uint n2 = n + half_n;
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < block_size; ++i) {
          T a       = io[n1][i];
          T b       = io[n2][i];
          io[n1][i] = (a + b);
          io[n2][i] = (a - b);
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T, uint N>
  static void hadamard (std::array<T*, N> io, uint block_size)
  {
    hadamard_no_norm (io, block_size);
    // normalize by 1/sqrtn, not using GCEM because this benefits from fixed
    // point constant literals.
    if constexpr (N == 2) {
      for (auto iov : io) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < block_size; ++i) {
          io[i] = (T) (io[i] * get_hadamard_ratio<N>());
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class... Ts>
  void run (xspan<T> io, Ts&&... args)
  {
    if constexpr (spec::is_allpass (Idx) || spec::is_delay (Idx)) {
      run_ap_or_delay<Idx> (io, std::forward<Ts> (args)...);
    }
    else if constexpr (spec::is_filter (Idx)) {
      if constexpr (spec::is_lowpass_filter (Idx)) {
        run_lp<Idx> (io, std::forward<Ts> (args)...);
      }
      else {
        run_hp<Idx> (io, std::forward<Ts> (args)...);
      }
    }
    else if constexpr (spec::is_quantizer (Idx)) {
      run_quantizer<Idx> (io, std::forward<Ts> (args)...);
    }
    else {
      static_assert (sizeof (T) != sizeof (T), "Invalid");
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx1, uint Idx2, class T, class... Ts>
  void run (xspan<T> io, Ts&&... args)
  {
    if constexpr (
      spec::is_allpass (Idx1)
      && (spec::is_allpass (Idx2) || spec::is_delay (Idx2))) {
      run_nested_ap2<Idx1, Idx2> (io, std::forward<Ts> (args)...);
    }
    else {
      static_assert (sizeof (T) != sizeof (T), "Invalid");
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T, class Lfo>
  static constexpr auto get_lfo_generator (Lfo&& v)
  {
    if constexpr (
      std::is_convertible_v<Lfo, std::function<T (uint)>>
      && !std::is_same_v<Lfo, std::nullptr_t>) {
      return std::forward<Lfo> (v);
    }
    else if constexpr (is_xspan<Lfo>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return nullptr;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Gain>
  static constexpr auto get_gain_generator (Gain&& v)
  {
    if constexpr (
      std::is_convertible_v<Gain, std::function<T (uint)>>
      && !std::is_same_v<Gain, std::nullptr_t>) {
      return std::forward<Gain> (v);
    }
    else if constexpr (is_xspan<Gain>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return [] (uint) { return (T) spec::get_gain (Idx); };
    }
  }
  //----------------------------------------------------------------------------
  template <uint N>
  static constexpr auto get_hadamard_ratio()
  {
    if constexpr (N == 4) {
      return 0.5_r;
    }
    else if constexpr (N == 8) {
      return 0.35355339059327373_r;
    }
    else if constexpr (N == 16) {
      return 0.25_r;
    }
    else if constexpr (N == 32) {
      return 0.17677669529663687_r;
    }
    else if constexpr (N == 64) {
      return 0.125_r;
    }
    else {
      static_assert (N == N, "Unimplemented");
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void decode_read (T& dst, s16 src)
  {
    if constexpr (std::is_same_v<T, fixpt_t>) {
      dst = fixpt_sto::from (src);
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
    for (uint i = 0; i < src[0].size(); ++i) {
      decode_read (dst[i], src[0][i]);
    }
    dst.cut_head (src[0].size());
    for (uint i = 0; i < src[1].size(); ++i) {
      decode_read (dst[i], src[1][i]);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (s16& dst, T src)
  {
    if constexpr (std::is_same_v<T, fixpt_t>) {
      dst = ((fixpt_sto) src).value();
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
    for (uint i = 0; i < dst[0].size(); ++i) {
      encode_write (dst[0][i], src[i]);
    }
    src.cut_head (dst[0].size());
    for (uint i = 0; i < dst[1].size(); ++i) {
      encode_write (dst[1][i], src[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void advance_pos()
  {
    constexpr auto sz = get_delay_size (Idx);

    ++_stage[Idx].pos;
    if (_stage[Idx].pos == sz) {
      _stage[Idx].pos = 0;
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T read_next()
  {
    constexpr auto spls = spec::get_delay_spls (Idx).to_int();
    return get<Idx, T> (spls);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push_one (T v)
  {
    encode_write (_stage[Idx].z[_stage[Idx].pos], v);
    advance_pos<Idx>();
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_lp (xspan<float> io, U g)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);
    float y1 = fetch_state<Idx> (float {});
    for (uint i = 0; i < io.size(); ++i) {
      float gv = (g == 0.f) ? (float) spec::get_gain (Idx) : g;
      float y  = y1 * gv + (1_r - gv) * io[i];
      io[i]    = y;
      y1       = y;
    }
    save_state<Idx> (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (xspan<fixpt_t> io, Func&& fn)
  {
    // decay with error-feedback/ fraction saving
    constexpr uint mask
      = lsb_mask<uint> (fixpt_acum_t::n_bits - fixpt_t::n_bits);

    s16 err = fetch_quantizer_err<Idx>();
    for (uint i = 0; i < io.size(); ++i) {
      auto v = (fixpt_acum_t) fn (io[i], i);
      v += fixpt_acum_t::from (err);
      io[i] = (fixpt_tr) v;
      err   = (s16) (v.value() & mask);
    }
    save_quantizer_err<Idx> (err);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (xspan<float> io, Func&& fn)
  {
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = fn (io[i], i);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_lp (xspan<fixpt_t> io, U g)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);
    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    constexpr uint mask
      = lsb_mask<uint> (fixpt_acum_t::n_frac - fixpt_t::n_frac);
    auto [y1v, errv] = fetch_state<Idx> (fixpt_t {});
    auto y1          = fixpt_t::from (y1v);
    auto err         = fixpt_acum_t::from (errv);

    for (uint i = 0; i < io.size(); ++i) {
      fixpt_t gvf = (g.value() == 0) ? spec::get_gain (Idx) : g;
      auto    gv  = (fixpt_acum_t) gvf;
      auto    x   = (fixpt_acum_t) io[i];
#if 0
// investigate why it does sound better without noise-shaping...
      auto    y   = gv * y1 + (1_r - gv) * x + err;
      y1          = (fixpt_t) y;
      err         = fixpt_acum_t::from (y.value() & mask);
#else
      auto y = gv * y1 + (1_r - gv) * x;
      y1     = (fixpt_t) y;
#endif
      io[i] = y1;
    }
    save_state<Idx> (y1.value(), err.value());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (xspan<T> io)
  {
    static_assert (spec::is_filter (Idx));
    run_lp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_hp (xspan<float> io, U g)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);
    float y1 = fetch_state<Idx> (float {});
    for (uint i = 0; i < io.size(); ++i) {
      float gv = (g == 0.f) ? (float) spec::get_gain (Idx) : g;
      float y  = io[i] * gv;
      y += (1_r - gv) * y1;
      io[i] = y;
      y1    = y;
    }
    save_state<Idx> (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_hp (xspan<fixpt_t> io, U g)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);
    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    constexpr uint mask
      = lsb_mask<uint> (fixpt_acum_t::n_frac - fixpt_t::n_frac);
    auto [y1v, errv] = fetch_state<Idx> (fixpt_t {});
    auto y1          = fixpt_t::from (y1v);
    auto err         = fixpt_acum_t::from (errv);

    for (uint i = 0; i < io.size(); ++i) {
      fixpt_t gvf = (g.value() == 0) ? spec::get_gain (Idx) : g;
      auto    gv  = (fixpt_acum_t) gvf;
      auto    x   = (fixpt_acum_t) io[i];
      auto    y   = gv * x + (1_r - gv) * y1 + err;
      y1          = (fixpt_t) y;
      err         = fixpt_acum_t::from (y.value() & mask);
      io[i]       = y1;
    }
    save_state<Idx> (y1.value(), err.value());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (xspan<T> io)
  {
    run_hp<Idx> (io, T {});
  }
  //----------------------------------------------------------------------------
  // run an allpass or pure delay.
  // Lfo and gain parameters can be either value generators (lambdas with
  // T(int) signature), spans or dummy types. When dummy types (no span and no
  // lambda types, e.g. nullptr) the parameters are defaulted (no delay
  // modulation, same allpass gain as in the specification)
  template <uint Idx, class T, class L, class G>
  void run_ap_or_delay (xspan<T> io, L lfo, G gain)
  {
    run_ap_or_delay_impl<Idx> (
      io,
      get_gain_generator<Idx, T> (std::forward<G> (gain)),
      get_lfo_generator<T> (std::forward<L> (lfo)));
  }
  //----------------------------------------------------------------------------
  // run an allpass or pure delay
  template <uint Idx, class T>
  void run_ap_or_delay (xspan<T> io)
  {
    run_ap_or_delay<Idx> (io, nullptr, nullptr);
  }
  //----------------------------------------------------------------------------
  // run 2 level nested allpass.
  // Lfo and gain parameters can be either value generators (lambdas with
  // T(int) signature), spans or dummy types. When dummy types (no span and no
  // lambda types, e.g. nullptr) the parameters are defaulted (no delay
  // modulation, same allpass gain as in the specification)
  template <
    uint Idx1,
    uint Idx2,
    class T,
    class L1,
    class G1,
    class L2,
    class G2>
  void run_nested_ap2 (xspan<T> io, L1 lfo1, G1 g1, L2 lfo2, G2 g2)
  {
    run_nested_ap_impl<Idx1, Idx2> (
      io,
      get_gain_generator<Idx1, T> (std::forward<G1> (g1)),
      get_lfo_generator<T> (std::forward<L1> (lfo1)),
      get_gain_generator<Idx2, T> (std::forward<G2> (g2)),
      get_lfo_generator<T> (std::forward<L2> (lfo2)));
  }
  //----------------------------------------------------------------------------
  // run a 2-level nested allpass (notice that more levels are possible on the
  // internal implementation, just not exposed outside)
  template <uint Idx1, uint Idx2, class T>
  void run_nested_ap2 (xspan<T> io)
  {
    run<Idx1, Idx2> (io, nullptr, nullptr, nullptr, nullptr);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<fixpt_t> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).add_sign();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).add_sign();

    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    // requires extra integer room on e.g the "d " calculation
    using fixpt_th      = decltype (fixpt_acum_t {}.resize<2, -2>());
    constexpr uint mask = lsb_mask<uint> (fixpt_th::n_frac - fixpt_t::n_frac);
    auto [y1v, errv]    = fetch_state<Idx> (fixpt_t {});
    auto y1             = fixpt_t::from (y1v);
    auto err            = fixpt_th::from (errv);

    for (uint i = 0; i < dst.size(); ++i) {
      auto fixpt_spls = (delay_spls + (lfo_gen (i) * delay_mod_spls));
      auto n_spls     = fixpt_spls.to_int();
      n_spls -= i;
      auto n_spls_frac = (fixpt_th) fixpt_spls.fractional();
      auto d           = n_spls_frac + 0.418_r; // this might exceed 1
      auto a           = (1_r - d) / (1_r + d); // 0.4104 to -1
      auto z0          = (fixpt_th) get<Idx, fixpt_t> (n_spls - 1);
      auto z1          = (fixpt_th) get<Idx, fixpt_t> (n_spls);
      auto y           = (z0 * a) + z1 - (a * y1) + err;
      y1               = (fixpt_t) y; // truncated
      err              = fixpt_th::from (y.value() & mask); // error
      dst[i]           = y1;
    }
    save_state<Idx> (y1.value(), (s16) err.value());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<float> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).to_floatp();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).to_floatp();

    float y1 = fetch_state<Idx> (float {});
    for (uint i = 0; i < dst.size(); ++i) {
      float fixpt_spls  = (delay_spls + (lfo_gen (i) * delay_mod_spls));
      uint  n_spls      = (uint) fixpt_spls;
      float n_spls_frac = fixpt_spls - n_spls;
      n_spls -= i;

      float d = n_spls_frac + 0.418_r;
      float a = (1_r - d) / (1_r + d); // 0.4104 to -1

      auto z0 = get<Idx, float> (n_spls - 1);
      auto z1 = get<Idx, float> (n_spls);
      y1      = z0 * a + z1 - a * y1;
      dst[i]  = y1;
    }
    save_state<Idx> (y1);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  T get (uint delay_spls)
  {
    constexpr auto sz = get_delay_size (Idx);

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
    static_assert (spec::is_allpass (Idx));
    auto u = in + yn * g;
    if constexpr (std::is_floating_point_v<T>) {
      if (abs (u) >= 1.) {
        assert (false);
      }
    }
    auto x = yn - u * g;
    if constexpr (std::is_floating_point_v<T>) {
      if (abs (x) >= 1.) {
        assert (false);
      }
      return std::make_tuple (x, u);
    }
    else {
#if 0
      // dither + noise shaping
      constexpr uint n_truncated = fixpt_acum_t::n_bits - T::n_bits;
      constexpr uint mask        = lsb_mask<uint> (n_truncated);
      auto [z_err_out, z_err_yn] = fetch_ap_err<Idx>();
      auto noise
        = tpdf_dither<n_truncated, fixpt_acum_t::scalar_type, 2> (_noise);
      // cast to something wide and known
      auto out = (fixpt_acum_t) x + fixpt_acum_t::from (z_err_out);
      auto yn  = (fixpt_acum_t) u + fixpt_acum_t::from (z_err_yn);
      auto d_out = out + fixpt_acum_t::from (noise[0]);
      auto d_yn  = yn + fixpt_acum_t::from (noise[1]);
      auto q_out = (T) d_out;
      auto q_yn  = (T) d_yn;
      auto err_out = out - q_out;
      auto err_yn  = yn - q_yn;
      save_ap_err<Idx> (err_out.value() & mask, err_yn.value() & mask);
      return std::make_tuple (q_out, q_yn);
#else
      // noise shaping only
      constexpr uint n_truncated = fixpt_acum_t::n_bits - T::n_bits;
      constexpr uint mask        = lsb_mask<uint> (n_truncated);
      auto [z_err_out, z_err_yn] = fetch_ap_err<Idx>();
      auto out   = (fixpt_acum_t) x + fixpt_acum_t::from (z_err_out);
      auto yn    = (fixpt_acum_t) u + fixpt_acum_t::from (z_err_yn);
      auto q_out = (T) out;
      auto q_yn  = (T) yn;
      save_ap_err<Idx> (out.value() & mask, yn.value() & mask);
      return std::make_tuple (q_out, q_yn);
#endif
    }
  }
  //----------------------------------------------------------------------------
  // Template function to run:
  // - delays
  // - modulated delays
  // - non-nested allpasses
  // - modulated non-nested allpasses.
  //
  template <uint Idx, class T, class GF, class LF>
  void run_ap_or_delay_impl (xspan<T> io, GF&& g_gen, LF&& lfo_gen)
  {
    static_assert (!spec::is_filter (Idx));

    constexpr auto sz    = get_delay_size (Idx);
    constexpr auto minsz = sz - spec::get_delay_mod_spls (Idx).to_int();

    auto szw = get_delay_size (Idx);

    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      assert (io.size() <= minsz);

      std::array<T, max_block_size> iocp_mem;
      xspan                         iocp {iocp_mem.data(), io.size()};
      xspan_memcpy (iocp, io);

      if constexpr (spec::has_modulated_delay (Idx)) {
        run_thiran<Idx> (io, std::forward<LF> (lfo_gen));
      }
      else {
        decode_read (io, get_read_buffers<Idx> (io.size(), 0));
      }
      if constexpr (spec::is_allpass (Idx)) {
        for (uint i = 0; i < io.size(); ++i) {
          auto [out, push] = run_allpass<Idx, T> (iocp[i], io[i], g_gen (i));
          io[i]            = out;
          iocp[i]          = push;
        }
      }
      encode_write (prepare_block_insertion<Idx> (io.size()), iocp.to_const());
    }
    else {
      // overlap, needs single-sample iteration
      for (uint i = 0; i < io.size(); ++i) {
        T delayed;
        if constexpr (spec::has_modulated_delay (Idx)) {
          // this is for completeness, modulations under the block size are
          // unlikely.
          delayed = io[i];
          run_thiran<Idx> (xspan {&delayed, 1}, [i, &lfo_gen] (uint) {
            return lfo_gen (i);
          });
        }
        else {
          delayed = read_next<Idx, T>();
        }
        if constexpr (spec::is_allpass (Idx)) {
          auto [out, push] = run_allpass<Idx, T> (io[i], delayed, g_gen (i));
          push_one<Idx> (push);
          io[i] = out;
        }
        else {
          push_one<Idx> (io[i]);
          io[i] = delayed;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run (arbitrarily) nested allpasses
  template <
    uint I,
    uint N,
    class Idxs,
    class T,
    class GainGen,
    class LfoGen,
    class... GainAndLfoGens>
  void run_nested_ap_impl_get_yn (
    std::array<xspan<T>, N>& yn,
    GainGen&&                _,
    LfoGen&&                 lfo_gen,
    GainAndLfoGens&&... fwd)
  {
    constexpr uint idx   = mp11::mp_at_c<Idxs, I>::value;
    constexpr uint sz    = get_delay_size (idx);
    constexpr uint minsz = sz - spec::get_delay_mod_spls (idx).to_int();

    // only block processing for now, non-block processing could be added.
    static_assert (minsz >= max_block_size);
    assert (yn.size() <= minsz);

    if constexpr (spec::has_modulated_delay (idx)) {
      run_thiran<idx> (yn[I], std::forward<LfoGen> (lfo_gen));
    }
    else {
      decode_read (yn[I], get_read_buffers<idx> (yn[I].size(), 0));
    }
    // recurse....
    if constexpr ((I + 1) < N) {
      run_nested_ap_impl_get_yn<I + 1, N, Idxs> (
        yn, std::forward<GainAndLfoGens> (fwd)...);
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run nested allpasses
  template <
    uint I,
    uint N,
    class T,
    class GainGen,
    class LfoGen,
    class... GainAndLfoGens>
  void run_nested_ap_impl_get_gains (
    std::array<T, N>& gains,
    uint              n_spl,
    GainGen&&         gain_gen,
    LfoGen&&          _,
    GainAndLfoGens&&... fwd)
  {
    gains[I] = gain_gen (n_spl);
    // recurse....
    if constexpr ((I + 1) < N) {
      run_nested_ap_impl_get_gains<I + 1, N> (
        gains, n_spl, std::forward<GainAndLfoGens> (fwd)...);
    }
  }
  //----------------------------------------------------------------------------
  // Helper template function to run nested allpasses
  template <class T, uint N, class... GainAndLfoGens>
  std::array<T, N> run_nested_ap_impl_get_gains (
    uint n_spl,
    GainAndLfoGens&&... glfo)
  {
    std::array<T, N> ret;
    run_nested_ap_impl_get_gains<0, N> (
      ret, n_spl, std::forward<GainAndLfoGens> (glfo)...);
    return ret;
  }
  //----------------------------------------------------------------------------
  // Template function to run arbitrarily nested allpass/plain delay
  // combinations with or gain and delay time modulation (the gain and time
  // modulations can be unity but have to be present). It only work for delay
  // sizes bigger than the block size for implementation "simplicity".
  template <uint... Idx, class T, class... GainAndLfoGens>
  void run_nested_ap_impl (xspan<T> io, GainAndLfoGens&&... gnlfo)
  {
    using Idxs          = mp_list<std::integral_constant<uint, Idx>...>;
    constexpr auto idx0 = mp11::mp_front<Idxs>::value;
    constexpr uint n    = mp11::mp_size<Idxs>::value;
    static_assert ((n * 2) == mp11::mp_size<mp_list<GainAndLfoGens...>>::value);

    array2d<T, max_block_size, n> yn_mem;
    std::array<xspan<T>, n>       yn;
    for (uint i = 0; i < n; ++i) {
      yn[i] = xspan {yn_mem[i].data(), io.size()};
    }
    run_nested_ap_impl_get_yn<0, n, Idxs> (
      yn, std::forward<GainAndLfoGens> (gnlfo)...);

    for (uint i = 0; i < io.size(); ++i) {
      auto g = run_nested_ap_impl_get_gains<T, n> (
        i, std::forward<GainAndLfoGens> (gnlfo)...);

      std::array<T, n> u;

      mp_foreach_idx (Idxs {}, [&] (auto j, auto stage_idx) {
        if constexpr (stage_idx.value == idx0) {
          static_assert (
            spec::is_allpass (idx0), "1st stage has to be allpass");
          auto [out, push]
            = run_allpass<stage_idx.value, T> (io[i], yn[0][i], g[0]);
          u[0]  = push;
          io[i] = out;
        }
        else {

          // subsequent stages of the lattice can be allpasses or plain delays
          if constexpr (spec::is_allpass (stage_idx.value)) {
            auto [out, push]
              = run_allpass<stage_idx.value, T> (u[0], yn[j][i], g[j]);
            u[0] = push;
            u[j] = out;
          }
          else {
            static_assert (spec::is_delay (stage_idx.value));
            // The optimizer would see the 0 mul anyways, just for
            // highlighting that plain delays can be added to the allpass
            // chain
            u[j] = yn[j][i];
          }
        }
      });
      for (uint j = 0; j < (n - 1); ++j) {
        yn[j][i] = u[j + 1];
      }
      yn[n - 1][i] = u[0];
    }
    mp_foreach_idx (Idxs {}, [&] (auto i, auto topo_stage) {
      encode_write (
        prepare_block_insertion<topo_stage.value> (io.size()),
        yn[i].to_const());
    });
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_size (uint i)
  {
    // + 1 is some legacy that I forgot. I vaguely remember that wasting one
    // element allowed saving code (TODO investigate)
    return spec::get_max_delay_spls (i) + 1;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_states_offset (uint i)
  {
    return get_delay_size (i);
  }
  //----------------------------------------------------------------------------
  static constexpr uint n_thiran_states    = 2;
  static constexpr uint n_filter_states    = 2;
  static constexpr uint n_allpass_states   = 2;
  static constexpr uint n_quantizer_states = 1;
  //----------------------------------------------------------------------------
  static constexpr uint get_stage_n_elems (uint i)
  {
    // ordering as appears here
    uint ret = get_delay_size (i);
    ret += ((uint) spec::has_modulated_delay (i)) * n_thiran_states;
    ret += ((uint) spec::is_filter (i)) * n_filter_states;
    ret += ((uint) spec::is_allpass (i)) * n_allpass_states;
    ret += ((uint) spec::is_quantizer (i)) * n_quantizer_states;
    return ret;
  }
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // fetch state for thiran interpolators or damp filters
  template <uint Idx>
  auto fetch_state (fixpt_t)
  {
    static_assert (spec::has_modulated_delay (Idx) || spec::is_filter (Idx));
    constexpr auto offset = get_states_offset (Idx);
    return std::make_tuple (_stage[Idx].z[offset], _stage[Idx].z[offset + 1]);
  }
  //----------------------------------------------------------------------------
  // fetch state for thiran interpolators or damp filters
  template <uint Idx>
  float fetch_state (float)
  {
    static_assert (spec::has_modulated_delay (Idx) || spec::is_filter (Idx));
    static_assert (sizeof (float) == 2 * sizeof (fixpt_t));
    float          ret;
    constexpr auto offset = get_states_offset (Idx);
    memcpy (&ret, &_stage[Idx].z[offset], sizeof ret);
    return ret;
  }
  //----------------------------------------------------------------------------
  // save state for thiran interpolators or damp filters
  template <uint Idx>
  void save_state (s16 v1, s16 v2)
  {
    static_assert (spec::has_modulated_delay (Idx) || spec::is_filter (Idx));
    constexpr auto offset     = get_states_offset (Idx);
    _stage[Idx].z[offset]     = v1;
    _stage[Idx].z[offset + 1] = v2;
  }
  //----------------------------------------------------------------------------
  // save state for thiran interpolators or damp filters
  template <uint Idx>
  void save_state (float v)
  {
    static_assert (spec::has_modulated_delay (Idx) || spec::is_filter (Idx));
    static_assert (sizeof (v) == 2 * sizeof (fixpt_t));
    constexpr auto offset = get_states_offset (Idx);
    memcpy (&_stage[Idx].z[offset], &v, sizeof v);
  }
  //----------------------------------------------------------------------------
  // fetch truncation error for allpass (fixed point only)
  template <uint Idx>
  auto fetch_ap_err()
  {
    static_assert (spec::is_allpass (Idx));
    constexpr auto offset = get_states_offset (Idx)
      + (spec::has_modulated_delay (Idx) ? n_thiran_states : 0);
    return std::make_tuple (_stage[Idx].z[offset], _stage[Idx].z[offset + 1]);
  }
  //----------------------------------------------------------------------------
  // store truncation error for allpass (fixed point only)
  template <uint Idx>
  void save_ap_err (s16 err1, s16 err2)
  {
    static_assert (spec::is_allpass (Idx));
    constexpr auto offset = get_states_offset (Idx)
      + (spec::has_modulated_delay (Idx) ? n_thiran_states : 0);
    _stage[Idx].z[offset]     = err1;
    _stage[Idx].z[offset + 1] = err2;
  }
  //----------------------------------------------------------------------------
  // fetch truncation error for decay (fixed point only)
  template <uint Idx>
  auto fetch_quantizer_err()
  {
    static_assert (spec::is_quantizer (Idx));
    constexpr auto offset = get_states_offset (Idx);
    return _stage[Idx].z[offset];
  }
  //----------------------------------------------------------------------------
  // store truncation error for decay (fixed point only)
  template <uint Idx>
  void save_quantizer_err (s16 err)
  {
    static_assert (spec::is_quantizer (Idx));
    constexpr auto offset = get_states_offset (Idx);
    _stage[Idx].z[offset] = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  std::array<xspan<s16>, 2> get_read_buffers (uint blocksize, uint neg_offset)
  {
    constexpr uint n_spls = spec::get_delay_spls (Idx).to_int();
    constexpr auto sz     = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= (blocksize + neg_offset));

    uint block1 = _stage[Idx].pos - n_spls - neg_offset;
    block1 += (block1 >= sz) ? sz : 0;
    uint end = _stage[Idx].pos - n_spls + blocksize - neg_offset - 1;
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
    constexpr auto sz = get_delay_size (Idx);
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
  using spec    = spec_access<Spec_array>;

  struct stage {
    s16* z {};
    uint pos {};
  };
  std::array<stage, spec::size()> _stage;
  lowbias32_hash<2>               _noise;
};
}}} // namespace artv::detail::lofiverb
//------------------------------------------------------------------------------
