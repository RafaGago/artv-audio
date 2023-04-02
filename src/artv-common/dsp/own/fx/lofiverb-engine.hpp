#pragma once

// internals of lofiverb.hpp

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>
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

using float16 = f16pack<5, -1, f16pack_dftz | f16pack_clamp>;
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

struct block_delay_data {
  fixpt_spls spls;
  fixpt_spls extra_spls; // In samples. To be able to do negative offsets
};

struct filter_data {
  fixpt_t g;
  bool    is_lowpass;
};

struct crossover_data {
  fixpt_t g; // fc
  fixpt_t g_lp; // lowpass gain
  fixpt_t g_hp; // lowpass highpass gain
};

struct quantizer_data {
  fixpt_t g;
};

using stage_data = std::variant<
  allpass_data,
  delay_data,
  block_delay_data,
  filter_data,
  crossover_data,
  quantizer_data>;
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_ap (
  u16   spls,
  float g   = 0.f,
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
// block delays have to be at least one block long. It has an additional
// parameter to overallocate memory to allow block processing "extra_spls".
//
// If the last element of a reverb loop is a delay of at least one block, all
// the outputs for a block can be fetched from existing previous work, before
// even starting processing the incoming sample batch.
//
// If those samples are to be used as feedback to sum with the current input
// batch, then the feedback samples need the last output from the past block to
// be added with the first incoming sample.
//
// So in total it is needed to fetch 1 block plus one old sample for feedback
// purposes, that's what the "extra_spls" parameter is for. This is used in
// conjunction with a value of 1 on "negative_offset" on the fetch_block
// function plus passing a buffer of the desired size + 1 spl.
static constexpr detail::lofiverb::stage_data make_block_delay (
  u16 spls,
  u16 extra_spls = 1)
{
  return block_delay_data {
    fixpt_spls::from_int (spls), fixpt_spls::from_int (extra_spls)};
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
static constexpr detail::lofiverb::stage_data make_crossover (
  float g    = 0,
  float g_lp = 0,
  float g_hp = 0)
{
  return crossover_data {
    fixpt_t::from_float (g),
    fixpt_t::from_float (g_lp),
    fixpt_t::from_float (g_hp)};
}
//------------------------------------------------------------------------------
static constexpr detail::lofiverb::stage_data make_quantizer()
{
  return quantizer_data {};
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
  static constexpr bool is_block_delay (uint i)
  {
    return std::holds_alternative<block_delay_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_delay (uint i)
  {
    return std::holds_alternative<delay_data> (values[i]) || is_block_delay (i);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_filter (uint i)
  {
    return std::holds_alternative<filter_data> (values[i]) || is_crossover (i);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_quantizer (uint i)
  {
    return std::holds_alternative<quantizer_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_lowpass (uint i)
  {
    if (is_filter (i) && !is_crossover (i)) {
      return std::get<filter_data> (values[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_highpass (uint i)
  {
    if (is_filter (i) && !is_crossover (i)) {
      return !std::get<filter_data> (values[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_crossover (uint i)
  {
    return std::holds_alternative<crossover_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_t get_gain (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (values[i]).g;
    }
    else if (is_crossover (i)) {
      return std::get<crossover_data> (values[i]).g;
    }
    else if (is_filter (i)) {
      return std::get<filter_data> (values[i]).g;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr auto get_crossover_data (uint i)
  {
    if (is_crossover (i)) {
      return std::get<crossover_data> (values[i]);
    }
    else {
      return nullptr;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls get_delay_spls (uint i)
  {
    if (is_allpass (i)) {
      return std::get<allpass_data> (values[i]).spls;
    }
    else if (is_block_delay (i)) {
      return std::get<block_delay_data> (values[i]).spls;
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
    else if (is_delay (i) && !is_block_delay (i)) {
      return std::get<delay_data> (values[i]).mod;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_extra_spls (uint i)
  {
    if (is_block_delay (i)) {
      return std::get<block_delay_data> (values[i]).extra_spls.to_int();
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
  static constexpr uint get_min_delay_spls (uint i)
  {
    return get_delay_spls (i).to_int() - get_delay_mod_spls (i).to_int();
  }
  //----------------------------------------------------------------------------
  static constexpr std::size_t size() { return values.size(); }
  //----------------------------------------------------------------------------
private:
  static constexpr auto values {Spec_array::values};
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
namespace state {
// these need to be outside of engine, as the constexpr function has to be fully
// defined before.
//------------------------------------------------------------------------------
struct empty {
  // struct clang_bug_workaround {};
  // clang_bug_workaround v;
};
//------------------------------------------------------------------------------
struct y1_float {
  float y1;
};
//------------------------------------------------------------------------------
struct y1_fixpt {
  fixpt_t y1;
  s16     y1_err;
};
//------------------------------------------------------------------------------
struct allpass_fixpt {
  s16 x_err;
  s16 u_err;
};
//------------------------------------------------------------------------------
struct allpass_and_y1_fixpt {
  fixpt_t y1;
  s16     y1_err;
  s16     x_err;
  s16     u_err;
};
//------------------------------------------------------------------------------
struct quantizer_fixpt {
  s16 err;
};
//------------------------------------------------------------------------------
template <class T, uint Idx, class SpecAccess>
static constexpr auto get_state_type()
{
  if constexpr (is_fixpt_v<T>) {
    if constexpr (SpecAccess::has_modulated_delay (Idx)) {
      if constexpr (SpecAccess::is_allpass (Idx)) {
        return allpass_and_y1_fixpt {};
      }
      else {
        return y1_fixpt {};
      }
    }
    else if constexpr (SpecAccess::is_allpass (Idx)) {
      return allpass_fixpt {};
    }
    else if constexpr (SpecAccess::is_filter (Idx)) {
      return y1_fixpt {};
    }
    else if constexpr (SpecAccess::is_quantizer (Idx)) {
      return quantizer_fixpt {};
    }
    else {
      return empty {};
    }
  }
  else {
    static_assert (std::is_same_v<T, float>);
    if constexpr (
      SpecAccess::has_modulated_delay (Idx) || SpecAccess::is_filter (Idx)) {
      return y1_float {};
    }
    else {
      return empty {};
    }
  }
}
//------------------------------------------------------------------------------
template <class T, class SpecAccess>
struct index_to_state_qfn {

  template <class Idx>
  using fn = decltype (get_state_type<T, Idx::value, SpecAccess>());
};
//------------------------------------------------------------------------------

} // namespace state

struct defaulted_tag {};
static constexpr defaulted_tag defaulted {};
//------------------------------------------------------------------------------
// A class to abstract 16-bit storage, queue access and common DSP operations
// when building reverbs based on allpass loops. Both on fixed and floating
// point.
template <class Spec_array, uint Max_block_size>
class engine {
public:
  //----------------------------------------------------------------------------
  using spec = spec_access<Spec_array>;
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
  //----------------------------------------------------------------------------
  static constexpr uint get_required_size()
  {
    uint ret = 0;
    for (uint i = 0; i < decltype (_stage) {}.size(); ++i) {
      ret += get_delay_size (i);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  constexpr void reset_memory (xspan<s16> mem)
  {
    for (uint i = 0; i < _stage.size(); ++i) {
      _stage[i].z = mem.cut_head (get_delay_size (i)).data();
    }
    memset (&_states, 0, sizeof _states);
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
    static_assert (spec::is_block_delay (Idx), "Only usable on block delays");
    static_assert (spec::get_delay_spls (Idx).to_int() >= max_block_size);
    assert (dst);
    decode_read (dst, get_read_buffers<Idx> (dst.size(), negative_offset));
  }

  // as above, but gets the contents added to "io"
  template <uint Idx, class T>
  void fetch_block_add (xspan<T> io, uint negative_offset = 0)
  {
    static_assert (spec::is_block_delay (Idx), "Only usable on block delays");
    static_assert (spec::get_delay_spls (Idx).to_int() >= max_block_size);
    assert (io);
    decode_read_add (io, get_read_buffers<Idx> (io.size(), negative_offset));
  }

  // see comment on fetch_block
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    static_assert (spec::is_block_delay (Idx), "Only usable on block delays");
    static_assert (spec::get_delay_spls (Idx).to_int() >= max_block_size);
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
#if 0
  // This doesn't optimize the datatype range, for ints it truncates more times
  // than necessary
  //----------------------------------------------------------------------------
  template <class T, uint N>
  static void hadamard_no_norm (std::array<T*, N> io, uint block_size)
  {
    static_assert (is_pow2 (N));
    if constexpr (N > 1) {
      constexpr int half_n = N / 2;
      // all this superfluous array copying stuff should be removed by the
      // optimizer.
      hadamard_no_norm (array_slice<0, half_n> (io), block_size);
      hadamard_no_norm (array_slice<half_n, half_n> (io), block_size);

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
#else
  // this might belong somewhere else
  template <class T>
  static void hadamard4 (std::array<T*, 4> io, uint block_size)
  {
    // these quantizations could use fraction saving too...
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < block_size; ++i) {
      auto y1 = io[0][i] - io[1][i];
      auto y2 = io[0][i] + io[1][i];
      auto y3 = io[2][i] - io[3][i];
      auto y4 = io[2][i] + io[3][i];

      io[0][i] = (T) ((y1 - y3) * 0.5_r);
      io[1][i] = (T) ((y2 - y4) * 0.5_r);
      io[2][i] = (T) ((y1 + y3) * 0.5_r);
      io[3][i] = (T) ((y2 + y4) * 0.5_r);
    }
  }
#endif
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class... Ts>
  void run (xspan<T> io, Ts&&... args)
  {
    if constexpr (spec::is_allpass (Idx) || spec::is_delay (Idx)) {
      constexpr uint n_args          = sizeof...(Ts);
      constexpr uint expected_n_args = 2;

      if constexpr (expected_n_args == n_args) {
        run_ap_or_delay<Idx> (io, std::forward<Ts> (args)...);
      }
      else {
        // add null on the remaining positions
        using filler_defaulted = mp11::
          mp_repeat_c<std::tuple<defaulted_tag>, expected_n_args - n_args>;
        std::apply (
          [=] (auto&&... targs) {
            run_ap_or_delay<Idx> (
              io, std::forward<decltype (targs)> (targs)...);
          },
          std::tuple_cat (
            std::forward_as_tuple (std::forward<Ts> (args)...),
            filler_defaulted {}));
      }
    }
    else if constexpr (spec::is_filter (Idx)) {
      if constexpr (spec::is_crossover (Idx)) {
        run_crossover<Idx> (io, std::forward<Ts> (args)...);
      }
      else if constexpr (spec::is_lowpass (Idx)) {
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
  // for arbitrarily nested allpasses or crossovers
  template <uint Idx1, uint Idx2, uint... Idxs, class T, class... Ts>
  void run (xspan<T> io, Ts&&... args)
  {
    if constexpr (
      spec::is_allpass (Idx1) && can_be_placed_on_nested_allpass<Idx2>
      && (... && can_be_placed_on_nested_allpass<Idxs>) ) {
      run_nested_ap<Idx1, Idx2, Idxs...> (io, std::forward<Ts> (args)...);
    }
    else if constexpr (
      sizeof...(Idxs) == 0 && spec::is_lowpass (Idx1)
      && spec::is_lowpass (Idx2)) {
      // 3-band crossover
      run_3band_crossover<Idx1, Idx2> (io, std::forward<Ts> (args)...);
    }
    else {
      static_assert (
        sizeof (T) != sizeof (T), "Invalid type on one of the indexes");
    }
  }
  //----------------------------------------------------------------------------
private:
  template <uint Idx>
  static constexpr bool can_be_placed_on_nested_allpass
    = spec::is_allpass (Idx) || spec::is_delay (Idx) || spec::is_lowpass (Idx)
    || spec::is_highpass (Idx);
  //----------------------------------------------------------------------------
  template <class T, class Lfo>
  static constexpr auto get_lfo_generator (Lfo&& v)
  {
    using Lfo_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Lfo>>;
    if constexpr (is_generator<T, Lfo_no_cv_ref>()) {
      return std::forward<Lfo_no_cv_ref> (v);
    }
    else if constexpr (is_array_subscriptable_v<Lfo_no_cv_ref>) {
      return [v] (uint i) { return v[i]; };
    }
    else {
      return defaulted_tag {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Gain>
  static constexpr auto get_gain_generator (Gain&& v)
  {
    using Gain_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Gain>>;
    if constexpr (is_generator<T, Gain_no_cv_ref>()) {
      return std::forward<Gain_no_cv_ref> (v);
    }
    else if constexpr (is_array_subscriptable_v<Gain_no_cv_ref>) {
      return [v] (uint i) { return v[i]; };
    }
    else if constexpr (std::is_same_v<T, Gain_no_cv_ref>) {
      return [v] (uint) { return v; };
    }
    else if constexpr (std::is_same_v<defaulted_tag, Gain_no_cv_ref>) {
      return [] (uint) { return (T) spec::get_gain (Idx); };
    }
    else {
      static_assert (
        sizeof (Gain_no_cv_ref) != sizeof (Gain_no_cv_ref), "unknown type");
    }
  }
  //----------------------------------------------------------------------------
  template <uint N>
  static constexpr auto get_hadamard_ratio()
  {
    // no constexpr 1/sqrt(N)
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
    if constexpr (is_fixpt_v<T>) {
      static_assert (T::n_bits == 16 && T::n_sign == 1);
      dst = T::from (src);
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
  void decode_read_add (xspan<T> dst, std::array<xspan<s16>, 2> src)
  {
    assert (dst.size() == (src[0].size() + src[1].size()));
    for (uint i = 0; i < src[0].size(); ++i) {
      T v;
      decode_read (v, src[0][i]);
      dst[i] = (T) (dst[i] + v);
    }
    dst.cut_head (src[0].size());
    for (uint i = 0; i < src[1].size(); ++i) {
      T v;
      decode_read (v, src[1][i]);
      dst[i] = (T) (dst[i] + v);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  void encode_write (s16& dst, T src)
  {
    if constexpr (is_fixpt_v<T>) {
      static_assert (T::n_bits == 16 && T::n_sign == 1);
      dst = src.value();
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
  T get (uint delay_spls)
  {
    constexpr auto sz = get_delay_size (Idx);

    assert (delay_spls <= sz);

    uint z = _stage[Idx].pos - delay_spls;
    z += (z >= sz) ? sz : 0;
    T ret;
    decode_read (ret, _stage[Idx].z[z]);
    return ret;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push_one (T v)
  {
    encode_write (_stage[Idx].z[_stage[Idx].pos], v);
    advance_pos<Idx>();
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (xspan<fixpt_t> io, Func&& fn)
  {
    // decay with error-feedback/ fraction saving
    constexpr uint mask
      = lsb_mask<uint> (fixpt_acum_t::n_bits - fixpt_t::n_bits);

    state::quantizer_fixpt& st  = std::get<Idx> (_states.fix);
    auto                    err = st.err;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto in        = fn (io[i], i);
      auto [v, err_] = round (in, err);
      io[i]          = v;
      err            = err_;
    }
    st.err = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (xspan<float> io, Func&& fn)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = fn (io[i], i);
    }
  }
  //----------------------------------------------------------------------------
  enum onepole_type { lp, hp, crossv };
  //----------------------------------------------------------------------------
  template <uint Idx, onepole_type Type, class Gain>
  void run_1pole (xspan<fixpt_t> io, Gain g, fixpt_t g_lp, fixpt_t g_hp)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);

    state::y1_fixpt& st  = std::get<Idx> (_states.fix);
    auto             y1  = st.y1;
    auto             err = st.y1_err;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto         x = (fixpt_acum_t) io[i];
      fixpt_acum_t gw;
      if constexpr (is_generator<fixpt_t, Gain>()) {
        gw = (fixpt_acum_t) g (i);
      }
      else if constexpr (is_array_subscriptable_v<Gain>) {
        gw = (fixpt_acum_t) g[i];
      }
      else {
        gw = (fixpt_acum_t) g;
      }
      auto y           = (1_r - g) * x + g * y1;
      auto [y1_, err_] = truncate (y, err);
      y1               = y1_;
      err              = err_;
      if constexpr (Type == onepole_type::lp) {
        io[i] = y1;
      }
      else if constexpr (Type == onepole_type::hp) {
        io[i] = (fixpt_t) (x - y1);
      }
      else {
        io[i] = (fixpt_t) (x * g_lp + (x - y1) * g_hp);
      }
    }
    st.y1     = y1;
    st.y1_err = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, onepole_type Type, class Gain>
  void run_1pole (xspan<float> io, Gain g, float g_lp, float g_hp)
  {
    static_assert (spec::is_filter (Idx));
    assert (io);

    state::y1_float& st = std::get<Idx> (_states.flt);
    float            y1 = st.y1;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      auto  x = io[i];
      float gv;
      if constexpr (is_generator<float, Gain>()) {
        gv = (float) g (i);
      }
      else if constexpr (is_array_subscriptable_v<Gain>) {
        gv = (float) g[i];
      }
      else {
        gv = (float) g;
      }
      y1 = (1.f - gv) * x + gv * y1;
      if constexpr (Type == onepole_type::lp) {
        io[i] = y1;
      }
      else if constexpr (Type == onepole_type::hp) {
        io[i] = (float) (x - y1);
      }
      else {
        io[i] = (float) (x * g_lp + (x - y1) * g_hp);
      }
    }
    st.y1 = y1;
  }

  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_lp (xspan<T> io, U&& g)
  {
    run_1pole<Idx, onepole_type::lp> (io, std::forward<U> (g), T {}, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (xspan<T> io)
  {
    run_lp<Idx> (io, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_hp (xspan<T> io, U&& g)
  {
    run_1pole<Idx, onepole_type::hp> (io, std::forward<U> (g), T {}, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (xspan<T> io)
  {
    run_hp<Idx> (io, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx1, uint Idx2, class T, class U>
  void run_3band_crossover (xspan<T> io, U&& f_lo, U&& g_lo, U&& f_hi, U&& g_hi)
  {
    std::array<T, max_block_size> hi, mid;

    assert (io.size() <= max_block_size);
    xspan_memdump (hi.data(), io);
    run_lp<Idx1> (io, f_hi); // io has lows + mids. highs removed
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      hi[i] = (T) (hi[i] - io[i]);
    }
    xspan_memdump (mid.data(), io);
    run_lp<Idx2> (io, f_lo); // io has lows. mids + highs removed
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < io.size(); ++i) {
      io[i] = (T) ((io[i] * g_lo) + (hi[i] * g_hi) + mid[i] - io[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_crossover (xspan<T> io, U&& g, U&& g_lp, U&& g_hp)
  {
    run_1pole<Idx, onepole_type::crossover> (io, g, g_lp, g_hp);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_crossover (xspan<T> io)
  {
    constexpr detail::lofiverb::crossover_data d
      = spec::get_crossover_data (Idx);
    run_crossover (io, (T) d.g, (T) d.g_lp, (T) d.g_hp);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<fixpt_t> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).add_sign();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).add_sign();

    using fixpt_th      = decltype (fixpt_acum_t {}.resize<2, -2>());
    constexpr uint mask = lsb_mask<uint> (fixpt_th::n_frac - fixpt_t::n_frac);

    auto& st  = std::get<Idx> (_states.fix);
    auto  y1  = st.y1;
    auto  err = st.y1_err;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < dst.size(); ++i) {
      auto fixpt_spls = (delay_spls + (lfo_gen (i) * delay_mod_spls));
      auto n_spls     = fixpt_spls.to_int();
      n_spls -= i;
      auto n_spls_frac = (fixpt_th) fixpt_spls.fractional();
      auto d           = n_spls_frac + 0.418_r; // this might exceed 1
      auto a           = (1_r - d) / (1_r + d); // 0.4104 to -1
      auto z0          = (fixpt_th) get<Idx, fixpt_t> (n_spls - 1);
      auto z1          = (fixpt_th) get<Idx, fixpt_t> (n_spls);
      auto y           = (z0 * a) + z1 - (a * y1);
      auto [y1_, err_] = truncate (y, err);
      y1               = y1_;
      err              = err_;
      dst[i]           = y1;
    }
    st.y1     = y1;
    st.y1_err = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class LF>
  void run_thiran (xspan<float> dst, LF&& lfo_gen)
  {
    static_assert (spec::has_modulated_delay (Idx));
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).to_floatp();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).to_floatp();

    auto& st = std::get<Idx> (_states.flt);
    float y1 = st.y1;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
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
    st.y1 = y1;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  auto run_allpass (T in, T yn, U g)
  {
    static_assert (spec::is_allpass (Idx));
    auto u = in + yn * g;
    auto x = yn - u * g;
    if constexpr (std::is_floating_point_v<T>) {
      return std::make_tuple (x, u);
    }
    else {
      auto& st          = std::get<Idx> (_states.fix);
      auto [q_x, x_err] = truncate (x, st.x_err);
      auto [q_u, u_err] = truncate (u, st.u_err);
      st.x_err          = x_err;
      st.u_err          = u_err;
      return std::make_tuple (q_x, q_u);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class FG>
  auto run_allpass (xspan<T> io, xspan<T> yn, FG&& g_gen)
  {
    static_assert (spec::is_allpass (Idx));
    if constexpr (std::is_floating_point_v<T>) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        auto [x, u] = run_allpass<Idx, T> (io[i], yn[i], g_gen (i));
        io[i]       = x;
        yn[i]       = u;
      }
    }
    else {
      auto& st    = std::get<Idx> (_states.fix);
      auto  x_err = st.x_err;
      auto  u_err = st.u_err;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < io.size(); ++i) {
        auto g             = g_gen (i);
        auto u             = io[i] + yn[i] * g;
        auto x             = yn[i] - u * g;
        auto [q_x, x_err_] = truncate (x, x_err);
        auto [q_u, u_err_] = truncate (u, u_err);
        x_err              = x_err_;
        u_err              = u_err_;
        io[i]              = q_x;
        yn[i]              = q_u;
      }
      st.x_err = x_err;
      st.u_err = u_err;
    }
  }
  //----------------------------------------------------------------------------
  // Template function to run:
  // - delays
  // - modulated delays
  // - non-nested allpasses
  // - modulated non-nested allpasses.
  //
  // Lfo and gain parameters can be either value generators (lambdas with
  // T(int) signature), spans or dummy types. Use "defaulted_tag" to pass
  // defaulted parameters (no delay modulation, same allpass gain as in the
  // specification)
  template <uint Idx, class T, class Lfo, class G>
  void run_ap_or_delay (xspan<T> io, Lfo&& lfo, G&& gain)
  {
    static_assert (!spec::is_filter (Idx));

    constexpr auto sz    = get_delay_size (Idx);
    constexpr auto minsz = sz - spec::get_delay_mod_spls (Idx).to_int();

    constexpr auto szw = get_delay_size (Idx);

    auto g_gen   = get_gain_generator<Idx, T> (std::forward<G> (gain));
    auto lfo_gen = get_lfo_generator<T> (std::forward<Lfo> (lfo));

    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      assert (io.size() <= minsz);

      std::array<T, max_block_size> z_mem;
      xspan                         z {z_mem.data(), io.size()};

      if constexpr (spec::has_modulated_delay (Idx)) {
        run_thiran<Idx> (z, lfo_gen);
      }
      else {
        decode_read (z, get_read_buffers<Idx> (z.size(), 0));
      }
      if constexpr (spec::is_allpass (Idx)) {
        run_allpass<Idx, T> (io, z, g_gen);
        encode_write (prepare_block_insertion<Idx> (io.size()), z.to_const());
      }
      else {
        encode_write (prepare_block_insertion<Idx> (io.size()), io.to_const());
        xspan_memcpy (io, z);
      }
    }
    else {
      // overlap, needs single-sample iteration
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
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
          io[i]            = out;
          push_one<Idx> (push);
        }
        else {
          push_one<Idx> (io[i]);
          io[i] = delayed;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <std::size_t N>
  static constexpr int get_previous_ap_pos (
    uint                       pos,
    std::array<uint, N> const& idxs)
  {
    if (pos == 0) {
      return -1;
    }
    for (uint i = pos - 1; i < pos; --i) {
      if (spec::is_allpass (idxs[i])) {
        return i;
      }
    }
    return -1;
  }
  //----------------------------------------------------------------------------
  template <std::size_t N>
  static constexpr std::array<uint, N + 1> get_arg_offsets (
    std::array<uint, N> idxs)
  {
    std::array<uint, N + 1> ret {};
    for (uint i = 0; i < idxs.size(); ++i) {
      if (spec::is_allpass (idxs[i])) {
        // lfo + gain parameter.
        ret[i + 1] = 2;
        continue;
      }
      else if (spec::is_delay (idxs[i])) {
        // lfo for modulation
        ret[i + 1] = 1;
        continue;
      }
      else if (spec::is_lowpass (idxs[i]) || spec::is_highpass (idxs[i])) {
        // gain parameter.
        ret[i + 1] = 1;
        continue;
      }
      else {
        assert (false); // Type not supported/added yet
      }
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <std::size_t N>
  static constexpr uint constexpr_accumulate (std::array<uint, N> v)
  {
    uint ret = 0;
    for (auto n : v) {
      ret += n;
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  // run arbitrarily nested allpass.
  // For each stage all its arguments have to be given, those are
  // consumed positionally.
  //
  // Not provided parameters will be defaulted. Use "defaulted_tag" instances to
  // fill with defaults until the one you want to provide.
  //----------------------------------------------------------------------------
  template <uint... Idx, class T, class... Ts>
  void run_nested_ap (xspan<T> io, Ts&&... argsp)
  {
    constexpr auto idxs        = make_array (Idx...);
    constexpr auto arg_offset  = get_arg_offsets (idxs);
    constexpr uint max_n_args  = constexpr_accumulate (arg_offset);
    constexpr int  last_ap_pos = get_previous_ap_pos (idxs.size(), idxs);

    static_assert (
      spec::is_allpass (idxs[0]), "1st stage has to be an allpass");
    static_assert (max_n_args >= sizeof...(Ts), "Excess arguments");
    using defaults = mp11::
      mp_repeat_c<std::tuple<defaulted_tag>, max_n_args - sizeof...(Ts)>;

    std::array<T, max_block_size> fwd;
    auto args = std::tuple_cat (std::forward_as_tuple (argsp...), defaults {});
    xspan_memdump (fwd.data(), io);

    mp_foreach_idx (
      mp_list<k_uint<Idx>...> {}, [&] (auto order, auto stage_idx) {
        bool constexpr is_serial = false; // not nested element // TBD

        if constexpr (spec::is_allpass (stage_idx) && !is_serial) {
          // having delays bigger than the block size allows to run each stage
          // once for all samples instead of each sample having to loop over
          // each stage.
          static_assert (
            spec::get_min_delay_spls (stage_idx) >= max_block_size,
            "Nested AP delays have to be GE than the block size");

          std::array<T, max_block_size> bwd;

          if constexpr (spec::has_modulated_delay (stage_idx)) {
            auto&& lfo_arg = std::get<arg_offset[order]> (args);
            auto   lfo_gen = get_lfo_generator<T> (lfo_arg);
            run_thiran<stage_idx> (xspan {bwd.data(), io.size()}, lfo_gen);
          }
          else {
            decode_read (
              xspan {bwd.data(), io.size()},
              get_read_buffers<stage_idx> (io.size(), 0));
          }
          auto&& gain_arg = std::get<arg_offset[order] + 1> (args);
          auto   gain_gen = get_gain_generator<stage_idx, T> (gain_arg);

          ARTV_LOOP_UNROLL_SIZE_HINT (16)
          for (uint i = 0; i < io.size(); ++i) {
            // See lattice form:
            // https://www.dsprelated.com/freebooks/pasp/Nested_Allpass_Filters.html
            auto [out, push]
              = run_allpass<stage_idx, T> (fwd[i], bwd[i], gain_gen (i));
            fwd[i] = push;
            bwd[i] = out;
          }
          if constexpr (order.value == 0) {
            // first allpass stage, backwards signal is the output
            xspan_memdump (io.data(), xspan {bwd.data(), io.size()});
          }
          else {
            // try to process the backwards signal with nested stages in
            // backwards order until another allpass is found.
            constexpr int ap_pos = get_previous_ap_pos (order.value, idxs);
            static_assert (ap_pos >= 0);
            constexpr uint n_steps = (order - ap_pos) - 1;
            using seq = add_offset_t<ap_pos, std::make_index_sequence<n_steps>>;

            // process backward signal with nested non-allpass elements
            mp11::mp_for_each<mp11::mp_from_sequence<seq>> (
              [&] (auto order_rev) {
                constexpr uint rev_idx       = idxs[order_rev];
                constexpr bool is_serial_rev = false; // TBD
                if constexpr (is_serial_rev) {
                  // skip... not belonging to the backwards path
                }
                if constexpr (
                  spec::is_highpass (rev_idx) || spec::is_lowpass (rev_idx)) {
                  auto gain = std::get<arg_offset[order_rev]> (args);
                  run<rev_idx> (xspan {bwd.data(), io.size()}, gain);
                }
                else if (spec::is_delay (rev_idx)) {
                  auto lfo_gen_rev = get_lfo_generator<T> (
                    std::get<arg_offset[order_rev]> (args));
                  run<rev_idx> (xspan {bwd.data(), io.size()}, lfo_gen_rev);
                }
                else {
                  static_assert (sizeof (order_rev), "Not implemented yet");
                }
              });
            // insert the backwards signal on the previous allpass...
            encode_write (
              prepare_block_insertion<idxs[ap_pos]> (io.size()),
              xspan {bwd.data(), io.size()}.to_const());
          }
        } // is_allpass
        else if constexpr ((order > last_ap_pos) || is_serial) {
          // elements after the last allpass can be processed on the forward
          // path so they are added to the last allpass on the final insertion
          // outside this foreach.
          //
          // If implementing serial elements inside the nestings they will have
          // to take this branch too, but instead of allpass they will have to
          // take a class like "serial<allpass> and be unwrapped"
          if constexpr (
            spec::is_highpass (stage_idx) || spec::is_lowpass (stage_idx)) {
            auto gain = std::get<arg_offset[order]> (args);
            run<stage_idx> (xspan {fwd.data(), io.size()}, gain);
          }
          else if constexpr (spec::is_delay (stage_idx)) {
            auto lfo_gen_rev
              = get_lfo_generator<T> (std::get<arg_offset[order]> (args));
            run<stage_idx> (xspan {fwd.data(), io.size()}, lfo_gen_rev);
          }
          else {
            static_assert (sizeof (stage_idx), "Not implemented yet");
          }
        }
        else {
          // skipping nested non-allpass elements, these are handled inside the.
          // allpass conditional.
        }
      });
    // the last allpass gets the fwd signal enqueued
    static_assert (last_ap_pos > 0);
    encode_write (
      prepare_block_insertion<idxs[last_ap_pos]> (io.size()),
      xspan {fwd.data(), io.size()}.to_const());
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_size (uint i)
  {
    uint v = spec::get_max_delay_spls (i);
    if (spec::is_block_delay (i)) {
      v += spec::get_delay_extra_spls (i);
    }
    return v;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_states_offset (uint i)
  {
    return get_delay_size (i);
  }
  //----------------------------------------------------------------------------
  // This might belong somewhere else if made more generic. Do when required.
  template <bool round, bool dither, class T>
  std::tuple<fixpt_t, s16> quantize (T spl, s16 err_prev)
  {
    static_assert (is_fixpt_v<T>);
    // Fraction-saving quantization
    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    // noise shaping only
    constexpr uint n_truncated   = fixpt_acum_t::n_bits - fixpt_t::n_bits;
    constexpr uint mask          = lsb_mask<uint> (n_truncated);
    constexpr bool bidirectional = round || dither;

    if constexpr (dither) {
      spl = fixpt_clamp (
        spl,
        fixpt_t::min() + 3_r * fixpt_t::epsilon(),
        fixpt_t::max() - 3_r * fixpt_t::epsilon());
    }
    else {
      spl = fixpt_clamp (
        spl,
        fixpt_t::min() + fixpt_t::epsilon(),
        fixpt_t::max() - fixpt_t::epsilon());
    }

    auto         out = (fixpt_acum_t) (spl + fixpt_acum_t::from (err_prev));
    fixpt_acum_t d_out;
    if constexpr (dither) {
      // not throughly tested...
      auto noise = tpdf_dither<n_truncated, fixpt_acum_t::scalar_type> (_noise);
      d_out      = out + fixpt_acum_t::from (noise[0]);
    }
    else {
      d_out = out;
    }
    assert_range (d_out);
    fixpt_t q_out;
    if constexpr (round) {
      q_out = (fixpt_tr) d_out;
    }
    else {
      q_out = (fixpt_t) d_out; // truncation
    }
    if constexpr (bidirectional) {
      fixpt_acum_t errv = out - q_out;
      auto         err  = (s16) (errv.value() & mask);
      return std::make_tuple ((fixpt_t) q_out, err);
    }
    else {
      auto err = (s16) (out.value() & mask);
      return std::make_tuple (q_out, err);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  std::tuple<fixpt_t, s16> round (T v, s16 err_prev)
  {
    return quantize<true, false> (v, err_prev);
  }
  //----------------------------------------------------------------------------
  template <class T>
  std::tuple<fixpt_t, s16> truncate (T v, s16 err_prev)
  {
    return quantize<false, false> (v, err_prev);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<!is_fixpt_v<T>>* = nullptr>
  void assert_range (T v)
  {
#ifndef NDEBUG
    assert (abs (v) < 1.f);
#endif
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<is_fixpt_v<T>>* = nullptr>
  void assert_range (T v)
  {
#ifndef NDEBUG
    assert_range (v.to_floatp());
#endif
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  std::array<xspan<s16>, 2> get_read_buffers (uint blocksize, uint neg_offset)
  {
    constexpr uint n_spls = spec::get_delay_spls (Idx).to_int();
    constexpr auto sz     = get_delay_size (Idx);
    assert (blocksize);
    assert (sz >= blocksize);

    uint block1 = _stage[Idx].pos - n_spls - neg_offset;
    block1 += (block1 >= sz) ? sz : 0;
    uint end = _stage[Idx].pos - n_spls + blocksize - neg_offset - 1;
    end += (end >= sz) ? sz : 0;

    uint block1sz, block2sz;
    if (block1 <= end) {
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
    if (block1 <= end) {
      // contiguous, no wraparound
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
  template <class T, class Candidate>
  static constexpr bool is_generator()
  {
    using C = std::remove_reference_t<Candidate>;
    return std::is_convertible_v<C, std::function<T (uint)>>
      && !std::is_same_v<C, std::nullptr_t>;
  }
  //----------------------------------------------------------------------------
  using indexes = mp11::mp_iota_c<spec::size()>;
  // TODO: maybe this class can be split on specializations for fixpt_t and then
  // float it might result in more code bloat?
  using states_flt_typelist
    = mp11::mp_transform_q<state::index_to_state_qfn<float, spec>, indexes>;
  using states_fixpt_typelist
    = mp11::mp_transform_q<state::index_to_state_qfn<fixpt_t, spec>, indexes>;

  using states_flt_tuple   = mp11::mp_rename<states_flt_typelist, std::tuple>;
  using states_fixpt_tuple = mp11::mp_rename<states_fixpt_typelist, std::tuple>;
  union states_union {
    states_flt_tuple   flt;
    states_fixpt_tuple fix;
  };
  //----------------------------------------------------------------------------
  struct stage {
    s16* z {};
    uint pos {};
  };
  std::array<stage, spec::size()> _stage {};
  states_union                    _states {};
  lowbias32_hash<1>               _noise {};
};
}}} // namespace artv::detail::lofiverb
//------------------------------------------------------------------------------
