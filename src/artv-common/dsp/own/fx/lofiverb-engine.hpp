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
#include "artv-common/misc/compiler.hpp"
#include "artv-common/misc/fixed_point.hpp"
#include "artv-common/misc/float.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/overaligned_allocator.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_approx_math.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/vec_util.hpp"
#include "artv-common/misc/xspan.hpp"
#include "boost/mp11/algorithm.hpp"

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

struct free_storage_data {
  uint count;
};

struct allpass_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
};

struct comb_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
};

// 1 tap delay
struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
};

struct variable_delay_data {
  fixpt_spls min_spls;
  fixpt_spls max_spls;
};

struct block_delay_data {
  fixpt_spls spls;
  fixpt_spls extra_spls; // In samples. To be able to do negative offsets
};

struct multitap_delay_data {
  std::array<fixpt_spls, 12> spls;
  uint                       count;
};

// A multitap delay summed in parallel
struct parallel_delay_data : public multitap_delay_data {
  std::array<fixpt_t, 12> g;
  uint                    n_outs;
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
  comb_data,
  delay_data,
  block_delay_data,
  multitap_delay_data,
  parallel_delay_data,
  filter_data,
  crossover_data,
  quantizer_data,
  variable_delay_data,
  free_storage_data>;
//------------------------------------------------------------------------------
static constexpr stage_data make_ap (u16 spls, float g = 0.f, u16 mod = 0)
{
  return allpass_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_comb (u16 spls, float g = 0.f, u16 mod = 0)
{
  return comb_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g)};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_delay (u16 spls, u16 mod = 0)
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
// batch, then the feedback samples need the last output from the past block
// to be added with the first incoming sample.
//
// So in total it is needed to fetch 1 block plus one old sample for feedback
// purposes, that's what the "extra_spls" parameter is for. This is used in
// conjunction with a value of 1 on "negative_offset" on the fetch_block
// function plus passing a buffer of the desired size + 1 spl.
static constexpr stage_data make_block_delay (u16 spls, u16 extra_spls = 1)
{
  return block_delay_data {
    fixpt_spls::from_int (spls), fixpt_spls::from_int (extra_spls)};
}
//------------------------------------------------------------------------------
template <class... Ts>
static constexpr stage_data make_multitap_delay (Ts&&... delay_spls)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert (nargs > 0);
  static_assert (nargs <= multitap_delay_data {}.spls.size());

  multitap_delay_data ret {};
  auto                tpl = std::forward_as_tuple (delay_spls...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    ret.spls[i] = fixpt_spls::from_int (std::get<i> (tpl));
  });
  ret.count = nargs;
  return ret;
}
//------------------------------------------------------------------------------
template <class... Ts>
static constexpr stage_data make_parallel_delay (
  uint n_outs,
  Ts&&... spl_gain_pair)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert ((nargs % 2) == 0);
  static_assert (nargs > 0);
  static_assert ((nargs / 2) <= parallel_delay_data {}.spls.size());

  parallel_delay_data ret {};
  auto                tpl = std::forward_as_tuple (spl_gain_pair...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    if constexpr ((i % 2) == 0) {
      ret.spls[i / 2] = fixpt_spls::from_int (std::get<i> (tpl));
    }
    else {
      ret.g[i / 2] = fixpt_t::from_float (std::get<i> (tpl));
    }
  });
  ret.count  = nargs / 2;
  ret.n_outs = n_outs;
  return ret;
}
//------------------------------------------------------------------------------
static constexpr stage_data make_hp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), false};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_lp (float g = 0)
{
  return filter_data {fixpt_t::from_float (g), true};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_crossover (
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
static constexpr stage_data make_quantizer()
{
  return quantizer_data {};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_free_storage (uint n_elems)
{
  return free_storage_data {n_elems};
}
//------------------------------------------------------------------------------
// regular delays user thiran for small interpolation. This uses quantization.
// TODO: lerp
static constexpr stage_data make_variable_delay (uint min_spls, uint max_spls)
{
  return variable_delay_data {
    fixpt_spls::from_int (min_spls), fixpt_spls::from_int (max_spls)};
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
  static constexpr bool is_comb (uint i)
  {
    return std::holds_alternative<comb_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_block_delay (uint i)
  {
    return std::holds_alternative<block_delay_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_1tap_delay (uint i)
  {
    return std::holds_alternative<delay_data> (values[i]) || is_block_delay (i);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_multitap_delay (uint i)
  {
    return std::holds_alternative<multitap_delay_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_parallel_delay (uint i)
  {
    return std::holds_alternative<parallel_delay_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_variable_delay (uint i)
  {
    return std::holds_alternative<variable_delay_data> (values[i]);
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
  static constexpr bool is_free_storage (uint i)
  {
    return std::holds_alternative<free_storage_data> (values[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_free_storage_count (uint i)
  {
    if (is_free_storage (i)) {
      return std::get<free_storage_data> (values[i]).count;
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  static constexpr xspan<fixpt_t const> get_gains (uint i)
  {
    if (is_parallel_delay (i)) {
      auto& v = std::get<parallel_delay_data> (values[i]);
      return {v.g.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (values[i]).g, 1};
    }
    else if (is_crossover (i)) {
      return {&std::get<crossover_data> (values[i]).g, 1};
    }
    else if (is_filter (i)) {
      return {&std::get<filter_data> (values[i]).g, 1};
    }
    else if (is_comb (i)) {
      return {&std::get<comb_data> (values[i]).g, 1};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_t get_gain (uint i)
  {
    auto gain = get_gains (i);
    if (gain.size() == 1) {
      return gain[0];
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
  static constexpr xspan<fixpt_spls const> get_delays_spls (uint i)
  {
    if (is_multitap_delay (i)) {
      auto& v = std::get<multitap_delay_data> (values[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_parallel_delay (i)) {
      auto& v = std::get<parallel_delay_data> (values[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (values[i]).spls, 1};
    }
    if (is_comb (i)) {
      return {&std::get<comb_data> (values[i]).spls, 1};
    }
    else if (is_block_delay (i)) {
      return {&std::get<block_delay_data> (values[i]).spls, 1};
    }
    else if (is_1tap_delay (i)) {
      return {&std::get<delay_data> (values[i]).spls, 1};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls get_delay_spls (uint i)
  {
    auto spls = get_delays_spls (i);
    if (spls.size() == 1) {
      return spls[0];
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
    if (is_comb (i)) {
      return std::get<comb_data> (values[i]).mod;
    }
    else if (is_1tap_delay (i) && !is_block_delay (i)) {
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
    if (is_parallel_delay (i) || is_multitap_delay (i)) {
      auto ds = get_delays_spls (i);
      return std::max_element (ds.begin(), ds.end())->to_int();
    }
    else if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (values[i]).max_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() + get_delay_mod_spls (i).to_int()
        + (has_modulated_delay (i) ? 1 : 0); // extra spl for thiran
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_min_delay_spls (uint i)
  {
    if (is_parallel_delay (i) || is_multitap_delay (i)) {
      auto ds = get_delays_spls (i);
      return std::min_element (ds.begin(), ds.end())->to_int();
    }
    if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (values[i]).min_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() - get_delay_mod_spls (i).to_int();
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_outs (uint i)
  {
    if (is_parallel_delay (i)) {
      return std::get<parallel_delay_data> (values[i]).n_outs;
    }
    else {
      return 1;
    }
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
// these need to be outside of engine, as the constexpr function has to be
// fully defined before.
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
struct quantizer_fixpt {
  s16 err;
};
//------------------------------------------------------------------------------
template <uint N>
struct quantizer_arr_fixpt {
  std::array<s16, N> err;
};
//------------------------------------------------------------------------------
struct allpass_and_y1_fixpt : public allpass_fixpt, public y1_fixpt {};
//------------------------------------------------------------------------------
struct quantizer_and_y1_fixpt : public quantizer_fixpt, public y1_fixpt {};
//------------------------------------------------------------------------------
template <class T, uint N>
struct free_storage {
  std::array<T, N> sto;
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
      else if constexpr (SpecAccess::is_comb (Idx)) {
        return quantizer_and_y1_fixpt {};
      }
      else {
        return y1_fixpt {};
      }
    }
    else if constexpr (SpecAccess::is_allpass (Idx)) {
      return allpass_fixpt {};
    }
    else if constexpr (SpecAccess::is_comb (Idx)) {
      return quantizer_fixpt {};
    }
    else if constexpr (SpecAccess::is_filter (Idx)) {
      return y1_fixpt {};
    }
    else if constexpr (SpecAccess::is_quantizer (Idx)) {
      return quantizer_fixpt {};
    }
    else if constexpr (SpecAccess::is_parallel_delay (Idx)) {
      return quantizer_arr_fixpt<SpecAccess::get_n_outs (Idx)> {};
    }
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<T, SpecAccess::get_free_storage_count (Idx)> {};
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
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<T, SpecAccess::get_free_storage_count (Idx)> {};
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

struct lofiverb_io_tag {};

struct defaulted_tag {};
static constexpr defaulted_tag defaulted {};

struct add_to_out_tag : public lofiverb_io_tag {};
static constexpr add_to_out_tag add_to {};

struct overwrite_out_tag : public lofiverb_io_tag {};
static constexpr overwrite_out_tag overwrite {};
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
  template <class T>
  static constexpr void span_add (
    xspan<T>       dst,
    xspan<T const> lhs,
    xspan<T const> rhs)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < rhs.size(); ++i) {
      dst[i] = (T) (lhs[i] + rhs[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <class T>
  static constexpr void span_add (xspan<T> lhs, xspan<T const> rhs)
  {
    span_add (lhs, lhs, rhs);
  }
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
  template <class T, uint Idx, uint... Idxs>
  auto get_gain_for_rt60 (float t_sec, uint srate)
  {
    constexpr uint vec_size    = 4;
    using V                    = vec<float, vec_size>;
    constexpr uint  n_vals     = 1 + sizeof...(Idxs);
    constexpr uint  n_vecs     = div_ceil (n_vals, vec_size);
    constexpr uint  n_rounded  = n_vecs * vec_size;
    constexpr float log2_milli = -9.965784284662087f; // log2(0.001)

    alignas (V) std::array<float, n_rounded> spls {
      {(float) spec::get_delay_spls (Idx),
       (float) spec::get_delay_spls (Idxs)...}};

    std::array<float, n_vals> ret_flt;

    for (uint i = 0; i < n_vecs; ++i) {
      auto exps = vec_load<V> (&spls[i * vec_size]);
      exps /= srate * t_sec;
      auto gains = vec_exp2_2dat (exps * vec_set<V> (log2_milli));
      if (i != (n_vecs - 1) || (n_rounded == n_vals)) {
        vec_store_unaligned (&ret_flt[i * vec_size], gains);
      }
      else {
        for (uint j = 0; j < (n_vals % vec_size); ++j) {
          ret_flt[i * vec_size + j] = gains[j];
        }
      }
    }
    if constexpr (std::is_same_v<T, float>) {
      return ret_flt;
    }
    else if constexpr (is_fixpt_v<T>) {
      return std::apply (
        [] (auto... x) { return make_array (T::from_float (x)...); }, ret_flt);
    }
    else {
      static_assert (sizeof (T) != sizeof (T), "unsupported type");
    }
  }
  //----------------------------------------------------------------------------
  // Pure delays with "size > blocksize" are placed before feedback loops to
  // enable block processing, so it is possible to fetch the future
  // feedbacks/outputs (minus the first) at once and then to push them all at
  // once. This is exactly what fetch/push accomplish.
  //
  // dst has to contain one element more at the head, where the previous
  // output will be placed. Once the feedback is applied to a current input,
  // this head feedback sample (last output) can be dropped.
  template <uint Idx, class T>
  void fetch_block (xspan<T> dst, uint negative_offset = 0)
  {
    static_assert (spec::is_block_delay (Idx), "Only usable on block delays");
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (dst);
    decode_read (dst, get_read_buffers<Idx> (dst.size(), negative_offset));
  }
  //----------------------------------------------------------------------------
  // fetches the output and feedback from an allpass. The output will
  // be able to be used directly. The feedback signal will have to be manually
  // pushed (see push).
  //
  // The intent is to be able to process further (no gain) the feedback
  // signal, e.g. filter it, before insertion.
  //
  // io[in, out = input samples, outoput samples at out
  // fb[out] signal, the feedback. Can be processed (e.g) filtered and the a
  // call to push is required
  template <uint Idx, class T, class Lfo, class G>
  void fetch_block (
    xspan<T>       out,
    xspan<T>       fb,
    xspan<T const> in,
    Lfo&&          lfo,
    G&&            gain)
  {
    static_assert (spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (out && fb && in);
    assert (fb.size() >= in.size());
    assert (out.size() >= in.size());

    auto g_gen = get_gain_generator<Idx, T> (std::forward<G> (gain));

    // no overlap, can run block-wise
    fetch_block_optmod<Idx> (
      xspan {fb.data(), in.size()}, std::forward<Lfo> (lfo));
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      if constexpr (spec::is_allpass (Idx)) {
        std::tie (out[i], fb[i])
          = run_allpass<Idx, T> (in[i], fb[i], g_gen (i));
      }
      else {
        std::tie (out[i], fb[i]) = run_comb<Idx, T> (in[i], fb[i], g_gen (i));
      }
    }
  }

  template <uint Idx, class T, class Lfo, class G>
  void fetch_block (xspan<T> io, xspan<T> fb, Lfo&& lfo, G&& gain)
  {
    static_assert (spec::is_allpass (Idx));
    fetch_block<Idx> (
      io, fb, io.to_const(), std::forward<Lfo> (lfo), std::forward<G> (gain));
  }
  //----------------------------------------------------------------------------
  // Fetch for combs, gets the delay modulated and with gain applied
  template <uint Idx, class T, class Lfo, class G>
  void fetch_block (xspan<T> fb, Lfo&& lfo, G&& gain)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (fb);
    fetch_block_optmod<Idx> (fb, std::forward<Lfo> (lfo));
    auto g_gen = get_gain_generator<Idx, T> (std::forward<G> (gain));
    run_quantizer<Idx> (fb.data(), fb, [&] (auto v, uint i) {
      return v * g_gen (i);
    });
  }
  // see comment on fetch_block overloads
  //----------------------------------------------------------------------------
  // Push for combs, gets the feedback signal (obtained from fetch_block and
  // maybe filtered) and the input, returns the sum.
  template <uint Idx, class T>
  void push (xspan<T> sum, xspan<T const> fb, xspan<T const> in)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (sum && fb && in);
    assert (fb.size() && in.size());
    span_add (sum, fb, in);
    encode_write (prepare_block_insertion<Idx> (in.size()), sum.to_const());
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void push (xspan<T const> src)
  {
    static_assert (spec::is_block_delay (Idx) || spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (src);
    encode_write (prepare_block_insertion<Idx> (src.size()), src);
  }
#if 0
  //----------------------------------------------------------------------------
  template <uint Start_idx, uint N, class T>
  void run_multichannel (std::array<T*, N> io, uint block_size)
  {
    mp11::mp_for_each<mp11::mp_iota_c<N>> ([&] (auto idx) {
      run<Start_idx + idx.value> (xspan {io[idx.value], block_size});
    });
  }

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
  template <uint Idx, class T>
  xspan<T> get_storage()
  {
    static_assert (spec::is_free_storage (Idx));
    if constexpr (std::is_same_v<T, fixpt_t>) {
      return {std::get<Idx> (_states.fix).sto};
    }
    if constexpr (std::is_same_v<T, float>) {
      return {std::get<Idx> (_states.flt).sto};
    }
    else {
      static_assert (sizeof (T) != sizeof (T), "Unimplemented");
      return {};
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class... Ts>
  void run (xspan<T> out, xspan<T const> in, Ts&&... args)
  {
    assert (out && in);
    assert (out.size() >= in.size());
    if constexpr (
      spec::is_allpass (Idx) || spec::is_comb (Idx)
      || spec::is_1tap_delay (Idx)) {
      constexpr uint n_args          = sizeof...(Ts);
      constexpr uint expected_n_args = 2;

      if constexpr (expected_n_args == n_args) {
        run_ap_comb_or_delay<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else {
        // add null on the remaining positions
        using filler_defaulted = mp11::
          mp_repeat_c<std::tuple<defaulted_tag>, expected_n_args - n_args>;
        std::apply (
          [&] (auto&&... targs) {
            run_ap_comb_or_delay<Idx> (
              out.data(), in, std::forward<decltype (targs)> (targs)...);
          },
          std::tuple_cat (
            std::forward_as_tuple (std::forward<Ts> (args)...),
            filler_defaulted {}));
      }
    }
    else if constexpr (spec::is_variable_delay (Idx)) {
      static_assert (
        sizeof...(args) == 1, "variable delays always require mod");
      run_variable_delay<Idx> (out, in, args...);
    }
    else if constexpr (spec::is_filter (Idx)) {
      if constexpr (spec::is_crossover (Idx)) {
        run_crossover<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else if constexpr (spec::is_lowpass (Idx)) {
        run_lp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else {
        run_hp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
    }
    else if constexpr (spec::is_quantizer (Idx)) {
      run_quantizer<Idx> (out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (sizeof (T) != sizeof (T), "Invalid");
    }
  }

  template <
    uint Idx,
    class T,
    class U,
    class Tag,
    class... Ts,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (xspan<T const> in, Tag t, U&& out, Ts&&... outs)
  {
    assert (in);
    if constexpr (spec::is_multitap_delay (Idx)) {
      run_multitap_delay<Idx> (
        in, t, std::forward<U> (out), std::forward<Ts> (outs)...);
    }
    else if constexpr (spec::is_parallel_delay (Idx)) {
      run_parallel_delays<Idx> (
        in, t, std::forward<U> (out), std::forward<Ts> (outs)...);
    }
  }

  template <
    uint Idx,
    class T,
    class Tag,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (xspan<T> io, Tag t)
  {
    run<Idx> (io.to_const(), t, io);
  }

  template <
    uint Idx,
    class T,
    class U,
    class... Ts,
    std::enable_if_t<!std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<U>>>* = nullptr>
  void run (xspan<T> io, U&& arg1, Ts&&... args)
  {
    run<Idx> (
      io, io.to_const(), std::forward<U> (arg1), std::forward<Ts> (args)...);
  }

  template <uint Idx, class T, class... Ts>
  void run (xspan<T> io)
  {
    run<Idx> (io, io.to_const());
  }
  //----------------------------------------------------------------------------
  // for arbitrarily nested allpasses or crossovers
  template <uint Idx1, uint Idx2, uint... Idxs, class T, class... Ts>
  void run (xspan<T> out, xspan<T const> in, Ts&&... args)
  {
    if constexpr (
      spec::is_allpass (Idx1) && can_be_placed_on_nested_allpass<Idx2>
      && (... && can_be_placed_on_nested_allpass<Idxs>) ) {
      run_nested_ap<Idx1, Idx2, Idxs...> (
        out.data(), in, std::forward<Ts> (args)...);
    }
    else if constexpr (
      sizeof...(Idxs) == 0 && spec::is_lowpass (Idx1)
      && spec::is_lowpass (Idx2)) {
      // 3-band crossover
      run_3band_crossover<Idx1, Idx2> (
        out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (
        sizeof (T) != sizeof (T), "Invalid type on one of the indexes");
    }
  }

  template <uint Idx1, uint Idx2, uint... Idxs, class T, class... Ts>
  void run (xspan<T> io, Ts&&... args)
  {
    run<Idx1, Idx2, Idxs...> (io, io.to_const(), std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
private:
  template <uint Idx>
  static constexpr bool can_be_placed_on_nested_allpass
    = spec::is_allpass (Idx) || spec::is_1tap_delay (Idx)
    || spec::is_lowpass (Idx) || spec::is_highpass (Idx);
  //----------------------------------------------------------------------------
  template <class T, class Lfo>
  static constexpr auto get_generic_generator (Lfo&& v)
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
  void run_quantizer (fixpt_t* out, xspan<fixpt_t const> in, Func&& fn)
  {
    // decay with error-feedback/ fraction saving
    state::quantizer_fixpt& st  = std::get<Idx> (_states.fix);
    auto                    err = st.err;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      auto [v, err_] = round (fn (in[i], i), err);
      out[i]         = v;
      err            = err_;
    }
    st.err = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (float* out, xspan<float const> in, Func&& fn)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      out[i] = fn (in[i], i);
    }
  }
  //----------------------------------------------------------------------------
  enum onepole_type { lp, hp, crossv };
  //----------------------------------------------------------------------------
  template <uint Idx, onepole_type Type, class Gain>
  void run_1pole (
    fixpt_t*             out,
    xspan<fixpt_t const> in,
    Gain                 g,
    fixpt_t              g_lp,
    fixpt_t              g_hp)
  {
    static_assert (spec::is_filter (Idx));
    assert (out && in);

    state::y1_fixpt& st  = std::get<Idx> (_states.fix);
    auto             y1  = st.y1;
    auto             err = st.y1_err;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      auto         x = (fixpt_acum_t) in[i];
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
        out[i] = y1;
      }
      else if constexpr (Type == onepole_type::hp) {
        out[i] = (fixpt_t) (x - y1);
      }
      else {
        out[i] = (fixpt_t) (x * g_lp + (x - y1) * g_hp);
      }
    }
    st.y1     = y1;
    st.y1_err = err;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, onepole_type Type, class Gain>
  void run_1pole (
    float*             out,
    xspan<float const> in,
    Gain               g,
    float              g_lp,
    float              g_hp)
  {
    static_assert (spec::is_filter (Idx));
    assert (in && out);

    state::y1_float& st = std::get<Idx> (_states.flt);
    float            y1 = st.y1;

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      auto  x = in[i];
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
        out[i] = y1;
      }
      else if constexpr (Type == onepole_type::hp) {
        out[i] = (float) (x - y1);
      }
      else {
        out[i] = (float) (x * g_lp + (x - y1) * g_hp);
      }
    }
    st.y1 = y1;
  }

  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_lp (T* out, xspan<T const> in, U&& g)
  {
    run_1pole<Idx, onepole_type::lp> (out, in, std::forward<U> (g), T {}, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_lp (T* out, xspan<T const> in)
  {
    run_lp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_hp (T* out, xspan<T const> in, U&& g)
  {
    run_1pole<Idx, onepole_type::hp> (out, in, std::forward<U> (g), T {}, T {});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T>
  void run_hp (T* out, xspan<T const> in)
  {
    run_hp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx1, uint Idx2, class T, class U>
  void run_3band_crossover (
    T*             out,
    xspan<T const> in,
    U              f_lo,
    U              g_lo,
    U              f_hi,
    U              g_hi)
  {
    block_array<T> hi, mid;

    assert (in.size() <= max_block_size);
    assert (in && out);
    run_lp<Idx1> (mid.data(), in, f_hi); // mid has lows + mids. highs removed
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      hi[i] = (T) (in[i] - mid[i]);
    }
    run_lp<Idx2> (out, xspan {mid.data(), in.size()}.to_const(), f_lo);
    // out has lows. mids + highs removed
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      out[i] = (T) ((out[i] * g_lo) + (hi[i] * g_hi) + mid[i] - out[i]);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_crossover (T* out, xspan<T const> in, U&& g, U&& g_lp, U&& g_hp)
  {
    run_1pole<Idx, onepole_type::crossover> (out, in, g, g_lp, g_hp);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class U>
  void run_crossover (T* out, xspan<T const> in)
  {
    constexpr detail::lofiverb::crossover_data d
      = spec::get_crossover_data (Idx);
    run_crossover (out, in, (T) d.g, (T) d.g_lp, (T) d.g_hp);
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
  template <uint Idx, class T, class U>
  auto run_comb (T in, T z, U g)
  {
    static_assert (spec::is_comb (Idx));
    auto u = in + z * g;
    if constexpr (std::is_floating_point_v<T>) {
      return std::make_tuple (z, u);
    }
    else {
      auto& st        = std::get<Idx> (_states.fix);
      auto [q_u, err] = truncate (u, st.err);
      st.err          = err;
      return std::make_tuple (z, q_u);
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Lfo>
  void fetch_block_optmod (xspan<T> dst, Lfo&& lfo)
  {
    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    assert (dst.size() <= minsz);
    if constexpr (minsz >= max_block_size) {
      if constexpr (spec::has_modulated_delay (Idx)) {
        auto lfo_gen = get_generic_generator<T> (std::forward<Lfo> (lfo));
        run_thiran<Idx> (dst, lfo_gen);
      }
      else {
        decode_read (dst, get_read_buffers<Idx> (dst.size(), 0));
      }
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
  void run_ap_comb_or_delay (T* out, xspan<T const> in, Lfo&& lfo, G&& gain)
  {
    // reminder, out and in can alias
    xspan          sout {out, in.size()};
    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      block_array<T> z_mem;
      xspan          z {z_mem.data(), in.size()};
      if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
        fetch_block<Idx> (
          sout,
          z,
          in.to_const(),
          std::forward<Lfo> (lfo),
          std::forward<G> (gain));
        encode_write (prepare_block_insertion<Idx> (in.size()), z.to_const());
      }
      else {
        fetch_block_optmod<Idx> (z, std::forward<Lfo> (lfo));
        encode_write (prepare_block_insertion<Idx> (in.size()), in.to_const());
        xspan_memcpy (sout, z);
      }
    }
    else {
      // overlap, needs single-sample iteration. (Notice that it could be
      // smarter and run in the smallest possible block size, as of now not
      // worth the complexity, as most delays are bigger than the block size)
      auto lfo_gen = get_generic_generator<T> (std::forward<Lfo> (lfo));
      auto g_gen   = get_gain_generator<Idx, T> (std::forward<G> (gain));

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        T z;
        if constexpr (spec::has_modulated_delay (Idx)) {
          // this is for completeness, modulations under the block size are
          // unlikely.
          z = in[i];
          run_thiran<Idx> (xspan {&z, 1}, [i, &lfo_gen] (uint) {
            return lfo_gen (i);
          });
        }
        else {
          z = read_next<Idx, T>();
        }
        if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
          if constexpr (spec::is_allpass (Idx)) {
            std::tie (out[i], z) = run_allpass<Idx, T> (in[i], z, g_gen (i));
          }
          else {
            std::tie (out[i], z) = run_comb<Idx, T> (in[i], z, g_gen (i));
          }
          push_one<Idx> (z);
        }
        else {
          push_one<Idx> (in[i]);
          out[i] = z;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class T, class... P>
  void run_multitap_delay (xspan<T const> in, Tag, T* out, P&&... outs)
  {
    constexpr bool out_overwrite = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_spls const> spls = spec::get_delays_spls (Idx);
    mp11::mp_for_each<mp11::mp_iota_c<spls.size()>> ([&] (auto i) {
      static_assert (spls[i].to_int() >= max_block_size);
    });
    // in and out can be aliased, but none of the outs shall be aliased with
    // in or within.
    block_array<T>                             tap0;
    std::array<T * artv_restrict, spls.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        return make_array<T * artv_restrict> (
          tap0.data(), static_cast<T*> (&args[0])...);
      },
      std::forward_as_tuple (outs...));

    assert (in.size() <= max_block_size);
    for (uint i = 0; i < tap_ptrs.size(); ++i) {
      uint spls_u = spls[i].to_int();
      if constexpr (out_overwrite) {
        auto dst = xspan {tap_ptrs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        decode_read (dst, get_read_buffers<Idx> (in.size(), 0, spls_u));
      }
      else {
        xspan dst {tap_ptrs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        block_array<T> tmp_mem;
        xspan          tmp {tmp_mem.data(), in.size()};
        decode_read (tmp, get_read_buffers<Idx> (in.size(), 0, spls_u));
        span_add (dst, tmp);
      }
    }
    encode_write (prepare_block_insertion<Idx> (in.size()), in.to_const());
    if constexpr (out_overwrite) {
      xspan_memdump (out, xspan {tap0.data(), in.size()});
    }
    else {
      span_add (out, xspan {tap0.data(), in.size()});
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Tag, class... Args>
  void run_parallel_delays (xspan<T const> in, Tag, Args&&... outs_arg)
  {
    constexpr bool out_overwrite     = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_t const> g = spec::get_gains (Idx);
    constexpr auto                 n_outs       = spec::get_n_outs (Idx);
    constexpr auto                 n_extra_buff = out_overwrite ? n_outs : 0;
    static_assert (n_outs <= g.size());
    static_assert (sizeof...(outs_arg) == n_outs);

    mp11::mp_repeat_c<std::tuple<block_array<T>>, g.size() - n_extra_buff>
      tapmem;
    std::apply (
      [&] (auto&&... args) {
        if constexpr (out_overwrite) {
          run_multitap_delay<Idx> (
            in, overwrite, &outs_arg[0]..., args.data()...);
        }
        else {
          run_multitap_delay<Idx> (in, overwrite, args.data()...);
        }
      },
      tapmem);
    std::array<T*, g.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        if constexpr (out_overwrite) {
          return std::array {&outs_arg[0]..., args.data()...};
        }
        else {
          return std::array {args.data()...};
        }
      },
      tapmem);
    std::array<T*, n_outs> outs = std::apply (
      [&] (auto&&... args) { return std::array {&args[0]...}; },
      std::forward_as_tuple (outs_arg...));

    if constexpr (is_fixpt_v<T>) {
      // multiply and accumulate at high resolution
      std::array<block_array<fixpt_acum_t>, n_outs> acum {};
      for (uint tap = 0; tap < tap_ptrs.size(); ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          acum[tap % n_outs][i] += tap_ptrs[tap][i] * (T) g[tap];
        }
      }
      // quantize
      state::quantizer_arr_fixpt<n_outs>& st  = std::get<Idx> (_states.fix);
      auto&                               err = st.err;

      for (uint out = 0; out < n_outs; ++out) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          auto [v, err_] = round (acum[out][i], err[out]);
          if constexpr (out_overwrite) {
            outs[out][i] = v;
          }
          else {
            outs[out][i] = (T) (v + outs[out][i]);
          }
          err[out] = err_;
        }
      }
    }
    else {
      for (uint tap = 0; tap < n_extra_buff; ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          // reminder: &out == &tap_ptrs[0]
          outs[tap][i] *= (T) g[tap];
        }
      }
      for (uint tap = n_extra_buff; tap < tap_ptrs.size(); ++tap) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          outs[tap % n_outs][i] += tap_ptrs[tap][i] * (T) g[tap];
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Range>
  void run_variable_delay (xspan<T> out, xspan<T const> in, Range&& gen)
  {
    constexpr uint min   = spec::get_min_delay_spls (Idx);
    constexpr uint max   = spec::get_max_delay_spls (Idx);
    constexpr uint range = max - min;
    auto           rgen  = get_generic_generator<T> (std::forward<Range> (gen));
    for (uint i = 0; i < in.size(); ++i) {
      T v;
      if constexpr (is_fixpt_v<T>) {
        auto gv = rgen (i);
        assert (gv >= T::from_int (0));
        v = get<Idx, T> ((uint) (min + range * gv.to_floatp()));
      }
      else {
        auto gv = rgen (i);
        assert (gv >= T {0.});
        v = get<Idx, T> ((uint) (min + range * gv));
      }
      push_one<Idx> (in[i]);
      out[i] = v; // they could be aliased...
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class T, class Range>
  void run_variable_delay (xspan<T> io, Range&& gen)
  {
    run_variable_delay<Idx, T> (io, io, std::forward<Range> (gen));
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
      else if (spec::is_1tap_delay (idxs[i])) {
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
  // Not provided parameters will be defaulted. Use "defaulted_tag" instances
  // to fill with defaults until the one you want to provide.
  //----------------------------------------------------------------------------
  template <uint... Idx, class T, class... Ts>
  void run_nested_ap (T* out, xspan<T const> in, Ts&&... argsp)
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

    block_array<T> fwd;
    auto args = std::tuple_cat (std::forward_as_tuple (argsp...), defaults {});
    xspan_memdump (fwd.data(), in);

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

          block_array<T> bwd;
          auto&&         lfo_arg = std::get<arg_offset[order]> (args);
          fetch_block_optmod<stage_idx> (
            xspan {bwd.data(), in.size()}, get_generic_generator<T> (lfo_arg));

          auto&& gain_arg = std::get<arg_offset[order] + 1> (args);
          auto   gain_gen = get_gain_generator<stage_idx, T> (gain_arg);

          ARTV_LOOP_UNROLL_SIZE_HINT (16)
          for (uint i = 0; i < in.size(); ++i) {
            // See lattice form:
            // https://www.dsprelated.com/freebooks/pasp/Nested_Allpass_Filters.html
            std::tie (bwd[i], fwd[i])
              = run_allpass<stage_idx, T> (fwd[i], bwd[i], gain_gen (i));
          }
          if constexpr (order.value == 0) {
            // first allpass stage: backwards signal is the output
            xspan_memdump (out, xspan {bwd.data(), in.size()});
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
                  run<rev_idx> (xspan {bwd.data(), in.size()}, gain);
                }
                else if (spec::is_1tap_delay (rev_idx)) {
                  auto lfo_gen_rev = get_generic_generator<T> (
                    std::get<arg_offset[order_rev]> (args));
                  run<rev_idx> (xspan {bwd.data(), in.size()}, lfo_gen_rev);
                }
                else {
                  static_assert (sizeof (order_rev), "Not implemented yet");
                }
              });
            // insert the backwards signal on the previous allpass...
            encode_write (
              prepare_block_insertion<idxs[ap_pos]> (in.size()),
              xspan {bwd.data(), in.size()}.to_const());
          }
        } // is_allpass
        else if constexpr ((order > last_ap_pos) || is_serial) {
          // elements after the last allpass can be processed on the forward
          // path so they are added to the last allpass on the final insertion
          // outside this foreach.
          //
          // If implementing serial elements inside the nestings they will
          // have to take this branch too, but instead of allpass they will
          // have to take a class like "serial<allpass> and be unwrapped"
          if constexpr (
            spec::is_highpass (stage_idx) || spec::is_lowpass (stage_idx)) {
            auto gain = std::get<arg_offset[order]> (args);
            run<stage_idx> (xspan {fwd.data(), in.size()}, gain);
          }
          else if constexpr (spec::is_1tap_delay (stage_idx)) {
            auto lfo_gen_rev
              = get_generic_generator<T> (std::get<arg_offset[order]> (args));
            run<stage_idx> (xspan {fwd.data(), in.size()}, lfo_gen_rev);
          }
          else {
            static_assert (sizeof (stage_idx), "Not implemented yet");
          }
        }
        else {
          // skipping nested non-allpass elements, these are handled inside
          // the. allpass conditional.
        }
      });
    // the last allpass gets the fwd signal enqueued
    static_assert (last_ap_pos > 0);
    encode_write (
      prepare_block_insertion<idxs[last_ap_pos]> (in.size()),
      xspan {fwd.data(), in.size()}.to_const());
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
  std::array<xspan<s16>, 2> get_read_buffers (
    uint blocksize,
    uint neg_offset,
    uint n_spls)
  {
    constexpr auto sz = get_delay_size (Idx);
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
  template <uint Idx>
  std::array<xspan<s16>, 2> get_read_buffers (uint blocksize, uint neg_offset)
  {
    return get_read_buffers<Idx> (
      blocksize, neg_offset, spec::get_delay_spls (Idx).to_int());
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
  template <class T>
  using block_array = std::array<T, max_block_size>;
  using indexes     = mp11::mp_iota_c<spec::size()>;
  // TODO: maybe this class can be split on specializations for fixpt_t and
  // then float it might result in more code bloat?
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
