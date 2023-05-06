#pragma once

// the engine to define all the algorithm blocks at compile time.

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
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
#include "boost/mp11/list.hpp"

namespace artv { namespace detail { namespace lofiverb {

//------------------------------------------------------------------------------
// truncating fixed point type for computation
using fixpt_tt = fixpt_s<1, 4, 27, fixpt_relaxed_frac_assign>;

using fixpt_t        = fixpt_tt;
using fixpt_sto16    = fixpt_s<1, 0, 15, 0>;
using fixpt_spls     = fixpt_m<0, 14, 0, 0>;
using fixpt_spls_mod = fixpt_m<0, 9, 0, 0>;

using float16 = f16pack<5, -1, f16pack_dftz | f16pack_clamp>;

static constexpr uint max_block_size = 32;
//------------------------------------------------------------------------------
enum class interpolation : u8 { thiran, linear };

struct free_storage_data {
  uint count;
};

struct allpass_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
  interpolation  interp;
};

struct comb_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  fixpt_t        g;
  interpolation  interp;
};

// 1 tap delay
struct delay_data {
  fixpt_spls     spls;
  fixpt_spls_mod mod; // In samples. If 0 the delay is not modulated
  interpolation  interp;
};

struct variable_delay_data {
  fixpt_spls min_spls;
  fixpt_spls max_spls;
  bool       mod_is_abs;
};

struct block_delay_data {
  fixpt_spls spls;
  fixpt_spls extra_spls; // In samples. To be able to do negative offsets
};

static constexpr uint max_multitap_n_elems = 12;

struct multitap_delay_data {
  std::array<fixpt_spls, max_multitap_n_elems> spls;
  uint                                         count;
};

struct multitap_mod_delay_data : public multitap_delay_data {
  std::array<fixpt_spls_mod, max_multitap_n_elems> mod;
  interpolation                                    interp;
};

// A multitap delay summed in parallel
struct parallel_delay_data : public multitap_delay_data {
  std::array<fixpt_t, max_multitap_n_elems> g;
  uint                                      n_outs;
};

struct parallel_mod_delay_data : public multitap_mod_delay_data {
  std::array<fixpt_t, max_multitap_n_elems> g;
  uint                                      n_outs;
};

struct filter_data {
  fixpt_t g;
  bool    is_lowpass;
};

struct crossover_data {
  static constexpr uint                max_n_bands = 2;
  std::array<fixpt_t, max_n_bands>     g; // fc
  std::array<fixpt_t, max_n_bands + 1> band_gain;
  uint                                 n_bands;
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
  multitap_mod_delay_data,
  parallel_delay_data,
  parallel_mod_delay_data,
  filter_data,
  crossover_data,
  quantizer_data,
  variable_delay_data,
  free_storage_data>;
//------------------------------------------------------------------------------
static constexpr stage_data make_ap (
  u16           spls,
  float         g      = 0.f,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return allpass_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g),
    interp};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_comb (
  u16           spls,
  float         g      = 0.f,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return comb_data {
    fixpt_spls::from_int (spls),
    fixpt_spls_mod::from_int (mod),
    fixpt_t::from_float (g),
    interp};
}
//------------------------------------------------------------------------------
static constexpr stage_data make_delay (
  u16           spls,
  u16           mod    = 0,
  interpolation interp = interpolation::thiran)
{
  return delay_data {
    fixpt_spls::from_int (spls), fixpt_spls_mod::from_int (mod), interp};
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
static constexpr stage_data make_multitap_mod_delay (
  uint n_outs,
  Ts&&... spls_mod_pair)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert ((nargs % 2) == 0);
  static_assert (nargs > 0);
  static_assert ((nargs / 2) <= parallel_delay_data {}.spls.size());

  multitap_mod_delay_data ret {};
  auto                    tpl = std::forward_as_tuple (spls_mod_pair...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    if constexpr ((i % 2) == 0) {
      ret.spls[i / 2] = fixpt_spls::from_int (std::get<i> (tpl));
    }
    else {
      ret.mod[i / 2] = fixpt_t::from_float (std::get<i> (tpl));
    }
  });
  ret.count = nargs / 2;
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
template <class... Ts>
static constexpr stage_data make_mod_parallel_delay (
  uint          n_outs,
  interpolation interp,
  Ts&&... spl_gain_mod_triplet)
{
  constexpr uint nargs = sizeof...(Ts);
  static_assert ((nargs % 3) == 0);
  static_assert (nargs > 0);
  static_assert ((nargs / 3) <= parallel_mod_delay_data {}.spls.size());

  parallel_mod_delay_data ret {};
  auto                    tpl = std::forward_as_tuple (spl_gain_mod_triplet...);
  mp11::mp_for_each<mp11::mp_iota_c<nargs>> ([&] (auto i) {
    constexpr uint rem = i % 3;
    if constexpr (rem == 0) {
      ret.spls[i / 3] = fixpt_spls::from_int (std::get<i> (tpl));
    }
    else if constexpr (rem == 1) {
      ret.g[i / 3] = fixpt_t::from_float (std::get<i> (tpl));
    }
    else {
      ret.mod[i / 3] = fixpt_t::from_int (std::get<i> (tpl));
    }
  });
  ret.count  = nargs / 3;
  ret.n_outs = n_outs;
  ret.interp = interp;
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
  float f    = 0,
  float g_lo = 0,
  float g_hi = 0)
{
  crossover_data r {};
  r.g[0]         = fixpt_t::from_float (f);
  r.band_gain[0] = fixpt_t::from_float (g_lo);
  r.band_gain[1] = fixpt_t::from_float (g_hi);
  r.n_bands      = 1;
  return r;
}
//------------------------------------------------------------------------------
static constexpr stage_data make_crossover2 (
  float f_lo  = 0,
  float f_hi  = 0,
  float g_lo  = 0,
  float g_mid = 0,
  float g_hi  = 0)
{
  crossover_data r {};
  r.g[0]         = fixpt_t::from_float (f_lo);
  r.g[1]         = fixpt_t::from_float (f_hi);
  r.band_gain[0] = fixpt_t::from_float (g_lo);
  r.band_gain[1] = fixpt_t::from_float (g_mid);
  r.band_gain[2] = fixpt_t::from_float (g_hi);
  r.n_bands      = 2;
  return r;
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
static constexpr stage_data make_variable_delay (
  uint min_spls,
  uint max_spls,
  bool abs_mod = false)
{
  return variable_delay_data {
    fixpt_spls::from_int (min_spls), fixpt_spls::from_int (max_spls), abs_mod};
}
//------------------------------------------------------------------------------
template <class Algorithm>
class spec_access {
public:
  //----------------------------------------------------------------------------
  static constexpr bool is_allpass (uint i)
  {
    return std::holds_alternative<allpass_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_comb (uint i)
  {
    return std::holds_alternative<comb_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_block_delay (uint i)
  {
    return std::holds_alternative<block_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_1tap_delay (uint i)
  {
    return std::holds_alternative<delay_data> (spec[i]) || is_block_delay (i);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_multitap_delay (uint i)
  {
    return std::holds_alternative<multitap_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_multitap_mod_delay (uint i)
  {
    return std::holds_alternative<multitap_mod_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_parallel_delay (uint i)
  {
    return std::holds_alternative<parallel_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_parallel_mod_delay (uint i)
  {
    return std::holds_alternative<parallel_mod_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_variable_delay (uint i)
  {
    return std::holds_alternative<variable_delay_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool uses_absolute_delay_modulation (uint i)
  {
    if (is_variable_delay (i)) {
      return std::get<variable_delay_data> (spec[i]).mod_is_abs;
    }
    else {
      return false;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_filter (uint i)
  {
    return std::holds_alternative<filter_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_quantizer (uint i)
  {
    return std::holds_alternative<quantizer_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_lowpass (uint i)
  {
    if (is_filter (i)) {
      return std::get<filter_data> (spec[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_highpass (uint i)
  {
    if (is_filter (i)) {
      return !std::get<filter_data> (spec[i]).is_lowpass;
    }
    return false;
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_crossover_n_bands (uint i)
  {
    if (std::holds_alternative<crossover_data> (spec[i])) {
      return std::get<crossover_data> (spec[i]).n_bands;
    }
    else {
      return 0;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr bool is_free_storage (uint i)
  {
    return std::holds_alternative<free_storage_data> (spec[i]);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_free_storage_count (uint i)
  {
    if (is_free_storage (i)) {
      return std::get<free_storage_data> (spec[i]).count;
    }
    return 0;
  }
  //----------------------------------------------------------------------------
  static constexpr crossover_data get_crossover_data (uint i)
  {
    if (get_crossover_n_bands (i)) {
      return std::get<crossover_data> (spec[i]);
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr xspan<fixpt_t const> get_gains (uint i)
  {
    if (is_parallel_delay (i) || is_parallel_mod_delay (i)) {
      auto& v = std::get<parallel_delay_data> (spec[i]);
      return {v.g.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (spec[i]).g, 1};
    }
    else if (get_crossover_n_bands (i)) {
      return {get_crossover_data (i).g};
    }
    else if (is_filter (i)) {
      return {&std::get<filter_data> (spec[i]).g, 1};
    }
    else if (is_comb (i)) {
      return {&std::get<comb_data> (spec[i]).g, 1};
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
  static constexpr xspan<fixpt_spls const> get_delays_spls (uint i)
  {
    if (is_multitap_delay (i)) {
      auto& v = std::get<multitap_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_multitap_mod_delay (i)) {
      auto& v = std::get<multitap_mod_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_parallel_delay (i)) {
      auto& v = std::get<parallel_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_parallel_mod_delay (i)) {
      auto& v = std::get<parallel_mod_delay_data> (spec[i]);
      return {v.spls.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (spec[i]).spls, 1};
    }
    if (is_comb (i)) {
      return {&std::get<comb_data> (spec[i]).spls, 1};
    }
    else if (is_block_delay (i)) {
      return {&std::get<block_delay_data> (spec[i]).spls, 1};
    }
    else if (is_1tap_delay (i)) {
      return {&std::get<delay_data> (spec[i]).spls, 1};
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
  static constexpr xspan<fixpt_spls_mod const> get_delays_mod_spls (uint i)
  {
    if (is_multitap_mod_delay (i)) {
      auto& v = std::get<multitap_mod_delay_data> (spec[i]);
      return {v.mod.data(), v.count};
    }
    else if (is_parallel_mod_delay (i)) {
      auto& v = std::get<parallel_mod_delay_data> (spec[i]);
      return {v.mod.data(), v.count};
    }
    else if (is_allpass (i)) {
      return {&std::get<allpass_data> (spec[i]).mod, 1};
    }
    else if (is_comb (i)) {
      return {&std::get<comb_data> (spec[i]).mod, 1};
    }
    else if (is_1tap_delay (i) && !is_block_delay (i)) {
      return {&std::get<delay_data> (spec[i]).mod, 1};
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr fixpt_spls_mod get_delay_mod_spls (uint i)
  {
    auto spls = get_delays_mod_spls (i);
    if (spls.size() == 1) {
      return spls[0];
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_extra_spls (uint i)
  {
    if (is_block_delay (i)) {
      return std::get<block_delay_data> (spec[i]).extra_spls.to_int();
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
      return std::get<variable_delay_data> (spec[i]).max_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() + get_delay_mod_spls (i).to_int()
        + ((has_modulated_delay (i)) ? 1 // extra tail spl for thiran/lerp
                                     : 0);
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
      return std::get<variable_delay_data> (spec[i]).min_spls.to_int();
    }
    else {
      return get_delay_spls (i).to_int() - get_delay_mod_spls (i).to_int();
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_delay_buffer_size (uint i)
  {
    return (uint) get_max_delay_spls (i) + (uint) get_delay_extra_spls (i);
  }
  //----------------------------------------------------------------------------
  static constexpr uint get_n_outs (uint i)
  {
    if (is_parallel_delay (i)) {
      return std::get<parallel_delay_data> (spec[i]).n_outs;
    }
    else {
      return 1;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr std::optional<interpolation> get_interp (uint i)
  {
    if (is_allpass (i) && has_modulated_delay (i)) {
      return std::get<allpass_data> (spec[i]).interp;
    }
    else if (is_comb (i) && has_modulated_delay (i)) {
      return std::get<comb_data> (spec[i]).interp;
    }
    else if (is_1tap_delay (i) && has_modulated_delay (i)) {
      return std::get<delay_data> (spec[i]).interp;
    }
    else if (is_parallel_mod_delay (i)) {
      return std::get<parallel_mod_delay_data> (spec[i]).interp;
    }
    else {
      return {};
    }
  }
  //----------------------------------------------------------------------------
  static constexpr std::size_t size() { return spec.size(); }
  //----------------------------------------------------------------------------
private:
  static constexpr auto spec {Algorithm::get_spec()};
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
template <uint N>
struct y1_float_arr {
  std::array<y1_float, N> arr;
};
//------------------------------------------------------------------------------
struct y1_fixpt {
  fixpt_t y1;
};
//------------------------------------------------------------------------------
template <uint N>
struct y1_fixpt_arr {
  std::array<y1_fixpt, N> arr;
};
//------------------------------------------------------------------------------
template <class T_err>
struct quantizer_fixpt {
  T_err err;
};
//------------------------------------------------------------------------------
template <class value_type, uint N>
struct free_storage {
  std::array<value_type, N> sto;
};
//------------------------------------------------------------------------------
template <class value_type, uint Idx, class SpecAccess>
static constexpr auto get_state_type()
{
  constexpr auto interp = SpecAccess::get_interp (Idx);
  if constexpr (is_fixpt_v<value_type>) {
    if constexpr (!!interp && *interp == interpolation::thiran) {
      if constexpr (
        SpecAccess::is_multitap_mod_delay (Idx)
        || SpecAccess::is_parallel_mod_delay (Idx)) {
        return y1_fixpt_arr<SpecAccess::get_delays_spls (Idx).size()> {};
      }
      else {
        return y1_fixpt {};
      }
    }
    else if constexpr (SpecAccess::is_filter (Idx)) {
      return y1_fixpt {};
    }
    else if constexpr (SpecAccess::get_crossover_n_bands (Idx)) {
      return y1_fixpt_arr<SpecAccess::get_crossover_n_bands (Idx)> {};
    }
    else if constexpr (SpecAccess::is_quantizer (Idx)) {
      using integral_type = typename value_type::value_type;
      return quantizer_fixpt<integral_type> {};
    }
    else if constexpr (SpecAccess::is_parallel_delay (Idx)) {
      return empty {};
    }
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<
        value_type,
        SpecAccess::get_free_storage_count (Idx)> {};
    }
    else {
      return empty {};
    }
  }
  else {
    static_assert (std::is_same_v<value_type, float>);
    if constexpr (!!interp && *interp == interpolation::thiran) {
      if constexpr (
        SpecAccess::is_parallel_mod_delay (Idx)
        || SpecAccess::is_multitap_mod_delay (Idx)) {
        return y1_float_arr<SpecAccess::get_delays_spls (Idx).size()> {};
      }
      else {
        return y1_float {};
      }
    }
    else if constexpr (SpecAccess::is_filter (Idx)) {
      return y1_float {};
    }
    else if constexpr (SpecAccess::get_crossover_n_bands (Idx)) {
      return y1_float_arr<SpecAccess::get_crossover_n_bands (Idx)> {};
    }
    else if constexpr (SpecAccess::is_free_storage (Idx)) {
      return free_storage<
        value_type,
        SpecAccess::get_free_storage_count (Idx)> {};
    }
    else {
      return empty {};
    }
  }
}
//------------------------------------------------------------------------------
template <class value_type, class SpecAccess>
struct index_to_state_qfn {

  template <class Idx>
  using fn = decltype (get_state_type<value_type, Idx::value, SpecAccess>());
};
//------------------------------------------------------------------------------

} // namespace state

struct quantization {
  //----------------------------------------------------------------------------
  // This might belong somewhere else if made more generic. Do when required.
  template <bool round, bool dither, class T_dst, class T>
  static std::tuple<T_dst, typename T_dst::value_type> quantize (
    T                          spl,
    typename T_dst::value_type err_prev,
    lowbias32_hash<1>*         noise_gen = nullptr)
  {
    using fixpt_src      = T;
    using fixpt_dst      = T_dst;
    using dst_value_type = typename T_dst::value_type;

    static_assert (is_fixpt_v<fixpt_dst>);
    static_assert (is_fixpt_v<fixpt_src>);
    static_assert (
      fixpt_src::n_frac > fixpt_dst::n_frac, "useless quantization");

    // Fraction-saving quantization
    // https://dsp.stackexchange.com/questions/66171/single-pole-iir-filter-fixed-point-design
    // https://dsp.stackexchange.com/questions/21792/best-implementation-of-a-real-time-fixed-point-iir-filter-with-constant-coeffic
    // noise shaping only
    constexpr uint n_truncated   = T::n_frac - fixpt_dst::n_frac;
    constexpr uint mask          = lsb_mask<uint> (n_truncated);
    constexpr bool bidirectional = round || dither;

    if constexpr (dither) {
      spl = fixpt_clamp (
        spl,
        fixpt_dst::min() + 3_r * fixpt_dst::epsilon(),
        fixpt_dst::max() - 3_r * fixpt_dst::epsilon());
    }
    else {
      spl = fixpt_clamp (
        spl,
        fixpt_dst::min() + fixpt_dst::epsilon(),
        fixpt_dst::max() - fixpt_dst::epsilon());
    }

    auto    out = (fixpt_t) (spl + fixpt_t::from (err_prev));
    fixpt_t d_out;
    if constexpr (dither) {
      // not throughly tested...
      assert (noise_gen);
      auto noise
        = tpdf_dither<n_truncated, fixpt_dst::scalar_type> (*noise_gen);
      d_out = out + fixpt_t::from (noise[0]);
    }
    else {
      d_out = out;
    }
    assert (d_out.to_floatp() < 1.f);
    fixpt_dst q_out;
    if constexpr (round) {
      q_out = (fixpt_dst) d_out.to_rounding();
    }
    else {
      q_out = (fixpt_dst) d_out; // truncation
    }
    if constexpr (bidirectional) {
      fixpt_t errv = out - q_out;
      auto    err  = (dst_value_type) (errv.value() & mask);
      // err -= (err > 0); // leakage
      // err += (err < 0); // leakage
      return std::make_tuple ((fixpt_dst) q_out, err);
    }
    else {
      auto err = (dst_value_type) (out.value() & mask);
      // err -= (err > 0); // leakage
      // err += (err < 0); // leakage
      return std::make_tuple (q_out, err);
    }
  }
  //----------------------------------------------------------------------------
  template <class T_dst, class T>
  static std::tuple<T_dst, typename T_dst::value_type> round (
    T                          v,
    typename T_dst::value_type err_prev)
  {
    return quantize<true, false, T_dst> (v, err_prev);
  }
  //----------------------------------------------------------------------------
  template <class T_dst, class T>
  static std::tuple<T_dst, typename T_dst::value_type> truncate (
    T                          v,
    typename T_dst::value_type err_prev)
  {
    return quantize<false, false, T_dst> (v, err_prev);
  }
  //----------------------------------------------------------------------------
};

namespace delay {

//------------------------------------------------------------------------------
// TODO: move to delay_line? this one provides the functions to get the segments
template <class T>
class static_buffer {
public:
  //----------------------------------------------------------------------------
  using value_type = T;
  //----------------------------------------------------------------------------
  void reset (xspan<T> mem)
  {
    xspan_memset (mem);
    _z    = mem.data();
    _size = mem.size();
    _pos  = 0;
  }
  //----------------------------------------------------------------------------
  void insert (value_type v)
  {
    _z[_pos] = v;
    advance_one();
  }
  //----------------------------------------------------------------------------
  // returns the memory for a block a insertion, updates current position
  std::array<xspan<value_type>, 2> insert_block (uint blocksize)
  {
    assert (blocksize);
    assert (_size >= blocksize);

    uint block1 = _pos;
    uint end    = block1 + blocksize - 1;
    end -= (end >= _size) ? _size : 0;
    _pos = end;
    advance_one();

    uint block1sz, block2sz;
    if (block1 <= end) {
      // contiguous, no wraparound
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = _size - block1;
      block2sz = blocksize - block1sz;
    }
    return {xspan {&_z[block1], block1sz}, xspan {&_z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
  value_type read (uint delay_spls)
  {
    assert (delay_spls <= _size);
    uint z = _pos - delay_spls;
    z += (z >= _size) ? _size : 0;
    return _z[z];
  }
  //----------------------------------------------------------------------------
  std::array<xspan<value_type>, 2> read_block (uint blocksize, uint del_spls)
  {
    assert (blocksize);
    assert (del_spls <= _size);
    assert (del_spls >= blocksize);
    assert (_size >= blocksize);

    uint block1 = _pos - del_spls;
    block1 += (block1 >= _size) ? _size : 0;
    uint end = _pos - del_spls + blocksize - 1;
    end += (end >= _size) ? _size : 0;

    uint block1sz, block2sz;
    if (block1 <= end) {
      // contiguous
      block1sz = blocksize;
      block2sz = 0;
    }
    else {
      // truncated
      block1sz = _size - block1;
      block2sz = blocksize - block1sz;
    }
    return {xspan {&_z[block1], block1sz}, xspan {&_z[0], block2sz}};
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  void advance_one()
  {
    ++_pos;
    if (_pos == _size) {
      _pos = 0;
    }
  }
  //----------------------------------------------------------------------------
  value_type* _z {};
  u32         _pos {};
  u32         _size {};
};
//------------------------------------------------------------------------------

template <class T, class C>
class static_encoded_buffer {
public:
  //----------------------------------------------------------------------------
  using value_type   = T;
  using encoder      = C;
  using storage_type = typename encoder::value_type;
  static constexpr bool has_encoding
    = !std::is_same_v<value_type, storage_type>;
  //----------------------------------------------------------------------------
  void reset (xspan<storage_type> mem) { _z.reset (mem); }
  //----------------------------------------------------------------------------
  value_type read (uint spl)
  {
    assert (spl);
    if constexpr (has_encoding) {
      return _encoder.decode (_z.read (spl));
    }
    else {
      return _z.read (spl);
    }
  }
  //----------------------------------------------------------------------------
  void push (value_type v)
  {
    if constexpr (has_encoding) {
      _z.insert (_encoder.encode (v));
    }
    else {
      _z.insert (v);
    }
  }
  //----------------------------------------------------------------------------
  void read_block (xspan<value_type> dst, uint del_spls)
  {
    std::array<xspan<storage_type>, 2> src
      = _z.read_block (dst.size(), del_spls);
    assert (dst.size() == (src[0].size() + src[1].size()));
    if constexpr (has_encoding) {
      _encoder.decode (dst.data(), src[0].data(), src[0].size());
      dst.cut_head (src[0].size());
      _encoder.decode (dst.data(), src[1].data(), src[1].size());
    }
    else {
      xspan_memdump (dst.data(), src[0]);
      dst.cut_head (src[0].size());
      xspan_memdump (dst.data(), src[1]);
    }
  }
  //----------------------------------------------------------------------------
  void push_block (xspan<value_type const> src)
  {
    std::array<xspan<storage_type>, 2> dst = _z.insert_block (src.size());
    assert (src.size() == (dst[0].size() + dst[1].size()));
    if constexpr (has_encoding) {
      _encoder.encode (dst[0].data(), src.data(), dst[0].size());
      src.cut_head (dst[0].size());
      _encoder.encode (dst[1].data(), src.data(), dst[1].size());
    }
    else {
      xspan_memcpy (dst[0], src);
      src.cut_head (dst[0].size());
      xspan_memcpy (dst[1], src);
    }
  }
  //----------------------------------------------------------------------------
private:
  static_buffer<storage_type> _z {};
  encoder                     _encoder {};
};
//------------------------------------------------------------------------------

template <class T>
struct no_enconding {
  using value_type = T;
};
//------------------------------------------------------------------------------
enum class data_type : uint { fixpt16, float16, float32 };
//------------------------------------------------------------------------------

template <class T_sto>
class fixpt_encoder {
public:
  using value_type  = T_sto;
  using extern_type = fixpt_t;
  //----------------------------------------------------------------------------
  constexpr value_type encode (extern_type v)
  {
    // companding and fraction-saving together seems to be resulting in too
    // much quality.
#if 1
    auto compand = extern_type::from_float (sqrt (abs (v.to_floatp())));
    compand      = v.value() < 0 ? -compand : compand;
#if 1
    auto [conv, err] = quantization::truncate<value_type> (compand, _err);
    _err             = err;
#else
    compand = fixpt_clamp (
      compand,
      value_type::min() + value_type::epsilon(),
      value_type::max() - value_type::epsilon());
    auto conv = (value_type) compand;
#endif
#else
    auto [conv, err] = quantization::round<value_type> (v, _err);
    _err             = err;
#endif
    return conv;
  }
  //----------------------------------------------------------------------------
  void encode (value_type* dst, extern_type const* src, uint n)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < n; ++i) {
      dst[i] = encode (src[i]);
    }
  }
  //----------------------------------------------------------------------------
  constexpr extern_type decode (value_type u)
  {
    auto r = ((extern_type) u);
#if 1
    r *= r;
    return u.value() < 0 ? -r : r;
#else
    return r;
#endif
  }
  //----------------------------------------------------------------------------
  void decode (extern_type* dst, value_type const* src, uint n)
  {
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < n; ++i) {
      dst[i] = decode (src[i]);
    }
  }
  //----------------------------------------------------------------------------
private:
  typename value_type::value_type _err {};
};

//------------------------------------------------------------------------------
template <data_type Te>
class delay_buffer;

template <>
class delay_buffer<data_type::fixpt16>
  : public static_encoded_buffer<fixpt_t, fixpt_encoder<fixpt_sto16>> {
public:
  using value_type   = fixpt_t;
  using storage_type = fixpt_sto16;
};

using float16_encoder = f16pack<5, -1, f16pack_dftz | f16pack_clamp>;

template <>
class delay_buffer<data_type::float16>
  : public static_encoded_buffer<float, float16_encoder> {
public:
  // TODO: companding and using a higher resolution type as input?
  using value_type = float;
  using storage_type =
    typename static_encoded_buffer<float, float16_encoder>::storage_type;
};

template <>
class delay_buffer<data_type::float32>
  : public static_encoded_buffer<float, no_enconding<float>> {
public:
  using value_type   = float;
  using storage_type = float;
};
//------------------------------------------------------------------------------
template <data_type Te, class SpecAccess>
struct index_to_delay_buffer_qfn {

  template <class Idx>
  using fn = std::conditional_t<
    SpecAccess::get_delay_buffer_size (Idx::value) != 0,
    delay_buffer<Te>,
    std::monostate>;
};
//------------------------------------------------------------------------------
} // namespace delay

struct lofiverb_io_tag {};

struct defaulted_tag {};
static constexpr defaulted_tag defaulted {};

struct add_to_out_tag : public lofiverb_io_tag {};
static constexpr add_to_out_tag add_to {};

struct overwrite_out_tag : public lofiverb_io_tag {};
static constexpr overwrite_out_tag overwrite {};

//------------------------------------------------------------------------------
template <uint... Idxs>
struct stage_list {}; // stage list, used for avoiding the template keyword
//------------------------------------------------------------------------------
// A class whose only purpose is to decouple the index template parameters on
// "engine" to try to reduce code-bloat and probably compile time;
class algo_engine_untemplated {
public:
  //----------------------------------------------------------------------------
  template <class T, std::size_t N>
  static auto get_gain_for_rt60 (
    std::array<float, N> t_spls,
    float                t_sec,
    uint                 srate)
  {
    constexpr uint vec_size    = 4;
    using V                    = vec<float, vec_size>;
    constexpr uint  n_vals     = N;
    constexpr uint  n_vecs     = div_ceil (n_vals, vec_size);
    constexpr uint  n_rounded  = n_vecs * vec_size;
    constexpr float log2_milli = -9.965784284662087f; // log2(0.001)

    alignas (V) std::array<float, n_rounded> spls {};
    for (uint i = 0; i < t_spls.size(); ++i) {
      spls[i] = (float) t_spls[i];
    }

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
      static_assert (N != N, "unsupported type");
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class S>
  auto run_allpass (T in, T yn, U g, S& state)
  {
    auto u = in + yn * g;
    auto x = yn - u * g;
    return std::make_tuple (x, u);
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class S>
  auto run_comb (T in, T z, U g, S& state)
  {
    auto u = in + z * g;
    return std::make_tuple (z, u);
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
  template <class Candidate>
  static constexpr bool is_uint_generator()
  {
    using C = std::remove_reference_t<Candidate>;
    return std::is_convertible_v<C, std::function<uint (uint)>>
      && !std::is_same_v<C, std::nullptr_t>;
  }
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
  template <class T, class Lfo>
  static constexpr auto get_sample_generator (Lfo&& v)
  {
    using Lfo_no_cv_ref = std::remove_cv_t<std::remove_reference_t<Lfo>>;
    if constexpr (is_generator<T, Lfo_no_cv_ref>()) {
      return std::forward<Lfo_no_cv_ref> (v);
    }
    else if constexpr (is_uint_generator<Lfo_no_cv_ref>()) {
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
  template <class T, class Gain, class U>
  static constexpr auto get_gain_generator (Gain&& v, U spec_gain)
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
      return [=] (uint) { return (T) spec_gain; };
    }
    else {
      static_assert (
        sizeof (Gain_no_cv_ref) != sizeof (Gain_no_cv_ref), "unknown type");
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class S, class Func>
  void run_quantizer (T* out, xspan<T const> in, S& state, Func&& fn)
  {
    if constexpr (is_fixpt_v<T>) {
      // decay with error-feedback/ fraction saving
      auto err = state.err;
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        auto [v, err_] = quantization::round<T> (fn (in[i], i), err);
        out[i]         = v;
        err            = err_;
      }
      state.err = err;
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        out[i] = fn (in[i], i);
      }
    }
  }
  //----------------------------------------------------------------------------
  enum onepole_type { lp, hp };
  //----------------------------------------------------------------------------
  template <onepole_type Type, class T, class Gain, class State>
  void run_1pole (T* out, xspan<T const> in, Gain g, State& st)
  {
    assert (out && in);

    auto y1 = st.y1;
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      auto x = in[i];
      T    gv;
      if constexpr (is_generator<T, Gain>()) {
        gv = (T) g (i);
      }
      else if constexpr (is_array_subscriptable_v<Gain>) {
        gv = (T) g[i];
      }
      else {
        gv = (T) g;
      }
      y1 = (1_r - gv) * x + gv * y1;
      if constexpr (Type == onepole_type::lp) {
        out[i] = y1;
      }
      else {
        static_assert (Type == onepole_type::hp);
        out[i] = (T) (x - y1);
      }
    }
    st.y1 = y1;
  }
  //----------------------------------------------------------------------------
  template <class T, std::size_t N, class S>
  void run_crossover (
    T*                   out,
    xspan<T const>       in,
    std::array<T, N>     freq,
    std::array<T, N + 1> gain,
    S&                   state)
  {
    constexpr uint                      n_bands = N;
    std::array<block_array<T>, n_bands> band_mem;

    xspan_memdump (band_mem[0].data(), in);
    std::array<T*, n_bands + 1> band;
    for (uint i = 0; i < n_bands; ++i) {
      band[i] = band_mem[i].data();
    }
    band[n_bands] = out;
    // compute bands, scale each one except the lowpass
    for (uint b = 0; b < n_bands; ++b) {
      // starting from higher to lower freq.
      uint idx = n_bands - 1 - b;
      run_1pole<onepole_type::lp, T> (
        band[b + 1], xspan {band[b], in.size()}, freq[idx], state.arr[b]);

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        // finishing the previous band by subtracting this one and gain
        // scaling
        band[b][i] = (T) (band[b][i] - band[b + 1][i]);
        band[b][i] = (T) (band[b][i] * gain[idx + 1]);
      }
    }
    // scale lowpass band and sum
    std::apply (
      [out, &band, &gain, &in] (auto&&... bands) {
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          band[n_bands][i] = (T) (band[n_bands][i] * gain[0]);
          out[i]           = (T) (bands[i] + ...);
        }
      },
      band);
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class V, class Y, class D, class LF>
  void run_thiran (
    xspan<T> dst,
    U        delay_spls,
    V        delay_mod_spls,
    D&       delay,
    Y&       state_y1,
    LF&&     lfo_gen)
  {
    if constexpr (is_fixpt_v<T>) {
      auto y1 = state_y1;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        auto dlfo        = lfo_gen (i).to_dynamic();
        auto n_spls_full = (delay_spls + (dlfo * delay_mod_spls));
        uint n_spls      = n_spls_full.to_int();
        n_spls -= i;
        auto n_spls_frac = n_spls_full.fractional();
        auto d           = n_spls_frac + 0.418_r;
        auto a           = ((1_r - d) / (1_r + d)); // 0.4104 to -0.172
        auto z0          = delay.read (n_spls - 1);
        auto z1          = delay.read (n_spls);
        y1               = z0 * a + z1 - a * y1;
        dst[i]           = y1;
      }
      state_y1 = y1;
    }
    else {
      T y1 = state_y1;

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        T n_spls_full
          = delay_spls.to_floatp() + (lfo_gen (i) * delay_mod_spls.to_floatp());
        uint n_spls      = (uint) n_spls_full;
        T    n_spls_frac = n_spls_full - n_spls;
        n_spls -= i;

        T d = n_spls_frac + 0.418_r;
        T a = (1_r - d) / (1_r + d); // 0.4104 to -1

        auto z0 = delay.read (n_spls - 1);
        auto z1 = delay.read (n_spls);
        y1      = z0 * a + z1 - a * y1;
        dst[i]  = y1;
      }
      state_y1 = y1;
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class U, class V, class D, class LF>
  static void run_lerp (
    xspan<T> dst,
    U        delay_spls,
    V        delay_mod_spls,
    D&       delay,
    LF&&     lfo_gen)
  {
    static uint prev = 0;
    if constexpr (is_fixpt_v<T>) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        auto dlfo        = lfo_gen (i).to_dynamic();
        auto n_spls_full = (delay_spls + (dlfo * delay_mod_spls));
        auto n_spls      = n_spls_full.to_int();
        n_spls -= i;
        auto n_spls_frac = n_spls_full.fractional();
        auto zb          = delay.read (n_spls + 1);
        auto z           = delay.read (n_spls);
        static_assert (n_spls_frac.n_int == 0);
        dst[i] = (fixpt_t) (z * (1_r - n_spls_frac) + zb * n_spls_frac);
      }
    }
    else {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < dst.size(); ++i) {
        T n_spls_full
          = delay_spls.to_floatp() + (lfo_gen (i) * delay_mod_spls.to_floatp());
        uint  n_spls      = (uint) n_spls_full;
        float n_spls_frac = n_spls_full - n_spls;
        n_spls -= i;
        auto zb = delay.read (n_spls + 1);
        auto z  = delay.read (n_spls);
        dst[i]  = z * (1.f - n_spls_frac) + zb * n_spls_frac;
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class T, class D>
  static void run_multitap_delay (
    xspan<T * artv_restrict> outs,
    xspan<T const>           in,
    xspan<fixpt_spls const>  spls,
    bool                     out_overwrite,
    D&                       delay)
  {
    assert (spls.size() == outs.size());
    assert (in.size() <= max_block_size);

    for (uint i = 0; i < outs.size(); ++i) {
      uint spls_u = spls[i].to_int();
      if (out_overwrite) {
        auto dst = xspan {outs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        delay.read_block (xspan {dst.data(), in.size()}, spls_u);
      }
      else {
        xspan dst {outs[i], in.size()};
        assert (dst.data() != in.data()); // unwanted aliasing with input
        block_array<T> tmp_mem;
        xspan          tmp {tmp_mem.data(), in.size()};
        delay.read_block (xspan {tmp.data(), in.size()}, spls_u);
        span_add<T> (dst, tmp);
      }
    }
    delay.push_block (in);
  }
  //----------------------------------------------------------------------------
  template <class T, class D, class... Lfo>
  void run_thiran_mod_multitap_delay (
    xspan<T * artv_restrict>    outs,
    xspan<T const>              in,
    xspan<fixpt_spls const>     spls,
    xspan<fixpt_spls_mod const> mod_spls,
    bool                        out_overwrite,
    D&                          delay,
    xspan<T* const>             mod,
    xspan<T>                    y1)
  {
    assert (spls.size() == outs.size());
    assert (in.size() <= max_block_size);

    for (uint i = 0; i < outs.size(); ++i) {
      xspan dst {outs[i], in.size()};
      assert (dst.data() != in.data()); // unwanted aliasing with input
      block_array<T> tmp_mem;
      xspan out {out_overwrite ? dst.data() : tmp_mem.data(), in.size()};

      if (mod[i]) {
        run_thiran (out, spls[i], mod_spls[i], delay, y1[i], [=] (uint j) {
          return mod[i][j];
        });
      }
      else {
        run_thiran (out, spls[i], mod_spls[i], delay, y1[i], [] (uint) {
          return T {};
        });
      }
      if (!out_overwrite) {
        span_add<T> (dst, out);
      }
    }
    delay.push_block (in);
  }
  //----------------------------------------------------------------------------
  template <class T, class D, class... Lfo>
  static void run_lerp_mod_multitap_delay (
    xspan<T * artv_restrict>    outs,
    xspan<T const>              in,
    xspan<fixpt_spls const>     spls,
    xspan<fixpt_spls_mod const> mod_spls,
    bool                        out_overwrite,
    D&                          delay,
    xspan<T* const>             mod)
  {
    assert (spls.size() == outs.size());
    assert (in.size() <= max_block_size);

    for (uint i = 0; i < outs.size(); ++i) {
      xspan dst {outs[i], in.size()};
      assert (dst.data() != in.data()); // unwanted aliasing with input
      block_array<T> tmp_mem;
      xspan out {out_overwrite ? dst.data() : tmp_mem.data(), in.size()};
      if (mod[i]) {
        run_lerp (out, spls[i], mod_spls[i], delay, [=] (uint j) {
          return mod[i][j];
        });
      }
      else {
        run_lerp (out, spls[i], mod_spls[i], delay, [] (uint) { return T {}; });
      }
      if (!out_overwrite) {
        span_add<T> (dst, out);
      }
    }
    delay.push_block (in);
  }
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
    span_add<T> (lhs, lhs, rhs);
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
  template <class T>
  using block_array = std::array<T, max_block_size>;
  //----------------------------------------------------------------------------
};

// A class to abstract 16-bit storage, queue access and common DSP operations
// when building reverbs based on allpass loops. Both on fixed and floating
// point.
template <class Algorithm, delay::data_type Data_type, uint Max_block_size>
class algo_engine : private algo_engine_untemplated {
  //----------------------------------------------------------------------------
public:
  //----------------------------------------------------------------------------
  using value_type   = typename delay::delay_buffer<Data_type>::value_type;
  using storage_type = typename delay::delay_buffer<Data_type>::storage_type;
  using algorithm    = Algorithm;
  using spec         = spec_access<Algorithm>;
  //----------------------------------------------------------------------------
  static constexpr uint max_block_size = Max_block_size;
  //----------------------------------------------------------------------------
  static constexpr uint get_required_bytes()
  {
    uint elems = 0;
    mp11::mp_for_each<mp11::mp_iota_c<spec::size()>> ([&] (auto i) {
      elems += spec::get_delay_buffer_size (i);
    });
    return elems * sizeof (storage_type);
  }
  //----------------------------------------------------------------------------
  void reset_memory (xspan<u8> mem) // needs malloc/new alignment
  {
    assert (mem.size() >= get_required_bytes());
    // xspan_memset (mem);
    auto elems = mem.cast<storage_type>();
    mp11::mp_for_each<mp11::mp_iota_c<spec::size()>> ([&] (auto i) {
      constexpr uint n_elems = spec::get_delay_buffer_size (i);
      if constexpr (n_elems > 0) {
        std::get<i.value> (_stages).delay.reset (elems.cut_head (n_elems));
      }
    });
  }
  //----------------------------------------------------------------------------
  template <uint Idx, uint... Idxs>
  auto get_gain_for_rt60 (stage_list<Idx, Idxs...>, float t_sec, uint srate)
  {
    return base::get_gain_for_rt60<value_type> (
      make_array (
        (float) spec::get_delay_spls (Idx),
        (float) spec::get_delay_spls (Idxs)...),
      t_sec,
      srate);
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
  //
  // Note: Regular delays and allpass delays have gained some of the abilities
  // of block delays, the difference is that block delays can have overcapacity.
  template <uint Idx, std::enable_if_t<spec::is_block_delay (Idx)>* = nullptr>
  void fetch_block (
    stage_list<Idx>,
    xspan<value_type> dst,
    uint              negative_offset)
  {
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    std::get<Idx> (_stages).delay.read_block (
      dst, spec::get_delay_spls (Idx).to_int() + negative_offset);
  }
  //----------------------------------------------------------------------------
  // Fetch samples from the queue of a processed allpass or delay. The allpass
  // or delay has to be independenly processed.
  template <
    uint Idx,
    std::enable_if_t<
      spec::is_allpass (Idx) || spec::is_1tap_delay (Idx)>* = nullptr>
  void fetch (stage_list<Idx>, xspan<value_type> dst, uint spls)
  {
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (dst);
    assert (spls <= spec::get_max_delay_spls (Idx));
    assert (spls >= max_block_size);
    std::get<Idx> (_stages).delay.read_block (dst, spls);
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
  template <uint Idx, class Lfo, class G>
  void fetch (
    stage_list<Idx>,
    xspan<value_type>       out,
    xspan<value_type>       fb,
    xspan<value_type const> in,
    Lfo&&                   lfo,
    G&&                     gain)
  {
    static_assert (spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (out && fb && in);
    assert (fb.size() >= in.size());
    assert (out.size() >= in.size());

    auto g_gen = base::get_gain_generator<value_type> (
      std::forward<G> (gain), spec::get_gain (Idx));

    // no overlap, can run block-wise
    fetch_mod_en<Idx> (xspan {fb.data(), in.size()}, std::forward<Lfo> (lfo));
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      if constexpr (spec::is_allpass (Idx)) {
        std::tie (out[i], fb[i]) = base::run_allpass (
          in[i], fb[i], g_gen (i), std::get<Idx> (_stages).state);
      }
      else {
        std::tie (out[i], fb[i]) = base::run_comb (
          in[i], fb[i], g_gen (i), std::get<Idx> (_stages).state);
      }
    }
  }

  template <uint Idx, class Lfo, class G>
  void fetch (
    stage_list<Idx>,
    xspan<value_type> io,
    xspan<value_type> fb,
    Lfo&&             lfo,
    G&&               gain)
  {
    static_assert (spec::is_allpass (Idx));
    fetch<Idx> (
      io, fb, io.to_const(), std::forward<Lfo> (lfo), std::forward<G> (gain));
  }
  //----------------------------------------------------------------------------
  // Fetch for combs, gets the delay modulated and with gain applied
  template <uint Idx, class value_type, class Lfo, class G>
  void fetch (stage_list<Idx>, xspan<value_type> fb, Lfo&& lfo, G&& gain)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (fb);
    fetch_mod_en<Idx> (fb, std::forward<Lfo> (lfo));
    auto g_gen = base::get_gain_generator<value_type> (
      std::forward<G> (gain), spec::get_gain (Idx));
    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < fb.size(); ++i) {
      fb[i] *= g_gen (i);
    }
  }
  // see comment on fetch/fetch_block overloads
  //----------------------------------------------------------------------------
  // Push for combs, gets the feedback signal (obtained from fetch/fetch_block
  // and maybe filtered) and the input, returns the sum.
  template <uint Idx>
  void push (
    stage_list<Idx>,
    xspan<value_type>       sum,
    xspan<value_type const> fb,
    xspan<value_type const> in)
  {
    static_assert (spec::is_comb (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    assert (sum && fb && in);
    assert (fb.size() && in.size());
    base::span_add (sum, fb, in);
    std::get<Idx> (_stages).delay.push_block (xspan {sum.data(), in.size()});
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class value_type>
  void push (stage_list<Idx>, xspan<value_type const> src)
  {
    static_assert (spec::is_block_delay (Idx) || spec::is_allpass (Idx));
    static_assert (spec::get_min_delay_spls (Idx) >= max_block_size);
    std::get<Idx> (_stages).delay.push_block (src);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  xspan<value_type> get_storage (stage_list<Idx>)
  {
    static_assert (spec::is_free_storage (Idx));
    return {std::get<Idx> (_stages).state.sto};
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class... Ts>
  void run (
    stage_list<Idx>,
    xspan<value_type>       out,
    xspan<value_type const> in,
    Ts&&... args)
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
      run_variable_delay<Idx> (out, in, std::forward<Ts> (args)...);
    }
    else if constexpr (spec::is_filter (Idx)) {
      if constexpr (spec::is_lowpass (Idx)) {
        run_lp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
      else {
        run_hp<Idx> (out.data(), in, std::forward<Ts> (args)...);
      }
    }
    else if constexpr (spec::get_crossover_n_bands (Idx)) {
      run_crossover<Idx> (out.data(), in, std::forward<Ts> (args)...);
    }
    else if constexpr (spec::is_quantizer (Idx)) {
      run_quantizer<Idx> (out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (Idx != Idx, "Invalid");
    }
  }
  //----------------------------------------------------------------------------
  template <
    uint Idx,
    class U,
    class Tag,
    class... Ts,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (
    stage_list<Idx>,
    xspan<value_type const> in,
    Tag                     t,
    U&&                     out,
    Ts&&... outs)
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
    class Tag,
    std::enable_if_t<std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<Tag>>>* = nullptr>
  void run (stage_list<Idx> l, xspan<value_type> io, Tag t)
  {
    run (l, io.to_const(), t, io);
  }

  template <
    uint Idx,
    class U,
    class... Ts,
    std::enable_if_t<!std::is_base_of_v<
      lofiverb_io_tag,
      std::remove_reference_t<U>>>* = nullptr>
  void run (stage_list<Idx> l, xspan<value_type> io, U&& arg1, Ts&&... args)
  {
    run (
      l, io, io.to_const(), std::forward<U> (arg1), std::forward<Ts> (args)...);
  }

  template <uint Idx>
  void run (stage_list<Idx> l, xspan<value_type> io)
  {
    run (l, io, io.to_const());
  }
  //----------------------------------------------------------------------------
  // for arbitrarily nested allpasses or crossovers
  template <uint Idx1, uint Idx2, uint... Idxs, class... Ts>
  void run (
    stage_list<Idx1, Idx2, Idxs...>,
    xspan<value_type>       out,
    xspan<value_type const> in,
    Ts&&... args)
  {
    if constexpr (
      spec::is_allpass (Idx1) && can_be_placed_on_nested_allpass<Idx2>
      && (... && can_be_placed_on_nested_allpass<Idxs>) ) {
      run_nested_ap<Idx1, Idx2, Idxs...> (
        out.data(), in, std::forward<Ts> (args)...);
    }
    else {
      static_assert (Idx1 != Idx1, "Invalid type on one of the indexes");
    }
  }

  template <uint Idx1, uint Idx2, uint... Idxs, class... Ts>
  void run (
    stage_list<Idx1, Idx2, Idxs...> l,
    xspan<value_type>               io,
    Ts&&... args)
  {
    run (l, io, io.to_const(), std::forward<Ts> (args)...);
  }
  //----------------------------------------------------------------------------
private:
  template <uint Idx>
  static constexpr bool can_be_placed_on_nested_allpass = spec::is_allpass (Idx)
    || spec::is_1tap_delay (Idx) || spec::is_lowpass (Idx)
    || spec::is_highpass (Idx) || spec::get_crossover_n_bands (Idx) != 0;
  //----------------------------------------------------------------------------
  template <uint Idx>
  value_type get (uint delay_spls)
  {
    return std::get<Idx> (_stages).delay.read (delay_spls);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  value_type read_next()
  {
    constexpr auto spls = spec::get_delay_spls (Idx).to_int();
    return get<Idx> (spls);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Func>
  void run_quantizer (value_type* out, xspan<value_type const> in, Func&& fn)
  {
    base::run_quantizer (
      out, in, std::get<Idx> (_stages).state, std::forward<Func> (fn));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_lp (value_type* out, xspan<value_type const> in, U&& g)
  {
    base::run_1pole<base::onepole_type::lp> (
      out, in, std::forward<U> (g), std::get<Idx> (_stages).state);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_lp (value_type* out, xspan<value_type const> in)
  {
    run_lp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_hp (value_type* out, xspan<value_type const> in, U&& g)
  {
    base::run_1pole<base::onepole_type::hp> (
      out, in, std::forward<U> (g), std::get<Idx> (_stages).state);
  }
  //----------------------------------------------------------------------------
  template <uint Idx>
  void run_hp (value_type* out, xspan<value_type const> in)
  {
    run_hp<Idx> (out, in, spec::get_gain (Idx));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class U>
  void run_crossover (value_type* out, xspan<value_type const> in)
  {
    static constexpr auto                 d = spec::get_crossover_data (Idx);
    std::array<value_type, d.n_bands>     g;
    std::array<value_type, d.n_bands + 1> bg;
    for (uint i = 0; i < d.n_bands; ++i) {
      g[i]  = d.g[i];
      bg[i] = d.bg[i];
    }
    bg[d.n_bands] = d.bg[d.n_bands];
    base::run_crossover (out, in, g, bg, std::get<Idx> (_stages).state);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class P1, class P2, class P3>
  void run_crossover (
    value_type*             out,
    xspan<value_type const> in,
    P1&&                    f_lo,
    P2&&                    g_lo,
    P3&&                    g_hi)
  {
    using vt = value_type;
    base::run_crossover (
      out,
      in,
      make_array (arith_cast<vt> (f_lo)),
      make_array (arith_cast<vt> (g_lo), arith_cast<vt> (g_hi)),
      std::get<Idx> (_stages).state);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class P1, class P2, class P3, class P4, class P5>
  void run_crossover (
    value_type*             out,
    xspan<value_type const> in,
    P1&&                    f_lo,
    P2&&                    g_lo,
    P3&&                    f_mid,
    P4&&                    g_mid,
    P5&&                    g_hi)
  {
    using vt = value_type;
    base::run_crossover (
      out,
      in,
      make_array (arith_cast<vt> (f_lo), arith_cast<vt> (f_mid)),
      make_array (
        arith_cast<vt> (g_lo), arith_cast<vt> (g_mid), arith_cast<vt> (g_hi)),
      std::get<Idx> (_stages).state);
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Lfo>
  void fetch_mod_en (xspan<value_type> dst, Lfo&& lfo)
  {
    constexpr auto delay_spls     = spec::get_delay_spls (Idx).add_sign();
    constexpr auto delay_mod_spls = spec::get_delay_mod_spls (Idx).add_sign();

    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    assert (dst.size() <= minsz);
    if constexpr (minsz >= max_block_size) {
      constexpr std::optional<interpolation> interp = spec::get_interp (Idx);
      if constexpr (!interp) {
        std::get<Idx> (_stages).delay.read_block (
          dst, (uint) spec::get_delay_spls (Idx));
      }
      else {
        auto lfo_gen
          = base::get_generic_generator<value_type> (std::forward<Lfo> (lfo));
        auto& stg = std::get<Idx> (_stages);
        if constexpr (*interp == interpolation::linear) {
          base::run_lerp (dst, delay_spls, delay_mod_spls, stg.delay, lfo_gen);
        }
        else {
          static_assert (*interp == interpolation::thiran);
          base::run_thiran (
            dst, delay_spls, delay_mod_spls, stg.delay, stg.state.y1, lfo_gen);
        }
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
  // value_type(int) signature), spans or dummy types. Use "defaulted_tag" to
  // pass defaulted parameters (no delay modulation, same allpass gain as in the
  // specification)
  template <uint Idx, class Lfo, class G>
  void run_ap_comb_or_delay (
    value_type*             out,
    xspan<value_type const> in,
    Lfo&&                   lfo,
    G&&                     gain)
  {
    // reminder, out and in can alias
    xspan          sout {out, in.size()};
    constexpr auto minsz = spec::get_min_delay_spls (Idx);
    if constexpr (minsz >= max_block_size) {
      // no overlap, can run block-wise
      block_array z_mem;
      xspan       z {z_mem.data(), in.size()};
      if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
        fetch (
          stage_list<Idx> {},
          sout,
          z,
          in.to_const(),
          std::forward<Lfo> (lfo),
          std::forward<G> (gain));
        std::get<Idx> (_stages).delay.push_block (xspan {z.data(), in.size()});
      }
      else {
        fetch_mod_en<Idx> (z, std::forward<Lfo> (lfo));
        std::get<Idx> (_stages).delay.push_block (in);
        xspan_memcpy (sout, z);
      }
    }
    else {
      // overlap, needs single-sample iteration. (Notice that it could be
      // smarter and run in the smallest possible block size, as of now not
      // worth the complexity, as most delays are bigger than the block size)
      auto lfo_gen
        = base::get_generic_generator<value_type> (std::forward<Lfo> (lfo));
      auto g_gen = base::get_gain_generator<value_type> (
        std::forward<G> (gain), spec::get_gain (Idx));
      constexpr std::optional<interpolation> interp = spec::get_interp (Idx);

      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        value_type z;
        if constexpr (!interp) {
          z = read_next<Idx>();
        }
        else {
          // this is for completeness, modulations under the block size are
          // unlikely.
          z                         = in[i];
          constexpr auto delay_spls = spec::get_delay_spls (Idx).add_sign();
          constexpr auto delay_mod_spls
            = spec::get_delay_mod_spls (Idx).add_sign();
          auto& stg = std::get<Idx> (_stages);
          if constexpr (*interp == interpolation::linear) {
            base::run_lerp (
              xspan {&z, 1},
              delay_spls,
              delay_mod_spls,
              stg.delay,
              [i, &lfo_gen] (uint) { return lfo_gen (i); });
          }
          else {
            static_assert (*interp == interpolation::thiran);
            base::run_thiran (
              xspan {&z, 1},
              delay_spls,
              delay_mod_spls,
              stg.delay,
              stg.state.y1,
              [i, &lfo_gen] (uint) { return lfo_gen (i); });
          }
        }
        if constexpr (spec::is_allpass (Idx) || spec::is_comb (Idx)) {
          if constexpr (spec::is_allpass (Idx)) {
            std::tie (out[i], z) = base::run_allpass<value_type> (
              in[i], z, g_gen (i), std::get<Idx> (_stages).state);
          }
          else {
            std::tie (out[i], z) = base::run_comb<value_type> (
              in[i], z, g_gen (i), std::get<Idx> (_stages).state);
          }
          std::get<Idx> (_stages).delay.push (z);
        }
        else {
          std::get<Idx> (_stages).delay.push (in[i]);
          out[i] = z;
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  template <class Tpl>
  static auto to_ptr_array (Tpl&& tpl)
  {
    return std::apply (
      [&] (auto&&... args) { return std::array {&args[0]...}; },
      std::forward<Tpl> (tpl));
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class... P>
  void run_multitap_delay (
    xspan<value_type const> in,
    Tag,
    value_type* out,
    P&&... outs)
  {
    constexpr bool out_overwrite = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_spls const> spls = spec::get_delays_spls (Idx);
    mp11::mp_for_each<mp11::mp_iota_c<spls.size()>> ([&] (auto i) {
      static_assert (spls[i].to_int() >= max_block_size);
    });
    // in and out can be aliased, but none of the outs shall be aliased with
    // in or within.
    block_array                                         tap0;
    std::array<value_type * artv_restrict, spls.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        return make_array<value_type * artv_restrict> (
          tap0.data(), static_cast<value_type*> (&args[0])...);
      },
      std::forward_as_tuple (outs...));

    base::run_multitap_delay<value_type> (
      tap_ptrs, in, spls, out_overwrite, std::get<Idx> (_stages).delay);

    if constexpr (out_overwrite) {
      xspan_memdump (out, xspan {tap0.data(), in.size()});
    }
    else {
      base::span_add (out, xspan {tap0.data(), in.size()});
    }
  }
  //----------------------------------------------------------------------------
  template <class Type_list>
  static constexpr uint n_requiring_block_buffer()
  {
    uint n = 0;
    mp11::mp_for_each<mp11::mp_iota<mp11::mp_size<Type_list>>> ([&] (auto i) {
      using T = std::remove_cv_t<
        std::remove_reference_t<mp11::mp_at_c<Type_list, i>>>;
      n += !is_array_subscriptable_v<T> && !std::is_same_v<T, defaulted_tag>
        && !std::is_same_v<T, std::nullptr_t>;
    });
    return n;
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class... Ts>
  void run_multitap_mod_delay (
    xspan<value_type const> in,
    Tag,
    value_type* out,
    Ts&&... out_mod_pairs)
  {
    constexpr bool out_overwrite = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_spls const>     spls = spec::get_delays_spls (Idx);
    constexpr xspan<fixpt_spls_mod const> mod_spls
      = spec::get_delays_mod_spls (Idx);
    mp11::mp_for_each<mp11::mp_iota_c<spls.size()>> ([&] (auto i) {
      static_assert (spls[i].to_int() >= max_block_size);
    });

    static_assert (sizeof...(Ts) == (spls.size() * 2));

    using even
      = index_seq_mul_t<2, std::make_index_sequence<sizeof...(Ts) / 2>>;
    using odd = index_seq_add_t<1, even>;

    // in and out can be aliased, but none of the outs shall be aliased with
    // in or within.
    block_array                                         tap0;
    std::array<value_type * artv_restrict, spls.size()> tap_ptrs = std::apply (
      [&] (auto&&... args) {
        return make_array<value_type * artv_restrict> (
          tap0.data(), static_cast<value_type*> (&args[0])...);
      },
      forward_range_as_tuple (even {}, std::forward<Ts> (out_mod_pairs)...));

    // flatten the generators to buffers all at once (if required)
    auto mod_tpl
      = forward_range_as_tuple (odd {}, std::forward<Ts> (out_mod_pairs)...);
    constexpr uint n_blockbuffs = n_requiring_block_buffer<decltype (mod_tpl)>;
    std::array<block_array, n_blockbuffs>      mod_mem;
    std::array<value_type const*, spls.size()> mod;
    uint                                       mod_mem_idx = 0;

    mp11::mp_for_each<mp11::mp_iota_c<mod.size()>> ([&] (auto i) {
      auto p  = std::get<i> (mod_tpl);
      using T = std::remove_cv_t<std::remove_reference_t<decltype (p)>>;
      if constexpr (is_array_subscriptable_v<T>) {
        mod[i] = &p[0];
      }
      else if constexpr (
        std::is_same_v<T, defaulted_tag> || std::is_same_v<T, std::nullptr_t>) {
        mod[i] = nullptr;
      }
      else {
        static_assert (base::is_generator<value_type, T>);
        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint j = 0; j < in.size(); ++j) {
          mod_mem[mod_mem_idx][j] = p (j);
        }
        mod[i] = mod_mem[mod_mem_idx].data();
        ++mod_mem_idx;
      }
    });

    auto& stage = std::get<Idx> (_stages);
    if constexpr (*spec::get_interp (Idx) == interpolation::thiran) {
      base::run_thiran_mod_multitap_delay<value_type> (
        tap_ptrs,
        in,
        spls,
        mod_spls,
        out_overwrite,
        stage.delay,
        xspan {mod},
        stage.state.arr);
    }
    else {
      base::run_lerp_mod_multitap_delay<value_type> (
        tap_ptrs, in, spls, mod_spls, out_overwrite, stage.delay, xspan {mod});
    }
    if constexpr (out_overwrite) {
      xspan_memdump (out, xspan {tap0.data(), in.size()});
    }
    else {
      base::span_add (out, xspan {tap0.data(), in.size()});
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Tag, class... Args>
  void run_parallel_delays (xspan<value_type const> in, Tag, Args&&... outs_arg)
  {
    constexpr bool out_overwrite     = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_t const> g = spec::get_gains (Idx);
    constexpr auto                 n_outs       = spec::get_n_outs (Idx);
    constexpr auto                 n_out_reused = out_overwrite ? n_outs : 0;
    static_assert (n_outs <= g.size());
    static_assert (sizeof...(outs_arg) == n_outs);

    // TODO: would an implementation adding in place, aka, not reusing
    // "run_multitap_delay" be better because of the reduced number of
    // buffers/cache?

    // prepare arguments for "run_multitap_delay"'
    std::array<value_type*, n_outs> outs
      = to_ptr_array (std::forward_as_tuple (outs_arg...));
    std::array<block_array, g.size() - n_out_reused> tap_tmp;
    std::array<value_type*, g.size()>                tap_ptrs;
    for (uint i = 0; i < n_out_reused; ++i) {
      tap_ptrs[i] = outs[i];
    }
    for (uint i = 0; i < tap_tmp.size(); ++i) {
      tap_ptrs[n_out_reused + i] = tap_tmp[i].data();
    }
    // call
    std::apply (
      [&] (auto&&... args) {
        run_multitap_delay<Idx> (in, overwrite, args...);
      },
      tap_ptrs);
    // output overwrite aware sum
    for (uint tap = 0; tap < n_out_reused; ++tap) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        // reminder: &out == &tap_ptrs[0]
        outs[tap][i] *= (value_type) g[tap];
      }
    }
    for (uint tap = n_out_reused; tap < tap_ptrs.size(); ++tap) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        outs[tap % n_outs][i] += tap_ptrs[tap][i] * (value_type) g[tap];
      }
    }
  }
  //----------------------------------------------------------------------------
  // args contains n_outs followed by all the modulations
  template <uint Idx, class Tag, class... Ts>
  void run_parallel_mod_delays (xspan<value_type const> in, Tag, Ts&&... args)
  {
    constexpr bool out_overwrite        = !std::is_same_v<Tag, add_to_out_tag>;
    constexpr xspan<fixpt_t const>  g   = spec::get_gains (Idx);
    constexpr xspan<fixpt_spls_mod> mod = spec::get_delays_mod_spls (Idx);
    constexpr auto                  n_outs       = spec::get_n_outs (Idx);
    constexpr auto                  n_out_reused = out_overwrite ? n_outs : 0;
    static_assert (n_outs <= g.size());
    static_assert (sizeof...(args) == (n_outs + g.size()));

    // TODO: would an implementation adding in place, aka, not reusing
    // "run_multitap_mod_delay" be better because of the reduced number of
    // buffers/cache?

    // prepare output pointers to call "run_multitap_mod_delay"
    std::array<value_type*, n_outs> outs = to_ptr_array (
      forward_range_as_tuple<0, n_outs> (std::forward<Ts> (args)...));
    std::array<block_array, g.size() - n_out_reused> tap_tmp;
    std::array<value_type*, g.size()>                tap_ptrs;
    for (uint i = 0; i < n_out_reused; ++i) {
      tap_ptrs[i] = outs[i];
    }
    for (uint i = 0; i < tap_tmp.size(); ++i) {
      tap_ptrs[n_out_reused + i] = tap_tmp[i].data();
    }
    // interleave modulator arguments to match "run_multitap_mod_delay"'s
    // signature
    auto mods_tpl = forward_range_as_tuple<n_outs> (std::forward<Ts> (args)...);
    // intersperse, inserts value_type* at odd positions, add something to erase
    using args_tmp1 = mp11::mp_push_front<decltype (mods_tpl), void>;
    using args_tmp2 = mp11::mp_intersperse<args_tmp1, value_type*>;
    // erase trailing dummy void type
    using args_tpl = mp11::mp_pop_front<args_tmp2>;

    args_tpl mtm_delay_args;
    mp11::mp_for_each<mp11::mp_iota_c<n_outs>> ([&] (auto i) {
      if constexpr ((i % 2) == 0) {
        std::get<i> (mtm_delay_args) = tap_ptrs[i / 2];
      }
      else {
        std::get<i> (mtm_delay_args) = std::get<i / 2> (mods_tpl);
      }
    });
    // do the call itself
    std::apply (
      [&] (auto&&... args) {
        run_multitap_mod_delay<Idx> (in, overwrite, args...);
      },
      mtm_delay_args);

    for (uint tap = 0; tap < n_out_reused; ++tap) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        // reminder: &out == &tap_ptrs[0]
        outs[tap][i] *= (value_type) g[tap];
      }
    }
    for (uint tap = n_out_reused; tap < tap_ptrs.size(); ++tap) {
      ARTV_LOOP_UNROLL_SIZE_HINT (16)
      for (uint i = 0; i < in.size(); ++i) {
        outs[tap % n_outs][i] += tap_ptrs[tap][i] * (value_type) g[tap];
      }
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Range>
  void run_variable_delay (
    xspan<value_type>       out,
    xspan<value_type const> in,
    Range&&                 gen)
  {
    constexpr uint min   = spec::get_min_delay_spls (Idx);
    constexpr uint max   = spec::get_max_delay_spls (Idx);
    constexpr uint range = max - min;
    auto           rgen
      = base::get_sample_generator<value_type> (std::forward<Range> (gen));

    static_assert (min >= max_block_size);

    ARTV_LOOP_UNROLL_SIZE_HINT (16)
    for (uint i = 0; i < in.size(); ++i) {
      value_type v;
      uint       spls;
      auto       gv = rgen (i);
      if constexpr (spec::uses_absolute_delay_modulation (Idx)) {
        spls = (uint) gv;
      }
      else if constexpr (is_fixpt_v<value_type>) {
        spls = (uint) (min + range * gv.to_floatp());
      }
      else {
        spls = (uint) (min + range * gv);
      }
      assert (spls >= min);
      assert (spls <= max);
      v = get<Idx> (spls);
      std::get<Idx> (_stages).delay.push (in[i]);
      out[i] = v; // they could be aliased...
    }
  }
  //----------------------------------------------------------------------------
  template <uint Idx, class Range>
  void run_variable_delay (xspan<value_type> io, Range&& gen)
  {
    run_variable_delay<Idx> (io, io, std::forward<Range> (gen));
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
    uint                    arg = 0;
    std::array<uint, N + 1> ret {};
    for (uint i = 0; i < idxs.size(); ++i) {
      if (spec::is_allpass (idxs[i])) {
        // lfo + gain parameter.
        arg += 2;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::is_1tap_delay (idxs[i])) {
        // lfo for modulation
        arg += 1;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::is_lowpass (idxs[i]) || spec::is_highpass (idxs[i])) {
        // gain parameter.
        arg += 1;
        ret[i + 1] = arg;
        continue;
      }
      else if (spec::get_crossover_n_bands (idxs[i]) != 0) {
        // [freq, gain]... gain
        arg += (spec::get_crossover_n_bands (idxs[i]) * 2) + 1;
        ret[i + 1] = arg;
        continue;
      }
      else {
        assert (false); // Type not supported/added yet
      }
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
  template <uint... Idx, class... Ts>
  void run_nested_ap (
    value_type*             out,
    xspan<value_type const> in,
    Ts&&... argsp)
  {
    constexpr auto idxs        = make_array (Idx...);
    constexpr auto arg_offset  = get_arg_offsets (idxs);
    constexpr uint max_n_args  = arg_offset.back();
    constexpr int  last_ap_pos = get_previous_ap_pos (idxs.size(), idxs);

    static_assert (
      spec::is_allpass (idxs[0]), "1st stage has to be an allpass");
    static_assert (max_n_args >= sizeof...(Ts), "Excess arguments");
    using defaults = mp11::
      mp_repeat_c<std::tuple<defaulted_tag>, max_n_args - sizeof...(Ts)>;

    block_array fwd;
    auto args = std::tuple_cat (std::forward_as_tuple (argsp...), defaults {});
    xspan_memdump (fwd.data(), in);

    using idxlist = mp_list<k_uint<Idx>...>;

    mp_foreach_idx (idxlist {}, [&] (auto order, auto stage_idx) {
      bool constexpr is_serial = false; // not nested element // TBD

      if constexpr (spec::is_allpass (stage_idx) && !is_serial) {
        // having delays bigger than the block size allows to run each stage
        // once for all samples instead of each sample having to loop over
        // each stage.
        static_assert (
          spec::get_min_delay_spls (stage_idx) >= max_block_size,
          "Nested AP delays have to be GE than the block size");

        block_array bwd;
        fetch_mod_en<stage_idx> (
          xspan {bwd.data(), in.size()}, std::get<arg_offset[order]> (args));

        auto&& gain_arg = std::get<arg_offset[order] + 1> (args);
        auto   gain_gen = base::get_gain_generator<value_type> (
          gain_arg, spec::get_gain (stage_idx));

        ARTV_LOOP_UNROLL_SIZE_HINT (16)
        for (uint i = 0; i < in.size(); ++i) {
          // See lattice form:
          // https://www.dsprelated.com/freebooks/pasp/Nested_Allpass_Filters.html
          std::tie (bwd[i], fwd[i]) = base::run_allpass<value_type> (
            fwd[i], bwd[i], gain_gen (i), std::get<stage_idx> (_stages).state);
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

          // process backward signal with nested non-allpass elements
          mp11::mp_for_each<mp11::mp_iota_c<n_steps>> ([&] (auto i) {
            constexpr uint order_rev     = order - i - 1; // reverse iter
            constexpr uint rev_idx       = idxs[order_rev];
            constexpr bool is_serial_rev = false; // TBD
            if constexpr (is_serial_rev) {
              // skip... not belonging to the backwards path
            }
            else if constexpr (
              spec::is_highpass (rev_idx) || spec::is_lowpass (rev_idx)) {
              auto gain = std::get<arg_offset[order_rev]> (args);
              run (stage_list<rev_idx> {}, xspan {bwd.data(), in.size()}, gain);
            }
            else if constexpr (spec::is_1tap_delay (rev_idx)) {
              auto lfo_gen_rev = base::get_generic_generator<value_type> (
                std::get<arg_offset[order_rev]> (args));
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                lfo_gen_rev);
            }
            else if constexpr (spec::get_crossover_n_bands (rev_idx) == 1) {
              static constexpr uint a1 = arg_offset[order_rev];
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                std::get<a1> (args),
                std::get<a1 + 1> (args),
                std::get<a1 + 2> (args));
            }
            else if constexpr (spec::get_crossover_n_bands (rev_idx) == 2) {
              static constexpr uint a1 = arg_offset[order_rev];
              run (
                stage_list<rev_idx> {},
                xspan {bwd.data(), in.size()},
                std::get<a1> (args),
                std::get<a1 + 1> (args),
                std::get<a1 + 2> (args),
                std::get<a1 + 3> (args),
                std::get<a1 + 4> (args));
            }
            else {
              static_assert (sizeof (order_rev), "Not implemented yet");
            }
          });
          // insert the backwards signal on the previous allpass...
          std::get<idxs[ap_pos]> (_stages).delay.push_block (
            xspan {bwd.data(), in.size()});
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
          run (stage_list<stage_idx> {}, xspan {fwd.data(), in.size()}, gain);
        }
        else if constexpr (spec::is_1tap_delay (stage_idx)) {
          auto lfo_gen_rev = base::get_generic_generator<value_type> (
            std::get<arg_offset[order]> (args));
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            lfo_gen_rev);
        }
        else if constexpr (spec::get_crossover_n_bands (stage_idx) == 1) {
          static constexpr uint a1 = arg_offset[order];
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            std::get<a1> (args),
            std::get<a1 + 1> (args),
            std::get<a1 + 2> (args));
        }
        else if constexpr (spec::get_crossover_n_bands (stage_idx) == 2) {
          static constexpr uint a1 = arg_offset[order];
          run (
            stage_list<stage_idx> {},
            xspan {fwd.data(), in.size()},
            std::get<a1> (args),
            std::get<a1 + 1> (args),
            std::get<a1 + 2> (args),
            std::get<a1 + 3> (args),
            std::get<a1 + 4> (args));
        }
        else {
          static_assert (sizeof (stage_idx), "Not implemented yet");
        }
      }
      else {
        // skipping nested non-allpass elements, these are handdebugled inside
        // the. allpass conditional.
      }
    });
    // the last allpass gets the fwd signal enqueued
    static_assert (last_ap_pos > 0 || !spec::is_allpass (idxs[1]));
    std::get<idxs[last_ap_pos]> (_stages).delay.push_block (
      xspan {fwd.data(), in.size()});
  }
  //----------------------------------------------------------------------------
  // cast handling ratios
  template <class T, class U>
  T arith_cast (U v)
  {
    if constexpr (is_ratio_v<U>) {
      return T {} + v;
    }
    else {
      return static_cast<T> (v);
    }
  }
  //----------------------------------------------------------------------------
  using block_array = std::array<value_type, max_block_size>;
  using base        = algo_engine_untemplated;

  using indexes = mp11::mp_iota_c<spec::size()>;

  using states_typelist = std::conditional_t<
    is_fixpt_v<value_type>,
    mp11::mp_transform_q<state::index_to_state_qfn<fixpt_t, spec>, indexes>,
    mp11::mp_transform_q<state::index_to_state_qfn<float, spec>, indexes>>;
  using states_tuple = mp11::mp_rename<states_typelist, std::tuple>;

  using delays_typelist = mp11::
    mp_transform_q<delay::index_to_delay_buffer_qfn<Data_type, spec>, indexes>;

  template <class S, class D>
  struct stage_type {
    S state;
    D delay;
  };

  using stages_tuple = mp_mix<stage_type, states_tuple, delays_typelist>;

  stages_tuple _stages {};
};
}}} // namespace artv::detail::lofiverb
//------------------------------------------------------------------------------
