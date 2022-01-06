#pragma once

// Time reversed filters. All on this header is based on the paper "A New
// Reverse IIR Algorithm" from Martin Vicanek.
// https://www.vicanek.de/articles/ReverseIIR.pdf

#include <complex>
#include <type_traits>

#include "artv-common/dsp/own/parts/filters/poles_zeros.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

namespace artv {
//------------------------------------------------------------------------------
static uint get_reversed_pole_n_stages (
  std::complex<double> pole,
  double               snr_db)
{
  auto   pole_mod = std::abs (pole); // modulus
  double stages   = std::log2 (-snr_db / (20. * std::log10 (pole_mod)));
  return (uint) std::ceil (stages);
}
//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------
struct t_rev_pole_stage_op {
  template <class T>
  static T run (T c, T in, T& state)
  {
    T& z1 = state; // just documentation...
    T  y  = c * in + z1;
    z1    = in;
    return y;
  }
};
//------------------------------------------------------------------------------
struct t_rev_zero_stage_op {
  template <class T>
  static T run (T c, T in, T& state)
  {
    // This requires no delay line to perform time reversal.
    T& z1 = state; // just documentation...
    T  y  = in - c * z1;
    z1    = in;
    return y;
  }
};
//------------------------------------------------------------------------------
template <class Stage_op, bool is_complex>
struct t_rev_single {
  //----------------------------------------------------------------------------
  template <class V>
  using value_type = std::conditional_t<is_complex, vec_complex<V>, V>;
  //----------------------------------------------------------------------------
  enum coeffs {
    n_coeffs = is_complex ? 2 : 1,
  };
  //----------------------------------------------------------------------------
  enum coeffs_int { n_coeffs_int };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    // sum of powers = (2^n+1) - 1.
    constexpr uint factor = is_complex ? 2 : 1;
    return ((1u << n_stages) - 1u) * factor;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, value_type<V> pole_or_zero)
  {
    assert (co.size() >= n_coeffs);
    auto c = (value_type<V>*) co.data();
    *c     = pole_or_zero;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static value_type<V> tick (
    crange<const V> co,
    crange<V>       st,
    value_type<V>   in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    auto c       = *((value_type<V>*) &co[0]);
    auto st_ptr  = (value_type<V>*) &st[0];
    uint st_size = 1;
    auto y       = in;

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = st_size - 1;
      uint pos  = (st_size + sample_idx) & mask;

      y = Stage_op::run (c, y, st_ptr[pos]);

      st_ptr += st_size;
      st_size *= 2;
      c *= c;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<value_type<V>> tick (
    crange<const V>       co,
    crange<V>             st,
    crange<value_type<V>> io, // ins on call, outs when returning
    uint                  n_stages,
    uint                  sample_idx) // sample counter (external)
  {
    // block version, the intent with it is to run some iterations on the same
    // cache/delay line before moving to the next. This is a theoretically
    // better memory access pattern once the delay lines become separated enough
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    auto c       = *((value_type<V>*) co.data());
    auto st_ptr  = (value_type<V>*) &st[0];
    uint st_size = 1;

    for (uint s = 0; s < n_stages; ++s) {
      uint mask = st_size - 1;
      // process the block stage-wise
      for (uint i = 0; i < io.size(); ++i) {
        uint pos = (st_size + sample_idx + i) & mask;
        io[i]    = Stage_op::run (c, io[i], st_ptr[pos]);
      }
      st_ptr += st_size;
      st_size *= 2;
      c *= c;
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class Stage_op>
struct t_rev_conjugate {
  using base = t_rev_single<Stage_op, true>;

  enum coeffs {
    ratio = base::n_coeffs,
    n_coeffs,
  };
  enum coeffs_int { n_coeffs_int };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return base::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, vec_complex<V> pole_or_zero)
  {
    assert (co.size() >= n_coeffs);
    base::reset_coeffs (co, pole_or_zero);
    co[ratio] = pole_or_zero.re / pole_or_zero.im;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               in_re,
    uint            n_stages, // sample counter (external)
    uint            sample_idx)
  {
    assert (co.size() >= n_coeffs);
    vec_complex<V> y
      = base::tick (co, st, vec_complex<V> (in_re), n_stages, sample_idx);
    return y.re + co[ratio] * y.im;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages, // sample counter (external)
    uint            sample_idx)
  {
    assert (co.size() >= n_coeffs);

    std::array<vec_complex<V>, 64> io_c;
    auto                           ratio_v = co[ratio];

    for (uint offset = 0; offset < io.size(); offset += io_c.size()) {
      uint blocksize = std::min<uint> (io_c.size(), io.size() - offset);
      // interleave
      for (uint i = 0; i < blocksize; ++i) {
        io_c[i] = vec_complex<V> {io[offset + i]};
      }
      // process
      base::tick (
        co,
        st,
        make_crange (io_c.data(), blocksize),
        n_stages,
        sample_idx + offset);
      // deinterleave
      for (uint i = 0; i < blocksize; ++i) {
        io[offset + i] = io_c[i].re + ratio_v * io_c[i].im;
      }
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// Time reversed partial fraction expansion pole pair.
// A class to be used when the poles can be mutate between:
// -equal poles
// -real poles
// -complex conjugate poles
template <bool is_complex>
struct t_rev_pfe_pole_pair {
  //----------------------------------------------------------------------------
  using base = t_rev_single<t_rev_pole_stage_op, is_complex>;
  //----------------------------------------------------------------------------
  template <class T>
  using value_type = typename base::template value_type<T>;
  //----------------------------------------------------------------------------
  enum coeffs {
    k_coeff_first = base::n_coeffs * 2,
    n_coeffs      = k_coeff_first * 2,
  };
  //----------------------------------------------------------------------------
  enum coeffs_int { n_coeffs_int };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return base::get_n_states (n_stages) * 2;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>     co,
    value_type<V> pole1,
    value_type<V> pole2)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    auto co_ptr      = (value_type<V>*) co.data();
    co_ptr[co_pole1] = pole1;
    co_ptr[co_pole2] = pole2;
    co_ptr[co_k1]    = (T) 1 / ((T) 1 - pole2 / pole1);
    co_ptr[co_k2]    = (T) 1 / ((T) 1 - pole1 / pole2);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    static constexpr uint n_interleavings = 2;

    auto co_ptr = (value_type<V>*) co.data();
    auto c1     = co_ptr[co_pole1];
    auto c2     = co_ptr[co_pole2];
    auto k1     = co_ptr[co_k1];
    auto k2     = co_ptr[co_k2];
    auto z_ptr  = (value_type<V>*) &st[0];
    uint z_size = 1;
    auto y1     = value_type<V> {in};
    auto y2     = value_type<V> {in};

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = z_size - 1;
      uint pos  = ((z_size + sample_idx) & mask) * n_interleavings;

      y1 = t_rev_pole_stage_op::run (c1, y1, z_ptr[pos]);
      y2 = t_rev_pole_stage_op::run (c2, y2, z_ptr[pos + 1]);

      c1 *= c1;
      c2 *= c2;

      z_ptr += z_size * n_interleavings;
      z_size *= 2;
    }

    auto pfe = y1 * k1 + y2 * k2;
    if constexpr (is_complex) {
      return pfe.re;
    }
    else {
      return pfe;
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {

    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    static constexpr uint n_interleavings = 2;
    constexpr uint        pole1           = 0;
    constexpr uint        pole2           = 1;

    std::array<std::array<value_type<V>, 2>, 64> buff;

    auto co_ptr = (value_type<V>*) co.data();
    auto k1     = co_ptr[co_k1];
    auto k2     = co_ptr[co_k2];

    for (uint offset = 0; offset < io.size(); offset += buff.size()) {
      uint blocksize = std::min<uint> (buff.size(), io.size() - offset);

      // prepare new inputs
      for (uint i = 0; i < blocksize; ++i) {
        buff[i][pole1] = buff[i][pole2] = value_type<V> {io[offset + i]};
      }

      auto c1     = co_ptr[co_pole1];
      auto c2     = co_ptr[co_pole2];
      auto z_ptr  = (value_type<V>*) &st[0];
      uint z_size = 1;

      // process a block stage-wise
      for (uint s = 0; s < n_stages; ++s) {
        uint mask     = z_size - 1;
        uint idx_base = z_size + sample_idx + offset;

        for (uint i = 0; i < blocksize; ++i) {
          uint pos = ((idx_base + i) & mask) * n_interleavings;

          auto& y1 = buff[i][pole1];
          auto& y2 = buff[i][pole2];
          y1       = t_rev_pole_stage_op::run (c1, y1, z_ptr[pos]);
          y2       = t_rev_pole_stage_op::run (c2, y2, z_ptr[pos + 1]);
        }

        z_ptr += z_size * n_interleavings;
        z_size *= 2;
        c1 *= c1;
        c2 *= c2;
      }
      // store outs
      for (uint i = 0; i < blocksize; ++i) {
        auto pfe = buff[i][pole1] * k1 + buff[i][pole2] * k2;
        if constexpr (is_complex) {
          io[offset + i] = pfe.re;
        }
        else {
          io[offset + i] = pfe;
        }
      }
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
private:
  enum coeffs_in_ptr {
    co_pole1,
    co_pole2,
    co_k1,
    co_k2,
  };
};
#if 0
//------------------------------------------------------------------------------
// Time reversed partial fraction expansion pole pair.
// A class to be used when the poles can be mutate between:
// -equal poles
// -real poles
// -complex conjugate poles
template <class Stage_op, bool is_complex>
struct t_rev_naive_cascade_pair {
  //----------------------------------------------------------------------------
  using base = t_rev_single<Stage_op, is_complex>;

  template <class T>
  using value_type = typename base::template value_type<T>;
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = base::n_coeffs * 2 };
  //----------------------------------------------------------------------------
  enum coeffs_int { n_coeffs_int };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return base::get_n_states (n_stages) * 2;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>     co,
    value_type<V> zero_or_pole1,
    value_type<V> zero_or_pole2)
  {
    base::reset_coeffs (co, zero_or_pole1);
    co = co.shrink_head (base::n_coeffs);
    base::reset_coeffs (co, zero_or_pole2);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    auto out = base::tick (co, st, value_type<V> {in}, n_stages, sample_idx);
    co       = co.shrink_head (base::n_coeffs);
    st       = st.shrink_head (base::get_n_states (n_stages));
    out      = base::tick (co, st, out, n_stages, sample_idx);
    if constexpr (is_complex) {
      return out.re;
    }
    else {
      return out;
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    [[maybe_unused]] std::array<value_type<V>, 64> buff_mem;
    crange<value_type<V>>                          buff;
    if constexpr (is_complex) {
      buff = make_crange (buff_mem);
    }
    else {
      buff = io;
    }

    auto co2 = co.shrink_head (base::n_coeffs);
    auto st2 = st.shrink_head (base::get_n_states (n_stages));

    for (uint offset = 0; offset < io.size(); offset += buff.size()) {
      auto block = make_crange (
        &buff[0], std::min<uint> (buff.size(), io.size() - offset));

      if constexpr (is_complex) {
        for (uint i = 0; i < block.size(); ++i) {
          block[i] = vec_complex<V> {io[offset + i]};
        }
      }
      base::template tick<V> (co, st, block, n_stages, sample_idx + offset);
      base::template tick<V> (co2, st2, block, n_stages, sample_idx + offset);

      if constexpr (is_complex) {
        for (uint i = 0; i < block.size(); ++i) {
          io[offset + i] = block[i].re;
        }
      }
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
};
#endif
//------------------------------------------------------------------------------
} // namespace detail
//------------------------------------------------------------------------------
// t_rev = time reversed
using t_rev_rpole = detail::t_rev_single<detail::t_rev_pole_stage_op, false>;
using t_rev_cpole = detail::t_rev_single<detail::t_rev_pole_stage_op, true>;
// time reversed complex conjugate pole pair
using t_rev_ccpole_pair = detail::t_rev_conjugate<detail::t_rev_pole_stage_op>;
// The poles pairs run in parallel and use partial fraction expansion
using t_rev_rpole_pair = detail::t_rev_pfe_pole_pair<false>;
using t_rev_cpole_pair = detail::t_rev_pfe_pole_pair<true>;

//------------------------------------------------------------------------------
// realtime switches between two real poles or two complex conjugates from the
// same internal state.
class t_rev_pole_pair {
public:
  //----------------------------------------------------------------------------
  enum coeffs {
    co_c1,
    co_c2,
    co_k1,
    co_k2,
    co_ratio,
    n_coeffs,
  };
  enum coeffs_int { n_coeffs_int };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return t_rev_rpole::get_n_states (n_stages) * 3;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<V>      co,
    vec_complex<V> pole1,
    vec_complex<V> pole2)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);

    if (pole1.im[0] == 0) {
#ifndef NDEBUG
      for (uint i = 0; i < vec_traits<V>().size; ++i) {
        assert (pole1.im[i] == (T) 0. || pole1.im[i] == (T) -0.);
        assert (pole2.im[i] == (T) 0. || pole2.im[i] == (T) -0.);
        assert (
          (pole1.re[i] != pole2.re[i] || pole1.re[i] == (T) 0.
           || pole2.re[i] == (T) 0.)
          && "Handle the equal pole case (if possible)");
      }
#endif
      // two real poles, assuming no different Q.
      co[co_c1] = pole1.re;
      co[co_c2] = pole2.re;
    }
    else {
#ifndef NDEBUG
      for (uint i = 0; i < vec_traits<V>().size; ++i) {
        assert (pole1.im[i] != (T) 0. && pole1.im[i] != (T) -0.);
        assert (pole2.im[i] != (T) 0. && pole2.im[i] != (T) -0.);
      }
#endif
      // complex conjugate
      co[co_c1] = pole1.re;
      co[co_c2] = pole1.im;
    }
    // these are always calculated for easier smoothing
    co[co_ratio] = pole1.re / pole1.im;
    co[co_k1]    = (T) 1 / ((T) 1 - pole2.re / pole1.re);
    co[co_k2]    = (T) 1 / ((T) 1 - pole1.re / pole2.re);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               in,
    uint            n_stages,
    uint            sample_idx, // sample counter (external)
    bool            is_real) // external, as for now parameters can be smoothed
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    if (!is_real) {
      return tick_cc (co, st, in, n_stages, sample_idx);
    }
    else {
      return tick_pfe (co, st, in, n_stages, sample_idx);
    }
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx, // sample counter (external)
    bool            is_real) // external, as for now parameters can be smoothed
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    if (!is_real) {
      return tick_cc (co, st, io, n_stages, sample_idx);
    }
    else {
      return tick_pfe (co, st, io, n_stages, sample_idx);
    }
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick_pfe (
    crange<const V> co,
    crange<V>       st,
    V               in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    using T = vec_value_type_t<V>;

    auto c1     = co[co_c1];
    auto c2     = co[co_c1];
    auto k1     = co[co_k1];
    auto k2     = co[co_k2];
    auto z_ptr  = &st[0];
    uint z_size = 1;
    auto y1     = in;
    auto y2     = in;

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = z_size - 1;
      uint pos  = ((z_size + sample_idx) & mask) * n_lanes;

      // Interleaving/lane 0: complex pole real or pole1 real part
      // Interleaving/lane 1: complex pole imaginary part
      // Interleaving/lane 2: pole2 real part
      y1 = detail::t_rev_pole_stage_op::run (c1, y1, z_ptr[pos + 0]);
      z_ptr[pos + 1] = vec_set<V> ((T) 0);
      y2 = detail::t_rev_pole_stage_op::run (c2, y2, z_ptr[pos + 2]);

      c1 *= c1;
      c2 *= c2;

      z_ptr += z_size * n_lanes;
      z_size *= 2;
    }

    auto pfe = y1 * k1 + y2 * k2;
    return pfe;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick_pfe (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    using T = vec_value_type_t<V>;

    constexpr uint pole1 = 0;
    constexpr uint pole2 = 1;

    std::array<std::array<V, 2>, 64> buff;

    auto k1 = co[co_k1];
    auto k2 = co[co_k2];

    for (uint offset = 0; offset < io.size(); offset += buff.size()) {
      uint blocksize = std::min<uint> (buff.size(), io.size() - offset);

      // prepare new inputs
      for (uint i = 0; i < blocksize; ++i) {
        buff[i][pole1] = buff[i][pole2] = io[offset + i];
      }

      auto c1     = co[co_c1];
      auto c2     = co[co_c2];
      auto z_ptr  = &st[0];
      uint z_size = 1;

      // process a block stage-wise
      for (uint s = 0; s < n_stages; ++s) {
        uint mask     = z_size - 1;
        uint idx_base = z_size + sample_idx + offset;

        for (uint i = 0; i < blocksize; ++i) {
          uint pos = ((idx_base + i) & mask) * n_lanes;

          auto& y1 = buff[i][pole1];
          auto& y2 = buff[i][pole2];
          // Interleaving/lane 0: complex pole real or pole1 real part
          // Interleaving/lane 1: complex pole imaginary part
          // Interleaving/lane 2: pole2 real part
          y1 = detail::t_rev_pole_stage_op::run (c1, y1, z_ptr[pos + 0]);
          z_ptr[pos + 1] = vec_set<V> ((T) 0);
          y2 = detail::t_rev_pole_stage_op::run (c2, y2, z_ptr[pos + 2]);
        }

        z_ptr += z_size * n_lanes;
        z_size *= 2;
        c1 *= c1;
        c2 *= c2;
      }
      // store outs
      for (uint i = 0; i < blocksize; ++i) {
        auto pfe       = buff[i][pole1] * k1 + buff[i][pole2] * k2;
        io[offset + i] = pfe;
      }
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick_cc (
    crange<const V> co,
    crange<V>       st,
    V               in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    using T = vec_value_type_t<V>;

    auto c = vec_load<vec_complex<V>> (&co[co_c1]);

    auto z_ptr  = &st[0];
    uint z_size = 1;
    auto y      = in;

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = z_size - 1;
      uint pos  = ((z_size + sample_idx) & mask) * n_lanes;

      // Interleaving/lane 0: complex pole real or pole1 real part
      // Interleaving/lane 1: complex pole imaginary part
      // Interleaving/lane 2: pole2 real part

      auto z = vec_load<vec_complex<V>> (&z_ptr[pos + 0]);
      y      = detail::t_rev_pole_stage_op::run (c, y, z);
      vec_store<vec_complex<V>> (&z_ptr[pos + 0], z);
      z_ptr[pos + 2] = z_ptr[pos];

      c *= c;
      z_ptr += z_size * n_lanes;
      z_size *= 2;
    }

    auto ratio = co[co_ratio];
    return y.re + ratio * y.im;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick_cc (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {

    std::array<vec_complex<V>, 64> io_c;
    auto                           ratio = co[co_ratio];

    for (uint offset = 0; offset < io.size(); offset += io_c.size()) {
      uint blocksize = std::min<uint> (io_c.size(), io.size() - offset);
      // interleave
      for (uint i = 0; i < blocksize; ++i) {
        io_c[i] = vec_complex<V> {io[offset + i]};
      }
      // process a block stage-wise
      auto c      = vec_load<vec_complex<V>> (&co[co_c1]);
      auto z_ptr  = &st[0];
      uint z_size = 1;

      for (uint s = 0; s < n_stages; ++s) {
        uint mask     = z_size - 1;
        uint idx_base = z_size + sample_idx + offset;

        for (uint i = 0; i < blocksize; ++i) {
          uint pos = ((idx_base + i) & mask) * n_lanes;

          // Interleaving/lane 0: complex pole real or pole1 real part
          // Interleaving/lane 1: complex pole imaginary part
          // Interleaving/lane 2: pole2 real part
          auto z  = vec_load<vec_complex<V>> (&z_ptr[pos + 0]);
          io_c[i] = detail::t_rev_pole_stage_op::run (c, io_c[i], z);
          vec_store<vec_complex<V>> (&z_ptr[pos + 0], z);
          z_ptr[pos + 2] = z_ptr[pos];
        }

        z_ptr += z_size * n_lanes;
        z_size *= 2;
        c *= c;
      }
      // deinterleave
      for (uint i = 0; i < blocksize; ++i) {
        io[offset + i] = io_c[i].re + ratio * io_c[i].im;
      }
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
  enum delay_lanes {
    lane_rpole1_or_cpole_re,
    lane_cpole_im,
    lane_rpole2,
    n_lanes
  };
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// one time reversed real pole one real zero filter
struct t_rev_rpole_rzero {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + t_rev_rpole::n_coeffs };
  //----------------------------------------------------------------------------
  enum coeffs_int {
    n_coeffs_int = rzero::n_coeffs_int + t_rev_rpole::n_coeffs_int
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return rzero::n_states + t_rev_rpole::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, V re_pole, V re_zero)
  {
    t_rev_rpole::reset_coeffs (co, re_pole);
    co = co.shrink_head (t_rev_rpole::n_coeffs);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               x,
    uint            n_stages,
    uint            sample_idx)
  {
    V out = t_rev_rpole::tick (co, st, x, n_stages, sample_idx);
    co    = co.shrink_head (t_rev_rpole::n_coeffs);
    st    = st.shrink_head (t_rev_rpole::get_n_states (n_stages));
    return rzero::tick (co, st, out);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx)
  {
    t_rev_rpole::tick (co, st, io, n_stages, sample_idx);
    co = co.shrink_head (t_rev_rpole::n_coeffs);
    st = st.shrink_head (t_rev_rpole::get_n_states (n_stages));

    for (uint i = 0; i < io.size(); ++i) {
      io[i] = rzero::tick (co, st, io[i]);
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
};
// time reversed complex conjugate poles pair + two equal real zeros filter ----
struct t_rev_ccpole_pair_rzero_eq_pair {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + t_rev_ccpole_pair::n_coeffs };
  //----------------------------------------------------------------------------
  enum coeffs_int {
    n_coeffs_int = rzero::n_coeffs_int + t_rev_ccpole_pair::n_coeffs_int
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return 2 * rzero::n_states + t_rev_ccpole_pair::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<V> co, vec_complex<V> pole, V re_zero)
  {
    t_rev_ccpole_pair::reset_coeffs (co, pole);
    co = co.shrink_head (t_rev_ccpole_pair::n_coeffs);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    uint numstates = get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (V) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const V> co,
    crange<V>       st,
    V               x,
    uint            n_stages,
    uint            sample_idx)
  {
    auto out = x;

    out = t_rev_ccpole_pair::tick (co, st, out, n_stages, sample_idx);
    co  = co.shrink_head (t_rev_ccpole_pair::n_coeffs);
    st  = st.shrink_head (t_rev_ccpole_pair::get_n_states (n_stages));

    for (uint i = 0; i < 2; ++i) {
      out = rzero::tick (co, st, out);
      // same zero location
      st = st.shrink_head (rzero::n_states);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static crange<V> tick (
    crange<const V> co,
    crange<V>       st,
    crange<V>       io, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx)
  {
    t_rev_ccpole_pair::tick (co, st, io, n_stages, sample_idx);
    co = co.shrink_head (t_rev_ccpole_pair::n_coeffs);
    st = st.shrink_head (t_rev_ccpole_pair::get_n_states (n_stages));

    for (uint j = 0; j < 2; ++j) {
      for (uint i = 0; i < io.size(); ++i) {
        io[i] = rzero::tick (co, st, io[i]);
      }
      st = st.shrink_head (rzero::n_states);
    }
    return io; // forwarding
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class Base, uint Max_stages>
struct make_max_stages_t_rev : private Base {
  using base = Base;

  using base::get_n_states;
  using base::n_coeffs;
  using base::n_coeffs_int;
  using base::reset_coeffs;
  using base::reset_states;
  using base::tick;

  static constexpr uint max_stages = Max_stages;
  static constexpr uint n_states   = base::get_n_states (max_stages);
};
//------------------------------------------------------------------------------
} // namespace artv
