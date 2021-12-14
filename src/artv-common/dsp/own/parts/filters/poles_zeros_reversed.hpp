#pragma once

#include <complex>

#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/simd_complex.hpp"

namespace artv {

static uint get_reversed_pole_n_stages (
  std::complex<double> pole,
  double               snr_db)
{
  auto   pole_mod = std::abs (pole); // modulus
  double stages   = std::log2 (-snr_db / (20. * std::log10 (pole_mod)));
  return (uint) std::ceil (stages);
}

//------------------------------------------------------------------------------
// time reversed real pole
struct t_rev_rpole {
  // Based on paper from Martin Vicanek "A New Reverse IIR Algorithm".
  // https://www.vicanek.de/articles.htm
  // file:///home/s0001192/Downloads/ReverseIIR.pdf
  //----------------------------------------------------------------------------
  enum coeffs {
    c,
    n_coeffs,
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    // sum of powers = (2^n+1) - 1.
    return (1u << n_stages) - 1u;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, V pole)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));

    vec_store (&co[c * traits.size], pole);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>,
    uint stages)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint n_stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 in,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * get_n_states (n_stages)));
    assert (n_stages >= 1);

    V    c_v    = vec_load<V> (&co[c * traits.size]);
    T*   z_ptr  = &st[0];
    uint z_size = 1;
    V    y      = in;

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = z_size - 1;
      uint pos  = ((z_size + sample_idx) & mask) * traits.size;

      auto z_y = vec_load<V> (&z_ptr[pos]);
      vec_store (&z_ptr[pos], y);
      y = y * c_v + z_y;
      z_ptr += z_size * traits.size;
      z_size *= 2;
      c_v *= c_v;
    }
    return y;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<V>                         io, // ins on call, outs when returning
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    // block version, the intent with it is to run some iterations on the same
    // cache lines before moving to the next.
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * get_n_states (n_stages)));
    assert (n_stages >= 1);

    V    c_v    = vec_load<V> (&co[c * traits.size]);
    T*   z_ptr  = &st[0];
    uint z_size = 1;

    for (uint s = 0; s < n_stages; ++s) {
      uint mask = z_size - 1;
      // process the block stage-wise
      for (uint i = 0; i < io.size(); ++i) {
        V    y   = io[i];
        uint pos = ((z_size + sample_idx + i) & mask) * traits.size;
        auto z_y = vec_load<V> (&z_ptr[pos]);
        vec_store (&z_ptr[pos], y);
        y     = y * c_v + z_y;
        io[i] = y;
      }
      z_ptr += z_size * traits.size;
      z_size *= 2;
      c_v *= c_v;
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// time reversed complex pole
struct t_rev_cpole {
  // Based on paper from Martin Vicanek "A New Reverse IIR Algorithm".
  // https://www.vicanek.de/articles.htm
  // file:///home/s0001192/Downloads/ReverseIIR.pdf
  //----------------------------------------------------------------------------
  enum coeffs {
    a, // re
    b, // im
    n_coeffs,
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return t_rev_rpole::get_n_states (n_stages) * 2;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, vec_complex<V> pole)
  {
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs * traits.size);

    V real = vec_real (pole);
    V imag = vec_imag (pole);

    vec_store (&co[a * traits.size], real);
    vec_store (&co[b * traits.size], imag);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>,
    uint)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint n_stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static vec_complex<V> tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    vec_complex<V>                    in,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * get_n_states (n_stages)));

    V a_v  = vec_load<V> (&co[a * traits.size]);
    V b_v  = vec_load<V> (&co[b * traits.size]);
    V y_re = vec_real (in);
    V y_im = vec_imag (in);

    T*   z_re_ptr = &st[0];
    T*   z_im_ptr;
    uint z_size = 1;

    for (uint i = 0; i < n_stages; ++i) {
      // recomputing imaginary buffer position (done)
      z_im_ptr = z_re_ptr + (z_size * traits.size);
      // computing sample position on each of the same-sized buffers
      uint mask = z_size - 1;
      uint pos  = ((z_size + sample_idx) & mask) * traits.size;
      // calculation
      auto z_re = vec_load<V> (&z_re_ptr[pos]);
      auto z_im = vec_load<V> (&z_im_ptr[pos]);
      vec_store (&z_re_ptr[pos], y_re);
      vec_store (&z_im_ptr[pos], y_im);
      auto y_re_cp = y_re;
      y_re         = (a_v * y_re_cp) - (b_v * y_im) + z_re;
      y_im         = (b_v * y_re_cp) + (a_v * y_im) + z_im;
      // adjusting real buffer position and buffers size for the next stage
      z_re_ptr = z_im_ptr + (z_size * traits.size);
      z_size *= 2;
      // as this is a serial algorithm, the initial assumption is that
      // precomputing the powers of 2 for each stage wouldn't make a lot of
      // difference because the pipelines are probably at low capacity. Favoring
      // less and fixed-size coefficient memory usage, as this is 4mul + 1add
      // that can probably be parallelized by the CPU, it has no dependencies on
      // the code above (TODO: measure).
      auto a_cp = a_v;
      a_v       = a_cp * a_cp - b_v * b_v;
      b_v       = (T) 2. * a_cp * b_v;
    }
    return vec_complex<V> (y_re, y_im);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<V>                         io_re, // ins on call, outs when returning
    crange<V>                         io_im, // ins on call, outs when returning
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    uint                              n_stages,
    uint                              sample_idx) // sample counter (external)
  {
    // block version, the intent with it is to run some iterations on the same
    // cache lines before moving to the next.

    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));
    assert (st.size() >= (traits.size * get_n_states (n_stages)));
    assert (io_re.size() <= io_im.size());

    V a_v = vec_load<V> (&co[a * traits.size]);
    V b_v = vec_load<V> (&co[b * traits.size]);

    T*   z_re_ptr = &st[0];
    T*   z_im_ptr;
    uint z_size = 1;

    for (uint i = 0; i < n_stages; ++i) {
      // recomputing imaginary buffer position (done)
      z_im_ptr = z_re_ptr + (z_size * traits.size);
      // computing sample position on each of the same-sized buffers
      uint mask = z_size - 1;

      for (uint i = 0; i < io_re.size(); ++i) {
        uint pos  = ((z_size + sample_idx + i) & mask) * traits.size;
        V    y_re = io_re[i];
        V    y_im = io_im[i];
        // calculation
        auto z_re = vec_load<V> (&z_re_ptr[pos]);
        auto z_im = vec_load<V> (&z_im_ptr[pos]);
        vec_store (&z_re_ptr[pos], y_re);
        vec_store (&z_im_ptr[pos], y_im);
        auto y_re_cp = y_re;
        io_re[i]     = (a_v * y_re_cp) - (b_v * y_im) + z_re;
        io_im[i]     = (b_v * y_re_cp) + (a_v * y_im) + z_im;
      }
      // adjusting real buffer position and buffers size for the next stage
      z_re_ptr = z_im_ptr + (z_size * traits.size);
      z_size *= 2;
      // as this is a serial algorithm, the initial assumption is that
      // precomputing the powers of 2 for each stage wouldn't make a lot of
      // difference because the pipelines are probably at low capacity. Favoring
      // less and fixed-size coefficient memory usage, as this is 4mul + 1add
      // that can probably be parallelized by the CPU, it has no dependencies on
      // the code above (TODO: measure).
      auto a_cp = a_v;
      a_v       = a_cp * a_cp - b_v * b_v;
      b_v       = (T) 2. * a_cp * b_v;
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
// time reversed complex conjugate pole pair
struct t_rev_ccpole_pair {
  // Based on paper from Martin Vicanek "A New Reverse IIR Algorithm".
  // https://www.vicanek.de/articles.htm
  // file:///home/s0001192/Downloads/ReverseIIR.pdf
  //
  // complex conjugate pair optimization.
  //----------------------------------------------------------------------------
  enum coeffs {
    ratio = t_rev_cpole::n_coeffs,
    n_coeffs,
  };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return t_rev_cpole::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (crange<vec_value_type_t<V>> co, vec_complex<V> pole)
  {
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= n_coeffs * traits.size);

    t_rev_cpole::reset_coeffs (co, pole);
    vec_store (&co[ratio * traits.size], vec_real (pole) / vec_imag (pole));
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>,
    uint)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint n_stages)
  {
    t_rev_cpole::reset_states (st, n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 in_re,
    uint                              n_stages, // sample counter (external)
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));

    vec_complex<V> y = t_rev_cpole::tick (
      co, st, vec_complex<V> (in_re), n_stages, sample_idx);

    return vec_real (y) + vec_load<V> (&co[ratio * traits.size]) * vec_imag (y);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<V>                         io, // ins on call, outs when returning
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    uint                              n_stages, // sample counter (external)
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    assert (co.size() >= (traits.size * n_coeffs));

    std::array<V, 32> imag;
    auto              ratio_v = vec_load<V> (&co[ratio * traits.size]);

    while (io.size() > 0) {
      uint blocksize = std::min (io.size(), imag.size());
      memset (&imag[0], 0, blocksize * sizeof imag[0]);
      t_rev_cpole::tick (
        make_crange (io.data(), blocksize),
        make_crange (imag),
        co,
        st,
        n_stages,
        sample_idx);

      for (uint i = 0; i < blocksize; ++i) {
        io[i] += ratio_v * imag[i];
      }
      sample_idx += blocksize;
      io = io.shrink_head (blocksize);
    }
  }
  //----------------------------------------------------------------------------
};
//----------------------------------------------------------------------------
// one time reversed real pole one real zero filter
struct t_rev_rpole_rzero {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + t_rev_rpole::n_coeffs };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return rzero::n_states + t_rev_rpole::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    V                           re_pole,
    V                           re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    t_rev_rpole::reset_coeffs (co, re_pole);
    co = co.shrink_head (t_rev_rpole::n_coeffs * traits.size);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint n_stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 x,
    uint                              n_stages,
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    V out = t_rev_rpole::tick (co, st, x, n_stages, sample_idx);
    co    = co.shrink_head (t_rev_rpole::n_coeffs * traits.size);
    st    = st.shrink_head (t_rev_rpole::get_n_states (n_stages) * traits.size);
    return rzero::tick (co, st, out);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<V>                         io, // ins on call, outs when returning
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    uint                              n_stages,
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    t_rev_rpole::tick (io, co, st, n_stages, sample_idx);
    co = co.shrink_head (t_rev_rpole::n_coeffs * traits.size);
    st = st.shrink_head (t_rev_rpole::get_n_states (n_stages) * traits.size);

    for (uint i = 0; i < io.size(); ++i) {
      io[i] = rzero::tick (co, st, io[i]);
    }
  }
  //----------------------------------------------------------------------------
};
// time reversed complex conjugate poles pair + two real zeros filter ----------
struct t_rev_ccpole_pair_rzero_pair {
  //----------------------------------------------------------------------------
  enum coeffs { n_coeffs = rzero::n_coeffs + t_rev_ccpole_pair::n_coeffs };
  //----------------------------------------------------------------------------
  static constexpr uint get_n_states (uint n_stages)
  {
    return 2 * rzero::n_states + t_rev_ccpole_pair::get_n_states (n_stages);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_coeffs (
    crange<vec_value_type_t<V>> co,
    vec_complex<V>              pole,
    V                           re_zero)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    t_rev_ccpole_pair::reset_coeffs (co, pole);
    co = co.shrink_head (t_rev_ccpole_pair::n_coeffs * traits.size);
    rzero::reset_coeffs (co, re_zero);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void fix_unsmoothable_coeffs (
    crange<vec_value_type_t<V>>,
    crange<vec_value_type_t<const V>>)
  {}
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<vec_value_type_t<V>> st, uint n_stages)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    uint numstates = traits.size * get_n_states (n_stages);
    assert (st.size() >= numstates);
    memset (st.data(), 0, sizeof (T) * numstates);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static V tick (
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    V                                 x,
    uint                              n_stages,
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    auto out = x;

    out = t_rev_ccpole_pair::tick (co, st, out, n_stages, sample_idx);
    co  = co.shrink_head (t_rev_ccpole_pair::n_coeffs * traits.size);
    st  = st.shrink_head (
      t_rev_ccpole_pair::get_n_states (n_stages) * traits.size);

    for (uint i = 0; i < 2; ++i) {
      out = rzero::tick (co, st, out);
      // same zero location
      st = st.shrink_head (rzero::n_states * traits.size);
    }
    return out;
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void tick (
    crange<V>                         io, // ins on call, outs when returning
    crange<const vec_value_type_t<V>> co,
    crange<vec_value_type_t<V>>       st,
    uint                              n_stages,
    uint                              sample_idx)
  {
    using T               = vec_value_type_t<V>;
    constexpr auto traits = vec_traits<V>();

    t_rev_ccpole_pair::tick (io, co, st, n_stages, sample_idx);
    co = co.shrink_head (t_rev_ccpole_pair::n_coeffs * traits.size);
    st = st.shrink_head (
      t_rev_ccpole_pair::get_n_states (n_stages) * traits.size);

    for (uint j = 0; j < 2; ++j) {
      for (uint i = 0; i < io.size(); ++i) {
        io[i] = rzero::tick (co, st, io[i]);
      }
      st = st.shrink_head (rzero::n_states * traits.size);
    }
  }
  //----------------------------------------------------------------------------
};
} // namespace artv
