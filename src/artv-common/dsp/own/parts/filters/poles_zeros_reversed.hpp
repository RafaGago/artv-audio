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
  static void reset_coeffs (crange<V> co, V pole)
  {
    assert (co.size() >= n_coeffs);
    co[c] = pole;
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

    V    c_v    = co[c];
    V*   z_ptr  = &st[0];
    uint z_size = 1;
    V    y      = in;

    for (uint i = 0; i < n_stages; ++i) {
      uint mask = z_size - 1;
      uint pos  = (z_size + sample_idx) & mask;

      auto z_y   = z_ptr[pos];
      z_ptr[pos] = y;
      y          = y * c_v + z_y;
      z_ptr += z_size;
      z_size *= 2;
      c_v *= c_v;
    }
    return y;
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
    // block version, the intent with it is to run some iterations on the same
    // cache/delay line before moving to the next. This is a theoretically
    // better memory access pattern once the delay lines become separated enough
    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (n_stages >= 1);

    V    c_v    = co[c];
    V*   z_ptr  = &st[0];
    uint z_size = 1;

    for (uint s = 0; s < n_stages; ++s) {
      uint mask = z_size - 1;
      // process the block stage-wise
      for (uint i = 0; i < io.size(); ++i) {
        V    y     = io[i];
        uint pos   = (z_size + sample_idx + i) & mask;
        auto z_y   = z_ptr[pos];
        z_ptr[pos] = y;
        y          = y * c_v + z_y;
        io[i]      = y;
      }
      z_ptr += z_size;
      z_size *= 2;
      c_v *= c_v;
    }
    return io; // forwarding
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
  static void reset_coeffs (crange<V> co, vec_complex<V> pole)
  {
    assert (co.size() >= n_coeffs);

    co[a] = vec_real (pole);
    co[b] = vec_imag (pole);
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
  static vec_complex<V> tick (
    crange<const V> co,
    crange<V>       st,
    vec_complex<V>  in,
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));

    V a_v  = co[a];
    V b_v  = co[b];
    V y_re = vec_real (in);
    V y_im = vec_imag (in);

    V*   z_re_ptr = &st[0];
    V*   z_im_ptr;
    uint z_size = 1;

    for (uint i = 0; i < n_stages; ++i) {
      // recomputing imaginary buffer position (done)
      z_im_ptr = z_re_ptr + z_size;
      // computing sample position on each of the same-sized buffers
      uint mask = z_size - 1;
      uint pos  = (z_size + sample_idx) & mask;
      // calculation
      auto z_re     = z_re_ptr[pos];
      auto z_im     = z_im_ptr[pos];
      z_re_ptr[pos] = y_re;
      z_im_ptr[pos] = y_im;
      auto y_re_cp  = y_re;
      y_re          = (a_v * y_re_cp) - (b_v * y_im) + z_re;
      y_im          = (b_v * y_re_cp) + (a_v * y_im) + z_im;
      // adjusting real buffer position and buffers size for the next stage
      z_re_ptr = z_im_ptr + z_size;
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
    crange<const V> co,
    crange<V>       st,
    crange<V>       io_re, // ins on call, outs when returning
    crange<V>       io_im, // ins on call, outs when returning
    uint            n_stages,
    uint            sample_idx) // sample counter (external)
  {
    // block version, the intent with it is to run some iterations on the same
    // cache/delay line before moving to the next. This is a theoretically
    // better memory access pattern once the delay lines become separated enough

    using T = vec_value_type_t<V>;

    assert (co.size() >= n_coeffs);
    assert (st.size() >= get_n_states (n_stages));
    assert (io_re.size() <= io_im.size());

    V a_v = co[a];
    V b_v = co[b];

    V*   z_re_ptr = &st[0];
    V*   z_im_ptr;
    uint z_size = 1;

    for (uint i = 0; i < n_stages; ++i) {
      // recomputing imaginary buffer position (done)
      z_im_ptr = z_re_ptr + z_size;
      // computing sample position on each of the same-sized buffers
      uint mask = z_size - 1;

      for (uint i = 0; i < io_re.size(); ++i) {
        uint pos  = (z_size + sample_idx + i) & mask;
        V    y_re = io_re[i];
        V    y_im = io_im[i];
        // calculation
        auto z_re     = z_re_ptr[pos];
        auto z_im     = z_im_ptr[pos];
        z_re_ptr[pos] = y_re;
        z_im_ptr[pos] = y_im;
        auto y_re_cp  = y_re;
        io_re[i]      = (a_v * y_re_cp) - (b_v * y_im) + z_re;
        io_im[i]      = (b_v * y_re_cp) + (a_v * y_im) + z_im;
      }
      // adjusting real buffer position and buffers size for the next stage
      z_re_ptr = z_im_ptr + z_size;
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
  static void reset_coeffs (crange<V> co, vec_complex<V> pole)
  {
    assert (co.size() >= n_coeffs);

    t_rev_cpole::reset_coeffs (co, pole);
    co[ratio] = vec_real (pole) / vec_imag (pole);
  }
  //----------------------------------------------------------------------------
  template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
  static void reset_states (crange<V> st, uint n_stages)
  {
    t_rev_cpole::reset_states (st, n_stages);
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

    vec_complex<V> y = t_rev_cpole::tick (
      co, st, vec_complex<V> (in_re), n_stages, sample_idx);

    return vec_real (y) + co[ratio] * vec_imag (y);
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

    std::array<V, 32> imag;

    while (io.size() > 0) {
      uint blocksize = std::min (io.size(), imag.size());
      memset (&imag[0], 0, blocksize * sizeof imag[0]);
      t_rev_cpole::tick (
        co,
        st,
        make_crange (io.data(), blocksize),
        make_crange (imag),
        n_stages,
        sample_idx);

      for (uint i = 0; i < blocksize; ++i) {
        io[i] += co[ratio] * imag[i];
      }
      sample_idx += blocksize;
      io = io.shrink_head (blocksize);
    }
    return io; // forwarding
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
} // namespace artv
