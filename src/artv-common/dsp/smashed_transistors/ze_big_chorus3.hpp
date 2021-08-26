#pragma once

// This one has a lot of manual work, as namespace parameters aren't supported
// by the generator. Just doing it because I liked this one a lot.

#include <array>
#include <cmath>
#include <vector>

#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace smashed_transistors {

class ze_big_chorus3 {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type = dsp_types::modulation;
  //----------------------------------------------------------------------------
  struct sl_algo_tag {};

  void set (sl_algo_tag, int v)
  {
    if (v == sl_algo) {
      return;
    }
    sl_algo = v;
    slider();
  }

  static constexpr auto get_parameter (sl_algo_tag)
  {
    // Original slider line: slider1:sl_algo=4<0,11,1{Series,Series of nested
    // I,Series of nested II,Parallel,Parallel series,Series in Parallel,Nested
    // in parallel I,Nested in parallel II,All nested,Series and nest,Series of
    // combs,Catacombs}>Algo
    return choice_param (
      4,
      make_cstr_array (
        "Series",
        "Series Nested 1",
        "Series Nested 2",
        "Parallel",
        "Parallel Series",
        "Series Parallel",
        "Nested Parallel1",
        "Nested Parallel2",
        "All Nested",
        "Series And Nest",
        "Comb Series",
        "Catacombs"));
  }
  //----------------------------------------------------------------------------
  struct sl_delay_tag {};

  void set (sl_delay_tag, float v)
  {
    if (v == sl_delay) {
      return;
    }
    sl_delay = v;
    slider();
  }

  static constexpr auto get_parameter (sl_delay_tag)
  {
    // Original slider line: slider3:sl_delay=15<0,100,0.00001>Static Delay (ms)
    return float_param ("ms", 0.0, 100.0, 15.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_dmod_tag {};

  void set (sl_dmod_tag, float v)
  {
    if (v == sl_dmod) {
      return;
    }
    sl_dmod = v;
    slider();
  }

  static constexpr auto get_parameter (sl_dmod_tag)
  {
    // Original slider line: slider4:sl_dmod=90<0,100,0.00001>Modulated Depth
    // (%)
    return float_param ("%", 0.0, 100.0, 90.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_dispstatic_tag {};

  void set (sl_dispstatic_tag, float v)
  {
    if (v == sl_dispstatic) {
      return;
    }
    sl_dispstatic = v;
    slider();
  }

  static constexpr auto get_parameter (sl_dispstatic_tag)
  {
    // Original slider line: slider6:sl_dispStatic=1<0,1,0.00001>Static Delay
    // disp.
    return float_param ("", 0.0, 1.0, 1.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_dispmod_tag {};

  void set (sl_dispmod_tag, float v)
  {
    if (v == sl_dispmod) {
      return;
    }
    sl_dispmod = v;
    slider();
  }

  static constexpr auto get_parameter (sl_dispmod_tag)
  {
    // Original slider line: slider7:sl_dispMod=1<0,1,0.00001>Modulated Delay
    // disp.
    return float_param ("", 0.0, 1.0, 1.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_rate_tag {};

  void set (sl_rate_tag, float v)
  {
    if (v == sl_rate) {
      return;
    }
    sl_rate = v;
    slider();
  }

  static constexpr auto get_parameter (sl_rate_tag)
  {
    // Original slider line: slider9:sl_rate=0.05<0,1,0.00001>Modulation Rate
    return float_param ("", 0.0, 1.0, 0.05, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_ratedisp_tag {};

  void set (sl_ratedisp_tag, float v)
  {
    if (v == sl_ratedisp) {
      return;
    }
    sl_ratedisp = v;
    slider();
  }

  static constexpr auto get_parameter (sl_ratedisp_tag)
  {
    // Original slider line: slider10:sl_ratedisp=0.1<0,1,0.00001>Rate disp.
    return float_param ("", 0.0, 1.0, 0.1, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_grate_tag {};

  void set (sl_grate_tag, float v)
  {
    if (v == sl_grate) {
      return;
    }
    sl_grate = v;
    slider();
  }

  static constexpr auto get_parameter (sl_grate_tag)
  {
    // Original slider line: slider12:sl_grate=0.1<0,6,0.00001>G ---- Mod Rate
    return float_param ("", 0.0, 6.0, 0.1, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_gratedisp_tag {};

  void set (sl_gratedisp_tag, float v)
  {
    if (v == sl_gratedisp) {
      return;
    }
    sl_gratedisp = v;
    slider();
  }

  static constexpr auto get_parameter (sl_gratedisp_tag)
  {
    // Original slider line: slider13:sl_gratedisp=1<0,1,0.00001>Rate disp.
    return float_param ("", 0.0, 1.0, 1.0, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_g1_tag {};

  void set (sl_g1_tag, float v)
  {
    if (v == sl_g1) {
      return;
    }
    sl_g1 = v;
    slider();
  }

  static constexpr auto get_parameter (sl_g1_tag)
  {
    // Original slider line: slider14:sl_g1=0.65<-0.95,0.95,0.00001>G Max
    return float_param ("", -0.95, 0.95, 0.65, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_g0_tag {};

  void set (sl_g0_tag, float v)
  {
    if (v == sl_g0) {
      return;
    }
    sl_g0 = v;
    slider();
  }

  static constexpr auto get_parameter (sl_g0_tag)
  {
    // Original slider line: slider15:sl_g0=0.75<-0.95,0.95,0.00001>G Min
    return float_param ("", -0.95, 0.95, 0.75, 1e-05);
  }
  //----------------------------------------------------------------------------
  struct sl_gain_tag {};

  void set (sl_gain_tag, float v)
  {
    if (v == sl_gain) {
      return;
    }
    sl_gain = v;
    slider();
  }

  static constexpr auto get_parameter (sl_gain_tag)
  {
    // Original slider line: slider18:sl_gain=0<-48,12>Gain (db)
    return float_param ("", -48.0, 12.0, 0.0, 0.1);
  }
  //----------------------------------------------------------------------------
  struct sl_drywet_tag {};

  void set (sl_drywet_tag, float v)
  {
    if (v == sl_drywet) {
      return;
    }
    sl_drywet = v;
    slider();
  }

  static constexpr auto get_parameter (sl_drywet_tag)
  {
    // Original slider line: slider19:sl_drywet=1<0,1>Dry/Wet
    return float_param ("", 0.0, 1.0, 1.0, 0.001);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    sl_algo_tag,
    sl_delay_tag,
    sl_dmod_tag,
    sl_dispstatic_tag,
    sl_dispmod_tag,
    sl_rate_tag,
    sl_ratedisp_tag,
    sl_grate_tag,
    sl_gratedisp_tag,
    sl_g1_tag,
    sl_g0_tag,
    sl_gain_tag,
    sl_drywet_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    jsfx_init_section();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process_block_replacing (std::array<T*, 2> chnls, uint samples)
  {
    if (_asl_algo != sl_algo || std::abs (sl_delay - _asl_delay) > 10.) {
      NAP_clear (_lines[0]);
      NAP_clear (_lines[1]);
      _asl_algo   = sl_algo;
      _coefSmooth = 1.;
    }
    else {
      _coefSmooth = coefOnePole (4.);
    }
    _asl_delay = sl_delay;
    // phase lock of quadrature lfos  A/B  C/D
    _lines[0].line[d_B].d.p = _lines[0].line[d_A].d.p;
    _lines[1].line[d_B].d.p = _lines[1].line[d_A].d.p;
    _lines[0].line[d_D].d.p = _lines[0].line[d_C].d.p;
    _lines[1].line[d_D].d.p = _lines[1].line[d_C].d.p;

    double (ze_big_chorus3::*algo) (delay_lines & dl, double x);

    // this switch is here to let the compiler remove function calls if
    // necessary. Placing it on another function may be interpreted as if
    // "sl_algo" can be changed during the loop, inhibiting inlining. Copying
    // the function pointer/ value could work, but I'm not so sure about how
    // deep the compiler goes with its analysis.

    switch (sl_algo) {
    case 0:
      algo = &ze_big_chorus3::NAP_alg1;
      break;
    case 1:
      algo = &ze_big_chorus3::NAP_alg2;
      break;
    case 2:
      algo = &ze_big_chorus3::NAP_alg3;
      break;
    case 3:
      algo = &ze_big_chorus3::NAP_alg4;
      break;
    case 4:
      algo = &ze_big_chorus3::NAP_alg5;
      break;
    case 5:
      algo = &ze_big_chorus3::NAP_alg6;
      break;
    case 6:
      algo = &ze_big_chorus3::NAP_alg7;
      break;
    case 7:
      algo = &ze_big_chorus3::NAP_alg8;
      break;
    case 8:
      algo = &ze_big_chorus3::NAP_alg9;
      break;
    case 9:
      algo = &ze_big_chorus3::NAP_alg10;
      break;
    case 10:
      algo = &ze_big_chorus3::NAP_alg11;
      break;
    case 11:
      algo = &ze_big_chorus3::NAP_alg12;
      break;
    default:
      algo = nullptr;
      break;
    }

    for (uint i = 0; i < samples; ++i) {
      ++_k;
      if (_k >= KRATE) {
        NAP_kProc (_lines[0]);
        NAP_kProc (_lines[1]);
        _k = 0;
      }

      double spl0 = chnls[0][i];
      double spl1 = chnls[1][i];
      spl0 *= _gain;
      spl1 *= _gain;
      double in0  = spl0;
      double in1  = spl1;
      spl0        = (this->*algo) (_lines[0], spl0);
      spl1        = (this->*algo) (_lines[1], spl1);
      chnls[0][i] = _gWet * softSat (spl0) + _gDry * in0;
      chnls[1][i] = _gWet * softSat (spl1) + _gDry * in1;
    }
  }
  //----------------------------------------------------------------------------
private:
  using heap_type = float;
  //----------------------------------------------------------------------------
  struct delay_line_intern { // find a better name...
    double     dp, dv, p, v, v0, v1mv0, _v;
    heap_type* t;
    heap_type* dt;
  };
  //----------------------------------------------------------------------------
  struct delay_line {
    delay_line_intern d;
    delay_line_intern g;

    double     u;
    double     v;
    double     x;
    heap_type* z;
  };
  //----------------------------------------------------------------------------
  struct delay_lines {
    std::array<delay_line, 7> line;
    uint                      i;
  };
  //----------------------------------------------------------------------------
  struct lfo_def {
    heap_type* t;
    heap_type* dt;
  };
  //----------------------------------------------------------------------------
  double coefOnePole (double fc)
  {
    double c = std::cos (2. * M_PI * fc / _plugcontext->get_sample_rate()) - 2.;
    return 1. + c + std::sqrt (c * c - 1.);
  }
  //----------------------------------------------------------------------------
  void realloc_mem()
  {
    _heap.clear();
    // 0.4 sec max delay (covers max modulation and dispersion)
    _delay_line_length = 1 + (uint) (0.4 * _plugcontext->get_sample_rate());
    uint total_mem     = _delay_line_length;
    total_mem *= _lines.size() * _lines[0].line.size(); // 2 * 7
    total_mem += _lfo.size() * 2 * rlfo_size;
    _heap.resize (total_mem);
    uint pos = 0;

    for (auto& lfo : _lfo) {
      lfo.t = &_heap[pos];
      pos += rlfo_size;
      lfo.dt = &_heap[pos];
      pos += rlfo_size;
    }

    for (auto& side : _lines) {
      for (delay_line& ln : side.line) {
        ln.z = &_heap[pos];
        pos += _delay_line_length;
      }
    }
  }
  //----------------------------------------------------------------------------
  void NAP_init_lr_all()
  {
    for (auto& side : _lines) {
      for (delay_line& ln : side.line) {
        ln.g.v = 0.707;
      }
    }
  }
  //----------------------------------------------------------------------------
  void RLFO_init2 (
    lfo_def& lfo,
    double   h,
    double   phi,
    double   i1,
    double   h1,
    double   i2,
    double   h2,
    double   i3,
    double   h3,
    double   i4,
    double   h4)
  {
    constexpr double dp = 2. * M_PI / rlfo_size;

    double p = 0.;

    for (uint i = 0; i < rlfo_size; ++i) {
      lfo.t[i] = std::sin (
        h * p + phi + i1 * M_PI * std::sin (h1 * p)
        + i2 * M_PI * std::cos (h2 * p) - i3 * M_PI * std::sin (h3 * p)
        - i4 * M_PI * std::cos (h4 * p));
      lfo.t[i] = 0.5 + 0.5 * lfo.t[i];
      p += dp;
    }

    for (uint i = 0; i < rlfo_size; ++i) {
      lfo.dt[i] = lfo.t[(i + 1) % rlfo_size] - lfo.t[i];
    }
  }
  //----------------------------------------------------------------------------
  double RLFO_kProc (delay_line_intern& l)
  {
    l.p += l.dp;
    l.p -= rlfo_size * (l.p >= rlfo_size);
    auto p0 = (u64) l.p;
    return l.t[p0] + (l.p - (double) (p0)) * l.dt[p0];
  }
  //----------------------------------------------------------------------------
  void RLFO_setRate (delay_line_intern& l, double rate)
  {
    rate *= (1. / 16.); // 16 cycles per table
    l.dp = rlfo_size * rate * KRATE / _plugcontext->get_sample_rate();
  }
  //----------------------------------------------------------------------------
  double NAP_interp (heap_type* z, double i)
  {
    // linear interpolation of the delay lines
    auto   i0 = (u64) i;
    double a  = i - (double) i0;
    i0 %= _delay_line_length;
    u64 i1 = (i0 + 1) % _delay_line_length;
    return z[i0] + a * (z[i1] - z[i0]);
  }
  //----------------------------------------------------------------------------
  void NAP_pDelay (delay_line& ln, uint i, double il)
  {
    ln.g._v += ln.g.dv;
    ln.g.v += _coefSmooth * (ln.g._v - ln.g.v);
    ln.d._v += ln.d.dv;
    ln.d.v += _coefSmooth * (ln.d._v - ln.d.v);
    ln.z[i] = ln.u;
    ln.v    = NAP_interp (ln.z, il - ln.d.v);
  }
  //----------------------------------------------------------------------------
  void NAP_pDelays (delay_lines& ln)
  {
    ln.i    = (ln.i + 1) % _delay_line_length;
    uint il = ln.i + _delay_line_length;
    for (delay_line& dl : ln.line) {
      NAP_pDelay (dl, ln.i, il);
    }
  }
  //----------------------------------------------------------------------------
  double softSat (double x)
  {
    x = std::min (1., std::max (-1., x));
    return x * (1.5 - 0.5 * x * x);
  }
  //----------------------------------------------------------------------------
  double p (delay_line& A, double x)
  {
    A.u = x - A.g.v * A.v;
    return A.v + A.g.v * A.u;
  }
  //----------------------------------------------------------------------------
  double p (delay_line& A, delay_line& B, double x)
  {
    A.x = x - B.g.v * B.v;
    A.u = A.x - A.g.v * A.v;

    B.u = A.v + A.g.v * A.u;
    return B.v + B.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  double p (delay_line& A, delay_line& B, delay_line& C, double x)
  {
    A.x = x - C.g.v * C.v;
    A.u = A.x - A.g.v * A.v;

    B.x = A.v + A.g.v * A.u;
    B.u = B.x - B.g.v * B.v;

    C.u = B.v + B.g.v * B.u;
    return C.v + C.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  double p (
    delay_line& A,
    delay_line& B,
    delay_line& C,
    delay_line& D,
    double      x)
  {
    A.x = x - D.g.v * D.v;
    A.u = A.x - A.g.v * A.v;

    B.x = A.v + A.g.v * A.u;
    B.u = B.x - B.g.v * B.v;

    C.x = B.v + B.g.v * B.u;
    C.u = C.x - C.g.v * C.v;

    D.u = C.v + C.g.v * C.u;
    return D.v + D.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  double p (
    delay_line& A,
    delay_line& B,
    delay_line& C,
    delay_line& D,
    delay_line& E,
    double      x)
  {
    A.x = x - E.g.v * E.v;
    A.u = A.x - A.g.v * A.v;

    B.x = A.v + A.g.v * A.u;
    B.u = B.x - B.g.v * B.v;

    C.x = B.v + B.g.v * B.u;
    C.u = C.x - C.g.v * C.v;

    D.x = C.v + C.g.v * C.u;
    D.u = D.x - D.g.v * D.v;

    E.u = D.v + D.g.v * D.u;
    return E.v + E.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  double p (
    delay_line& A,
    delay_line& B,
    delay_line& C,
    delay_line& D,
    delay_line& E,
    delay_line& F,
    double      x)
  {
    A.x = x - F.g.v * F.v;
    A.u = A.x - A.g.v * A.v;

    B.x = A.v + A.g.v * A.u;
    B.u = B.x - B.g.v * B.v;

    C.x = B.v + B.g.v * B.u;
    C.u = C.x - C.g.v * C.v;

    D.x = C.v + C.g.v * C.u;
    D.u = D.x - D.g.v * D.v;

    E.x = D.v + D.g.v * D.u;
    E.u = E.x - E.g.v * E.v;

    F.u = E.v + E.g.v * E.u;
    return F.v + F.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  double p (
    delay_line& A,
    delay_line& B,
    delay_line& C,
    delay_line& D,
    delay_line& E,
    delay_line& F,
    delay_line& G,
    double      x)
  {
    A.x = x - G.g.v * G.v;
    A.u = A.x - A.g.v * A.v;

    B.x = A.v + A.g.v * A.u;
    B.u = B.x - B.g.v * B.v;

    C.x = B.v + B.g.v * B.u;
    C.u = C.x - C.g.v * C.v;

    D.x = C.v + C.g.v * C.u;
    D.u = D.x - D.g.v * D.v;

    E.x = D.v + D.g.v * D.u;
    E.u = E.x - E.g.v * E.v;

    F.x = E.v + E.g.v * E.u;
    F.u = F.x - F.g.v * F.v;

    G.u = F.v + F.g.v * F.u;
    return G.v + G.g.v * A.x;
  }
  //----------------------------------------------------------------------------
  void k (delay_line_intern& dli)
  {
    double _v = dli.v0 + RLFO_kProc (dli) * dli.v1mv0;
    dli.dv    = (_v - dli.v) * _KRATE;
  }
  //----------------------------------------------------------------------------
  // jump
  void j (delay_line_intern& dli)
  {
    dli.v  = dli.v0 + RLFO_kProc (dli) * dli.v1mv0;
    dli.dv = 0.;
  }
  //----------------------------------------------------------------------------
  void NAP_j (delay_lines& dl)
  {
    for (delay_line& ln : dl.line) {
      j (ln.d);
      j (ln.g);
    }
  }
  //----------------------------------------------------------------------------
  void NAP_kProc (delay_lines& dl)
  {
    for (delay_line& ln : dl.line) {
      k (ln.d);
      k (ln.g);
    }
  }
  //----------------------------------------------------------------------------
  void NAP_clear (delay_lines& dl)
  {
    for (delay_line& ln : dl.line) {
      std::memset (ln.z, 0, _delay_line_length);
    }
    NAP_j (dl);
  }
  //----------------------------------------------------------------------------
  //                               Algorithms based on basic structures
  // this = r or l
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 1 [A]~[B*]~[C]~[D]~[E]~[F]~[G]                               series
  enum { d_A, d_B, d_C, d_D, d_E, d_F, d_G };

  double NAP_alg1 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_G],
      p (
        dl.line[d_F],
        p (
          dl.line[d_E],
          p (
            dl.line[d_D],
            p (dl.line[d_C], p (dl.line[d_B], p (dl.line[d_A], x)))))));
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 2 [[A]C*] ~ [[B]D] ~ [[E]F] ~ [G]                  series of nested
  double NAP_alg2 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_G],
      p (
        dl.line[d_E],
        dl.line[d_F],
        p (dl.line[d_B], dl.line[d_D], p (dl.line[d_A], dl.line[d_C], x))));
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 3 [[A][C]E]~[[B][D]F]~[G*]                         series of nested
  double NAP_alg3 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_G],
      p (
        dl.line[d_B],
        dl.line[d_D],
        dl.line[d_F],
        p (dl.line[d_A], dl.line[d_C], dl.line[d_E], x)));
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 4 [A*]-[B*]+[C*]+[D*]+[E*]+[F*]-[G*]                      parallel
  // added some different gains to limit beating effects
  double NAP_alg4 (delay_lines& dl, double x)
  {
    double y = 0.377
      * (p (dl.line[d_A], x) - 0.7 * p (dl.line[d_B], x)
         + 1.4 * p (dl.line[d_C], x) + p (dl.line[d_D], x)
         + 0.5 * p (dl.line[d_E], x) + p (dl.line[d_F], x)
         - 1.5 * p (dl.line[d_G], x));
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 5 ([A*]+[B]-[C*]-[D]+[E*]) ~ [F] ~ [G]         parallel and series
  double NAP_alg5 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_G],
      p (
        dl.line[d_F],
        0.447
          * (p (dl.line[d_A], x) + p (dl.line[d_B], x) - p (dl.line[d_C], x) - p (dl.line[d_D], x) + p (dl.line[d_E], x))));

    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 6 [A]~[C*] + [B]~[D*] - [E]~[F*] + [G*]         series in parallel
  double NAP_alg6 (delay_lines& dl, double x)
  {
    double y = p (dl.line[d_C], p (dl.line[d_A], x));
    y += p (dl.line[d_D], p (dl.line[d_B], x));
    y -= p (dl.line[d_F], p (dl.line[d_E], x));
    y += p (dl.line[d_G], x);
    y *= 0.5;
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 7 [[A]C*] + [[B]D*] - [[E]F*] + [G]             nested in parallel
  double NAP_alg7 (delay_lines& dl, double x)
  {
    double y = p (dl.line[d_A], dl.line[d_C], x);
    y += p (dl.line[d_B], dl.line[d_D], x);
    y -= p (dl.line[d_E], dl.line[d_F], x);
    y += p (dl.line[d_G], x);
    y *= 0.5;
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 8 ([[A]C*] - [[B]D*]) ~ [[E][F] G]              nested in parallel
  double NAP_alg8 (delay_lines& dl, double x)
  {
    double y = p (dl.line[d_A], dl.line[d_C], x);
    y -= p (dl.line[d_B], dl.line[d_D], x);
    y = p (dl.line[d_E], dl.line[d_F], dl.line[d_G], y);
    y *= 0.7;
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 9 [[A][B][C][D][E][F] G*]                               all nested
  double NAP_alg9 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_A],
      dl.line[d_B],
      dl.line[d_C],
      dl.line[d_D],
      dl.line[d_E],
      dl.line[d_F],
      dl.line[d_G],
      x);
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 10 [A] ~ [B*] ~ [[C][D][E][F] G]                   series and nest
  double NAP_alg10 (delay_lines& dl, double x)
  {
    double y = p (
      dl.line[d_C],
      dl.line[d_D],
      dl.line[d_E],
      dl.line[d_F],
      dl.line[d_G],
      p (dl.line[d_B], p (dl.line[d_A], x)));
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // 11 [A] +~ ([B*]+[C*]) ~ ([D]+[E]) ~ ([F]+[G])      series of combs
  double NAP_alg11 (delay_lines& dl, double x)
  {
    double y = x + p (dl.line[d_A], x);
    y        = p (dl.line[d_B], y) + p (dl.line[d_C], y);
    y        = p (dl.line[d_D], y) + p (dl.line[d_E], y);
    y        = p (dl.line[d_F], y) + p (dl.line[d_G], y);
    y *= 1. / sqrt (16.);
    NAP_pDelays (dl);
    return y;
  }
  // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
  // _ 12 [A]+~[B*]+~[C]+~[D]+~[E]+~[F]+~[G] Catacombs
  double NAP_alg12 (delay_lines& dl, double x)
  {
    double y = x + p (dl.line[d_A], x);
    y += p (dl.line[d_B], y);
    y += p (dl.line[d_C], y);
    y += p (dl.line[d_D], y);
    y += p (dl.line[d_E], y);
    y += p (dl.line[d_F], y);
    y += p (dl.line[d_G], y);
    y *= 1. / sqrt (128.);
    NAP_pDelays (dl);
    return y;
  }
  //----------------------------------------------------------------------------
  // Assuming all SCP functions are SCOPE related.
  //----------------------------------------------------------------------------
  void copyT (delay_line_intern& dl, lfo_def const& lfo)
  {
    dl.t  = lfo.t;
    dl.dt = lfo.dt;
  }
  //----------------------------------------------------------------------------
  void jsfx_init_section()
  {
    _asl_algo  = -1;
    _asl_delay = 15.;
    _k         = 0;

    realloc_mem();

    RLFO_init2 (_lfo[0], 16, 0.0, 0.5, 1, 0.05, 38, 0.02, 71, 0.01, 113);
    RLFO_init2 (_lfo[1], 16, 1.5, 0.5, 1, 0.05, 39, 0.02, 70, 0.01, 110);
    RLFO_init2 (_lfo[2], 16, 0.0, 0.5, 2, 0.05, 36, 0.02, 73, 0.01, 123);
    RLFO_init2 (_lfo[3], 16, 1.5, 0.5, 2, 0.05, 37, 0.02, 72, 0.01, 113);

    // lines 0 = L. Done with regex.
    copyT (_lines[0].line[d_A].d, _lfo[0]);
    copyT (_lines[0].line[d_B].d, _lfo[1]);
    copyT (_lines[1].line[d_A].d, _lfo[2]);
    copyT (_lines[1].line[d_B].d, _lfo[3]);
    copyT (_lines[0].line[d_C].d, _lfo[2]);
    copyT (_lines[0].line[d_D].d, _lfo[3]);
    copyT (_lines[1].line[d_C].d, _lfo[0]);
    copyT (_lines[1].line[d_D].d, _lfo[1]);
    copyT (_lines[0].line[d_E].d, _lfo[0]);
    copyT (_lines[1].line[d_E].d, _lfo[2]);
    copyT (_lines[0].line[d_F].d, _lfo[1]);
    copyT (_lines[1].line[d_F].d, _lfo[3]);
    copyT (_lines[0].line[d_G].d, _lfo[2]);
    copyT (_lines[1].line[d_G].d, _lfo[0]);
    copyT (_lines[0].line[d_A].g, _lfo[3]);
    copyT (_lines[0].line[d_B].g, _lfo[2]);
    copyT (_lines[0].line[d_C].g, _lfo[1]);
    copyT (_lines[0].line[d_D].g, _lfo[0]);
    copyT (_lines[0].line[d_E].g, _lfo[1]);
    copyT (_lines[0].line[d_F].g, _lfo[2]);
    copyT (_lines[0].line[d_G].g, _lfo[3]);
    copyT (_lines[1].line[d_A].g, _lfo[1]);
    copyT (_lines[1].line[d_B].g, _lfo[0]);
    copyT (_lines[1].line[d_C].g, _lfo[3]);
    copyT (_lines[1].line[d_D].g, _lfo[2]);
    copyT (_lines[1].line[d_E].g, _lfo[3]);
    copyT (_lines[1].line[d_F].g, _lfo[0]);
    copyT (_lines[1].line[d_G].g, _lfo[1]);

    NAP_init_lr_all();

    sl_algo       = get_parameter (sl_algo_tag {}).defaultv;
    sl_delay      = get_parameter (sl_delay_tag {}).defaultv;
    sl_dmod       = get_parameter (sl_dmod_tag {}).defaultv;
    sl_dispstatic = get_parameter (sl_dispstatic_tag {}).defaultv;
    sl_dispmod    = get_parameter (sl_dispmod_tag {}).defaultv;
    sl_rate       = get_parameter (sl_rate_tag {}).defaultv;
    sl_ratedisp   = get_parameter (sl_ratedisp_tag {}).defaultv;
    sl_grate      = get_parameter (sl_grate_tag {}).defaultv;
    sl_gratedisp  = get_parameter (sl_gratedisp_tag {}).defaultv;
    sl_g1         = get_parameter (sl_g1_tag {}).defaultv;
    sl_g0         = get_parameter (sl_g0_tag {}).defaultv;
    sl_gain       = get_parameter (sl_gain_tag {}).defaultv;
    sl_drywet     = get_parameter (sl_drywet_tag {}).defaultv;

    slider();
  }
  //----------------------------------------------------------------------------
  void setGRate (delay_line_intern& l, double disp)
  {
    return RLFO_setRate (l, sl_grate * (1. + sl_gratedisp * disp));
  }
  //----------------------------------------------------------------------------
  void setDelay (delay_line& l, double dispstatic, double dispmod)
  {
    double st = (sl_delay * 0.001 * (double) _plugcontext->get_sample_rate())
      * (1. + sl_dispstatic * dispstatic);
    double md = 0.5 * st * sl_dmod * 0.01 * (1. + sl_dispmod * dispmod);
    l.d.v0    = st - md; // [0 st]
    l.d.v1mv0 = 2. * md; // [st 2*st]
  }
  //----------------------------------------------------------------------------
  void setRate (delay_line_intern& l, double disp)
  {
    return RLFO_setRate (l, sl_rate * (1. + sl_ratedisp * disp));
  }
  //----------------------------------------------------------------------------
  void mulG (delay_line& l, double m)
  {
    l.g.v0 *= m;
    l.g.v1mv0 *= m;
  }
  //----------------------------------------------------------------------------
  void slider()
  {
    //                                    7 x 2 lfos for delay modulation
    setRate (_lines[0].line[d_A].d, 0.60);
    setRate (_lines[1].line[d_A].d, 0.63); // rather used
    setRate (_lines[0].line[d_B].d, 0.60);
    setRate (_lines[1].line[d_B].d, 0.63); // for chorus
    setRate (_lines[0].line[d_C].d, 0.17);
    setRate (_lines[1].line[d_C].d, 0.09);
    setRate (_lines[0].line[d_D].d, 0.17);
    setRate (_lines[1].line[d_D].d, 0.09);
    setRate (_lines[0].line[d_E].d, -0.33);
    setRate (_lines[1].line[d_E].d, -0.35); // rather used
    setRate (_lines[0].line[d_F].d, -0.52);
    setRate (_lines[1].line[d_F].d, -0.49); // for rev
    setRate (_lines[0].line[d_G].d, -0.75);
    setRate (_lines[1].line[d_G].d, -0.76);

    setDelay (_lines[0].line[d_A], -0.70, +0.33);
    setDelay (_lines[1].line[d_A], -0.75, +0.37);
    setDelay (_lines[0].line[d_B], -0.50, +0.05);
    setDelay (_lines[1].line[d_B], -0.57, 0.05);
    setDelay (_lines[0].line[d_C], -0.33, 0.35);
    setDelay (_lines[1].line[d_C], -0.41, 0.30);
    setDelay (_lines[0].line[d_D], 0.00, 0.00);
    setDelay (_lines[1].line[d_D], 0.00, 0.00);
    setDelay (_lines[0].line[d_E], 0.30, -0.27);
    setDelay (_lines[1].line[d_E], 0.32, -0.29); // for rev
    setDelay (_lines[0].line[d_F], 0.63, -0.60);
    setDelay (_lines[1].line[d_F], 0.80, -0.40); // => small variable delay
    setDelay (_lines[0].line[d_G], 0.98, -0.50);
    setDelay (_lines[1].line[d_G], 0.99, -0.60); //    big static delay

    setGRate (_lines[0].line[d_A].g, 0.90);
    setGRate (_lines[1].line[d_A].g, 0.85);
    setGRate (_lines[0].line[d_B].g, 0.35);
    setGRate (_lines[1].line[d_B].g, 0.47);
    setGRate (_lines[0].line[d_C].g, 0.25);
    setGRate (_lines[1].line[d_C].g, 0.17);
    setGRate (_lines[0].line[d_D].g, 0.00);
    setGRate (_lines[1].line[d_D].g, 0.00);
    setGRate (_lines[0].line[d_E].g, -0.13);
    setGRate (_lines[1].line[d_E].g, -0.17);
    setGRate (_lines[0].line[d_F].g, -0.39);
    setGRate (_lines[1].line[d_F].g, -0.33);
    setGRate (_lines[0].line[d_G].g, -0.57);
    setGRate (_lines[1].line[d_G].g, -0.63);
    // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    _lines[0].line[d_A].g.v0     = _lines[1].line[d_A].g.v0
      = _lines[0].line[d_B].g.v0 = _lines[1].line[d_B].g.v0
      = _lines[0].line[d_C].g.v0 = _lines[1].line[d_C].g.v0
      = _lines[0].line[d_D].g.v0 = _lines[1].line[d_D].g.v0
      = _lines[0].line[d_E].g.v0 = _lines[1].line[d_E].g.v0
      = _lines[0].line[d_F].g.v0 = _lines[1].line[d_F].g.v0
      = _lines[0].line[d_G].g.v0 = _lines[1].line[d_G].g.v0 = sl_g0;

    _lines[0].line[d_A].g.v1mv0     = _lines[1].line[d_A].g.v1mv0
      = _lines[0].line[d_B].g.v1mv0 = _lines[1].line[d_B].g.v1mv0
      = _lines[0].line[d_C].g.v1mv0 = _lines[1].line[d_C].g.v1mv0
      = _lines[0].line[d_D].g.v1mv0 = _lines[1].line[d_D].g.v1mv0
      = _lines[0].line[d_E].g.v1mv0 = _lines[1].line[d_E].g.v1mv0
      = _lines[0].line[d_F].g.v1mv0 = _lines[1].line[d_F].g.v1mv0
      = _lines[0].line[d_G].g.v1mv0 = _lines[1].line[d_G].g.v1mv0
      = sl_g1 - sl_g0;
    // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    _gain = (2. / 3.) * pow (2., sl_gain / 6.);

    _gWet = sqrt (sl_drywet);
    _gDry = sqrt (1 - sl_drywet);
    // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    //                                                g depending of algo
    // nested allpass needs lower g to avoid over reverberation
    switch (sl_algo) {
    case 1:
      mulG (_lines[0].line[d_C], 0.7);
      mulG (_lines[1].line[d_C], 0.7);
      mulG (_lines[0].line[d_D], 0.7);
      mulG (_lines[1].line[d_D], 0.7);
      mulG (_lines[0].line[d_F], 0.7);
      mulG (_lines[1].line[d_F], 0.7);
      break;
    case 2:
      mulG (_lines[0].line[d_E], 0.7);
      mulG (_lines[1].line[d_E], 0.7);
      mulG (_lines[0].line[d_F], 0.7);
      mulG (_lines[1].line[d_F], 0.7);
      break;
    case 6:
      mulG (_lines[0].line[d_C], 0.7);
      mulG (_lines[1].line[d_C], 0.7);
      mulG (_lines[0].line[d_D], 0.7);
      mulG (_lines[1].line[d_D], 0.7);
      mulG (_lines[0].line[d_F], 0.7);
      mulG (_lines[1].line[d_F], 0.7);
      break;
    case 7:
      mulG (_lines[0].line[d_C], 0.7);
      mulG (_lines[1].line[d_C], 0.7);
      mulG (_lines[0].line[d_D], 0.7);
      mulG (_lines[1].line[d_D], 0.7);
      mulG (_lines[0].line[d_G], 0.7);
      mulG (_lines[1].line[d_G], 0.7);
      break;
    case 8:
      mulG (_lines[0].line[d_G], 0.7);
      mulG (_lines[1].line[d_G], 0.7);
      break;
    case 9:
      mulG (_lines[0].line[d_G], 0.7);
      mulG (_lines[1].line[d_G], 0.7);
      break;
    case 10:
      mulG (_lines[0].line[d_B], -1.);
      mulG (_lines[1].line[d_B], -1.);
      mulG (_lines[0].line[d_D], -1.);
      mulG (_lines[1].line[d_D], -1.);
      mulG (_lines[0].line[d_F], -1.);
      mulG (_lines[1].line[d_F], -1.);
      break;
    default:
      break;
    }
  }
  //----------------------------------------------------------------------------
  static constexpr uint  rlfo_size = 4096;
  static constexpr uint  KRATE     = 32;
  static constexpr float _KRATE    = 1. / (float) KRATE;
  //----------------------------------------------------------------------------
  // sliders
  uint  sl_algo;
  float sl_delay;
  float sl_dispmod;
  float sl_dispstatic;
  float sl_dmod;
  float sl_drywet;
  float sl_g0;
  float sl_g1;
  float sl_gain;
  float sl_grate;
  float sl_gratedisp;
  float sl_rate;
  float sl_ratedisp;
  //----------------------------------------------------------------------------
  int    _asl_algo;
  double _asl_delay;
  double _coefSmooth;
  uint   _k = 0;
  double _gain;
  double _gWet;
  double _gDry;
  //----------------------------------------------------------------------------
  std::vector<heap_type>     _heap;
  std::array<delay_lines, 2> _lines;
  std::array<lfo_def, 4>     _lfo;
  plugin_context*            _plugcontext;
  uint                       _delay_line_length;
  //----------------------------------------------------------------------------
};
}} // namespace artv::smashed_transistors
