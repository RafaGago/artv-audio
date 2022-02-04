#pragma once

#pragma once

#include <numeric>

#include "artv-common/dsp/own/classes/delay_line.hpp"
#include "artv-common/dsp/own/classes/diffusion_matrix.hpp"

#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/bits.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {

// this may not belong here:

//------------------------------------------------------------------------------
static constexpr uint gcd (uint a, uint b)
{
  while (a != 0) {
    int b_prev = b;
    b          = a;
    a          = b_prev % a;
  }
  return b;
}
//------------------------------------------------------------------------------
static constexpr uint eulers_totient (uint x)
{
  uint y = x;

  for (uint p = 2; p * p <= x; ++p) {
    if (x % p == 0) {
      do {
        x /= p;
      } while (x % p == 0);
      y -= y / p;
    }
  }
  return (x > 1) ? (y - y / x) : y;
}
//------------------------------------------------------------------------------
class reverb {
public:
  //----------------------------------------------------------------------------
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;
  //----------------------------------------------------------------------------
  struct time_tag {};
  void set (time_tag, float v)
  {
    if (v == _fb_rt60_sec) {
      return;
    }
    _fb_rt60_sec = v;
    compute_feedback_gain();
  }

  static constexpr auto get_parameter (time_tag)
  {
    return float_param ("sec", 0.1f, 10.9f, 2.f, 0.01f, 0.5f);
  }

  //----------------------------------------------------------------------------
  struct predelay_tag {};
  void set (predelay_tag, float v) {}

  static constexpr auto get_parameter (predelay_tag)
  {
    return float_param ("sec", 0.01f, 2.f, 0.1f, 0.01f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct er_late_tag {};
  void set (er_late_tag, float v) {}

  static constexpr auto get_parameter (er_late_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct algorithm_tag {};
  void set (algorithm_tag, uint v)
  {
    if (v == _diff_mode) {
      return;
    }
    _diff_mode = v;
    reset_diffusor_times();
  }

  static constexpr auto get_parameter (algorithm_tag)
  {
    return choice_param (
      0, make_cstr_array ("1", "2", "3", "4", "5", "6", "7"), 64);
  }
  //----------------------------------------------------------------------------
  struct algorithm_param_tag {};
  void set (algorithm_param_tag, float v) {}

  static constexpr auto get_parameter (algorithm_param_tag)
  {
    return float_param ("%", 0.f, 100.f, 520.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct tail_gap_tag {};
  void set (tail_gap_tag, float v) {}

  static constexpr auto get_parameter (tail_gap_tag)
  {
    return float_param ("sec", 0.01f, 2.f, 0.1f, 0.01f, 0.6f);
  }
  //----------------------------------------------------------------------------
  struct tilt_tag {};
  void set (tilt_tag, float v) {}

  static constexpr auto get_parameter (tilt_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  struct width_tag {};
  void set (width_tag, float v) {}

  static constexpr auto get_parameter (width_tag)
  {
    return float_param ("%", 0.f, 100.f, 20.f, 0.01f);
  }
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    time_tag,
    predelay_tag,
    er_late_tag,
    algorithm_tag,
    algorithm_param_tag,
    tail_gap_tag,
    tilt_tag,
    width_tag>;

  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext  = &pc;
    _diff_mode    = 0;
    _fb_val       = decltype (_fb_val) {};
    _fb_base_gain = 0.3f;
    _fb_time_ms = 5.f; // this is constant as of now, echoes shouldn't be heard
    _fb_base_time = ms_to_spls<float> (_fb_time_ms);
    _fb_rt60_sec  = 0.;

    allocate();
    reset_diffusor_times();
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (crange<T*> outs, crange<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    static uint todo_rm = 0;

    std::array<std::array<float_x1, n_channels>, 32> io;

    for (uint offset = 0; offset < samples; offset += io.size()) {
      uint blocksize = std::min<uint> (io.size(), samples - offset);

      // interleave and create the channels
      for (uint i = 0; i < blocksize; ++i) {
        float_x1 l = vec_set<1> (ins[0][offset + i]);
        float_x1 r = vec_set<1> (ins[1][offset + i]);

        io[i] = make_array (l, r, r, r, r, r, l, l);
      }

      // diffusor stages
      for (uint df = 0; df < n_diffusors; ++df) {
        // delay lines
        for (uint i = 0; i < blocksize; ++i) {
          _diff[df].push (io[i]);
          for (uint ch = 0; ch < n_channels; ++ch) {
            float t   = _diff_base_times[df][ch];
            io[i][ch] = _diff[df].get<linear_interp> (t, ch);
          }
        }
#if 1
        // shuffle
        for (uint i = 0; i < blocksize; ++i) {
          auto io0 = io[i][0];
          for (uint c = 0; c < n_channels - 1; ++c) {
            io[i][c] = io[i][c + 1] * ((c & 1) ? 1.f : -1.f);
          }
          io[i].back() = io0;
        }
#endif
        // diffusion matrix
        for (uint i = 0; i < blocksize; ++i) {
          io[i] = hadamard_matrix<n_channels>::tick (io[i]);
        }
      }
#if 0
      // feedback loop
      for (uint i = 0; i < blocksize; ++i) {
        for (uint c = 0; c < n_channels; ++c) {
          io[i][c] += _fb_base_gain * _fb_val[c];
        }
        _fb_delay.push (io[i]);
        for (uint ch = 0; ch < n_channels; ++ch) {
          io[i][ch] = _fb_delay.get (_fb_base_time, ch);
        }
        _fb_val = householder_matrix<n_channels>::tick (io[i]);
      }
#endif

      if ((todo_rm & ((1 << 12) - 1)) == 0) {
        printf ("iter %u\n", todo_rm);
        printf (
          "%f %f %f %f\n", io[0][0][0], io[0][1][0], io[0][2][0], io[0][3][0]);
        printf (
          "%f %f %f %f\n", io[0][4][0], io[0][5][0], io[0][6][0], io[0][7][0]);
      }

      // deinterleave and output
      for (uint i = 0; i < blocksize; ++i) {
        outs[0][offset + i] = io[i][0][0];
        outs[1][offset + i] = io[i][1][0];
      }
      ++todo_rm;
    }
  }
  //----------------------------------------------------------------------------
private:
  static constexpr bool  interleaved         = true;
  static constexpr float max_lfo_depth_up_ms = 5.f;
  static constexpr float max_fb_delay_ms     = 2000.f;
  static constexpr uint  n_channels          = 8;
  static constexpr uint  n_diffusors         = 5;
  static constexpr uint  n_diffusor_modes    = 7;
  //----------------------------------------------------------------------------
  template <class T>
  T ms_to_spls (float ms)
  {
    float srate_spl_msec = (float) _plugcontext->get_sample_rate() * 0.001f;
    if constexpr (std::is_integral_v<T>) {
      return (T) std::ceil (ms * srate_spl_msec);
    }
    else {
      return ms * srate_spl_msec;
    }
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  std::array<T, N> ms_to_spls (std::array<float, N> ms)
  {
    std::array<T, N> ret;
    for (uint i = 0; i < N; ++i) {
      ret[i] = ms_to_spls<T> (ms[i]);
    }
    return ret;
  }
  //----------------------------------------------------------------------------
  template <size_t N>
  std::array<uint, N> to_next_pow2 (std::array<uint, N> in)
  {
    for (auto& v : in) {
      v = pow2_round_ceil (v);
    }
    return in;
  }
  //----------------------------------------------------------------------------
  static constexpr std::array<float, n_diffusors> get_max_diffusor_delays_ms()
  {
    std::array<float, n_diffusors> maxv {};

    for (uint m = 0; m < n_diffusor_modes; ++m) {
      auto diff = get_diffusor_times (m);
      for (uint i = 0; i < diff.size(); ++i) {
        for (auto v : diff[i]) {
          maxv[i] = std::max (maxv[i], v + max_lfo_depth_up_ms);
        }
      }
    }
    return maxv;
  }
  //----------------------------------------------------------------------------
  static constexpr std::array<std::array<float, n_channels>, n_diffusors>
  get_diffusor_times (uint mode)
  {
    switch (mode) {
    default:
      [[fallthrough]];
    case 0:
      return make_array (
        get_diff_array (0, 90, 0),
        get_diff_array (180, 90, 8),
        get_diff_array (180, 360, 16),
        get_diff_array (720, 360, 24),
        get_diff_array (720, 1440, 32));
    case 1:
      return make_array (
        get_diff_array (60, 0, 0),
        get_diff_array (60, 120, 8),
        get_diff_array (120, 240, 16),
        get_diff_array (240, 480, 24),
        get_diff_array (960, 480, 32));
    case 2:
      return make_array (
        get_diff_array (0, 50, 0),
        get_diff_array (50, 100, 8),
        get_diff_array (200, 100, 16),
        get_diff_array (400, 200, 24),
        get_diff_array (400, 800, 32));
    case 3:
      return make_array (
        get_diff_array (15, 0, 0),
        get_diff_array (15, 30, 8),
        get_diff_array (30, 60, 16),
        get_diff_array (60, 120, 24),
        get_diff_array (240, 120, 32));
    case 4:
      return make_array (
        get_diff_array (0, 7, 0),
        get_diff_array (14, 7, 8),
        get_diff_array (28, 14, 16),
        get_diff_array (28, 56, 24),
        get_diff_array (112, 56, 32));
    case 5:
      return make_array (
        get_diff_array (0.8, 0, 0),
        get_diff_array (1.6, 0.8, 8),
        get_diff_array (1.6, 3.2, 16),
        get_diff_array (6.4, 3.2, 24),
        get_diff_array (6.4, 12.8, 32));
    case 6:
      return make_array (
        get_diff_array (0.25, 0, 0),
        get_diff_array (0.5, 0.25, 8),
        get_diff_array (0.5, 1.0, 16),
        get_diff_array (1.0, 2.0, 24),
        get_diff_array (4.0, 2.0, 32));
      // adding or removing her requires adjusting "n_diffusor_modes" to match.
    }
  }
  //----------------------------------------------------------------------------
  static constexpr std::array<float, n_channels> get_diff_array (
    float range_start,
    float range_end,
    uint  rnd_start)
  {
    // Until the RNG's are constexpr if provided with a seed.
    // does this belong here?
    // for x in range(32):
    //    random.random()
    constexpr auto rnd = make_array<double> (
      0.9171434810105141,
      0.8569858412166442,
      0.5178699413011407,
      0.8658419727056448,
      0.09615608560228828,
      0.8657091878698523,
      0.8569333970393207,
      0.3780605117952399,
      0.26031208092491054,
      0.5635124119976632,
      0.9790658438505838,
      0.8562823856318246,
      0.21556298702180277,
      0.8600632971753791,
      0.662714633786504,
      0.2225621933588111,
      0.6457530747930535,
      0.7827105700278855,
      0.6705869303441022,
      0.5154710337106151,
      0.815454332575039,
      0.6179902227520485,
      0.7115313466684177,
      0.9378033055153567,
      0.21433529585823752,
      0.8701474992411431,
      0.7086038807361402,
      0.30052303721084295,
      0.28393219786694557,
      0.5983530311667046,
      0.20020536916058207,
      0.6392286472751323,
      0.37143886775293566,
      0.6898805855917455,
      0.1884387811019529,
      0.5686068227042015,
      0.9620012698662417,
      0.4707056753390745,
      0.5257648252025556,
      0.6742146878570825,
      0.7082473720416732,
      0.13154771079490413,
      0.881639016990223,
      0.5319574884475743,
      0.37221870621745656,
      0.29767888275867493,
      0.7901841537170252,
      0.9446750496773592,
      0.9225501410767506,
      0.9805160330674125,
      0.37215404486327974,
      0.8896940430361793,
      0.4460397289641458,
      0.5596925008309813,
      0.972865691753777,
      0.09152757909635534,
      0.8255157060110575,
      0.5475708401593774,
      0.8558832930841987,
      0.944978975726182,
      0.265799720765052,
      0.7421827049142025,
      0.7365250938202301,
      0.9102159605707534);

    auto rnd_factors = make_array<float> (
      (float) rnd[(rnd_start + 0) % rnd.size()],
      (float) rnd[(rnd_start + 1) % rnd.size()],
      (float) rnd[(rnd_start + 2) % rnd.size()],
      (float) rnd[(rnd_start + 3) % rnd.size()],
      (float) rnd[(rnd_start + 4) % rnd.size()],
      (float) rnd[(rnd_start + 5) % rnd.size()],
      (float) rnd[(rnd_start + 6) % rnd.size()],
      (float) rnd[(rnd_start + 7) % rnd.size()]);

    float step = (range_end - range_start) / (float) n_channels;

    return make_array (
      range_start + (0 * step) + rnd_factors[0],
      range_start + (1 * step) + rnd_factors[1],
      range_start + (2 * step) + rnd_factors[2],
      range_start + (3 * step) + rnd_factors[3],
      range_start + (4 * step) + rnd_factors[4],
      range_start + (5 * step) + rnd_factors[5],
      range_start + (6 * step) + rnd_factors[6],
      range_start + (7 * step) + rnd_factors[7]);
  }
  //----------------------------------------------------------------------------
  void allocate()
  {
    auto diff_max    = ms_to_spls<uint> (get_max_diffusor_delays_ms());
    diff_max         = to_next_pow2 (diff_max);
    uint max_spls    = std::accumulate (diff_max.begin(), diff_max.end(), 0);
    uint fb_max_spls = pow2_round_ceil (ms_to_spls<uint> (max_fb_delay_ms));
    max_spls += fb_max_spls;

    _mem.clear();
    _mem.resize (max_spls * n_channels);
    auto ptr = _mem.data();

    for (uint i = 0; i < diff_max.size(); ++i) {
      uint size = diff_max[i] * n_channels;
      _diff[i].reset (make_crange (ptr, size), n_channels);
      ptr += size;
    }

    uint size = fb_max_spls * n_channels;
    _fb_delay.reset (make_crange (ptr, size), n_channels);
    ptr += size;

    assert (ptr == (&_mem.back()) + 1);
  }
  //----------------------------------------------------------------------------
  void reset_diffusor_times()
  {
    _diff_base_times = get_diffusor_times (_diff_mode);
    for (auto& arr : _diff_base_times) {
      arr = ms_to_spls<float> (arr);
    }
  }
  //----------------------------------------------------------------------------
  void compute_feedback_gain()
  {
    // TODO: this will need to happen for each delay line
    float rate    = (float) _plugcontext->get_sample_rate() / _fb_base_time;
    _fb_base_gain = exp (M_LN10 * -60. / 20. * 1. / (_fb_rt60_sec * rate));
  }
  //----------------------------------------------------------------------------
  using reverb_delay_line = modulable_delay_line<float_x1, interleaved>;

  std::array<reverb_delay_line, n_diffusors>             _diff;
  std::array<std::array<float, n_channels>, n_diffusors> _diff_base_times;
  uint                                                   _diff_mode;

  reverb_delay_line                _fb_delay;
  float                            _fb_base_time;
  float                            _fb_rt60_sec;
  float                            _fb_time_ms;
  float                            _fb_base_gain;
  std::array<float_x1, n_channels> _fb_val;

  std::vector<float_x1> _mem;

  plugin_context* _plugcontext;
};
//------------------------------------------------------------------------------
} // namespace artv
