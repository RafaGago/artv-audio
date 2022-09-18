#pragma once

#include <algorithm>

#include "artv-common/dsp/own/classes/add_ducker.hpp"
#include "artv-common/dsp/own/classes/misc.hpp"
#include "artv-common/dsp/own/classes/pitch_shift.hpp"
#include "artv-common/dsp/own/classes/plugin_context.hpp"
#include "artv-common/dsp/own/classes/reverb_tools.hpp"
#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"
#include "artv-common/dsp/types.hpp"
#include "artv-common/juce/parameter_definitions.hpp"
#include "artv-common/juce/parameter_types.hpp"
#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/simd.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/filters/andy_svf.hpp"

namespace artv { namespace shabtronic {

//------------------------------------------------------------------------------
// this can probably be reused with another name
template <class V, uint N>
class fdn_verb_householder;

//------------------------------------------------------------------------------
template <class V>
class fdn_verb_householder<V, 4> {
public:
  static constexpr uint size = 4;
  //----------------------------------------------------------------------------
  void reset (xspan<V> mem) { _z.reset (mem, size); }
  //----------------------------------------------------------------------------
  void set (std::array<uint, size> t, V feedback)
  {
    set_times (t);
    set_feedback (feedback);
  }
  //----------------------------------------------------------------------------
  void set_times (std::array<uint, size> t)
  {
    for (uint v : t) {
      assert (v < _z.size());
    }
    _t = t;
  }
  //----------------------------------------------------------------------------
  void set_feedback (V feedback) { _fb = feedback; }
  //----------------------------------------------------------------------------
  V tick (V in)
  {
    using T = vec_value_type_t<V>;

    std::array<V, size> x, y;
    V                   factor = vec_set<V> ((T) 0);

    for (uint i = 0; i < size; ++i) {
      x[i] = _z.get (_t[i], i);
      factor += x[i];
    }
    factor *= vec_set<V> ((T) 0.5);
    for (uint i = 0; i < size; ++i) {
      y[i] = in + (x[i] - factor) * _fb;
    }
    _z.push (y);
    return factor;
  }
  //----------------------------------------------------------------------------
private:
  static_delay_line<V, true, true> _z;
  std::array<uint, size>           _t {};
  V                                _fb {};
};
//------------------------------------------------------------------------------
class fdn_verb : private add_ducker<double_x2> {
public:
  static constexpr dsp_types dsp_type  = dsp_types::reverb;
  static constexpr bus_types bus_type  = bus_types::stereo;
  static constexpr uint      n_inputs  = 1;
  static constexpr uint      n_outputs = 1;

public:
  //----------------------------------------------------------------------------
  // Snippet for parameter boilerplate in the authors framework....
  struct time_tag {};

  void set (time_tag, float v)
  {
    decltype (_rtime) time
      = v * 0.001 * (float) _plugcontext->get_sample_rate();
    if (time == _rtime) {
      return;
    }
    _rtime = time;
    set_delay_lengths (_rtime, _cascade, _stereo);
  }

  static constexpr auto get_parameter (time_tag)
  {
    // Original slider line: slider1:500<1,500>-Time
    return float_param ("ms", 1.0, 500.0, 500.0, 0.1);
  }
  //----------------------------------------------------------------------------
  struct feedback_tag {};

  void set (feedback_tag, float v)
  {
    // Original slider line: slider2:0.6<0,1>-Feedback
    decltype (_feedback) fb = v * 0.01f;
    if (fb == _feedback) {
      return;
    }
    _feedback = fb;
    _dlevel   = std::max (fb, 0.75f);
    _dlevel *= _dlevel;
    _dlevel *= _dlevel;
    _dlevel = 1.f / _dlevel;
    _hhmatrix[0].set_feedback (vec_set<1> (_feedback));
    _hhmatrix[1].set_feedback (vec_set<1> (_feedback));
  }

  static constexpr auto get_parameter (feedback_tag)
  {
    // Original slider line: slider2:0.6<0,1>-Feedback
    return float_param ("%", 0.0, 100.0, 60., 0.1);
  }
  //----------------------------------------------------------------------------
  struct density_tag {};

  void set (density_tag, float v)
  {
    decltype (_density) d = v * 0.01;
    if (_density == d) {
      return;
    }
    _density = d;
    for (auto& ap : _ap) {
      ap.set_gain (vec_set<float_x1> (d));
    }
  }

  static constexpr auto get_parameter (density_tag)
  {
    // Original slider line: slider3:0.6<0,1>-Density
    return float_param ("%", 0.0, 100.0, 60., 0.1);
  }
  //----------------------------------------------------------------------------
  struct cascade_tag {};

  void set (cascade_tag, float v)
  {
    decltype (_cascade) casc = v * 0.01;
    if (casc == _cascade) {
      return;
    }
    _cascade = casc;
    set_delay_lengths (_rtime, _cascade, _stereo);
  }

  static constexpr auto get_parameter (cascade_tag)
  {
    // Original slider line: slider4:0.7<0,1>-Cascade
    return float_param ("%", 0.0, 100.0, 70., 0.1);
  }
  //----------------------------------------------------------------------------
  struct mod_rate_tag {};

  void set (mod_rate_tag, float v)
  {
    if (_mod_hz == (decltype (_mod_hz)) v) {
      return;
    }
    _mod_hz = v;
    reset_mod();
  }

  static constexpr auto get_parameter (mod_rate_tag)
  {
    // Original slider line: slider5:0.414<0,10>-ModRate
    return float_param ("Hz", 0.0, 10.0, 0.414, 0.001);
  }
  //----------------------------------------------------------------------------
  struct mod_depth_tag {};

  void set (mod_depth_tag, float v)
  {
    decltype (_mod_depth) depth = v * 0.01;
    if (_mod_depth == depth) {
      return;
    }
    _mod_depth = depth;
    reset_mod();
  }

  static constexpr auto get_parameter (mod_depth_tag)
  {
    // Original slider line: slider6:0.75<0,1>-Depth
    return float_param ("%", 0.0, 100.0, 75., 0.1);
  }
  //----------------------------------------------------------------------------
  struct stereo_tag {};

  void set (stereo_tag, float v)
  {
    decltype (_stereo) st = v * 0.01;
    if (st == _stereo) {
      return;
    }
    _stereo = st;
    set_delay_lengths (_rtime, _cascade, _stereo);
  }

  static constexpr auto get_parameter (stereo_tag)
  {
    // Original slider line: slider7:0.25<-1,1>-Stereo
    return float_param ("%", -100.0, 100.0, 25., 0.1);
  }
  //--------------------------------------------------------------------------
  struct smooth_tag {};

  void set (smooth_tag, float v) { _smooth = v; }

  static constexpr auto get_parameter (smooth_tag)
  {
    // Original slider line: slider8:4<3,32>-Smooth
    return float_param ("", 3.0, 32.0, 4.0, 0.1);
  }
  //--------------------------------------------------------------------------
  struct hipass_tag {};

  void set (hipass_tag, float v)
  {
    decltype (_hp) f = midi_note_to_hz (v);
    if (f == _hp) {
      return;
    }
    _hp = f;
    _filters.reset_coeffs<filter_hp> (
      vec_set<1> (f), vec_set<1> (0.25f), _t_spl);
  }

  static constexpr auto get_parameter (hipass_tag)
  {
    // Original slider line: slider9:10<10,2000>-Hipass
    return frequency_parameter (10.0, 2000.0, 10.0);
  }
  //--------------------------------------------------------------------------
  struct lopass_tag {};

  void set (lopass_tag, float v)
  {
    decltype (_lp) f = midi_note_to_hz (v);
    if (f == _lp) {
      return;
    }
    _lp = f;
    _filters.reset_coeffs<filter_lp> (
      vec_set<1> (f), vec_set<1> (0.25f), _t_spl);
  }

  static constexpr auto get_parameter (lopass_tag)
  {
    // Original slider line: slider10:20000<200,20000>-Lopass
    return frequency_parameter (200.0, 20000.0, 20000.0);
  }
  //--------------------------------------------------------------------------
  struct mix_tag {};

  void set (mix_tag, float v) { _mix = v * 0.01; }

  static constexpr auto get_parameter (mix_tag)
  {
    // Original slider line: slider11:0.5<0,1>-Mix
    return float_param ("%", 0.0, 100.0, 100., 0.1);
  }
  //--------------------------------------------------------------------------
  struct pre_shift_tag {};

  void set (pre_shift_tag, float v)
  {
    if (v == _pre_shift) {
      return;
    }
    _pre_shift = v;
    _fb_shift.set_reader (_fb_shift_read, v);
  }

  static constexpr auto get_parameter (pre_shift_tag)
  {
    // Original slider line: slider12:0<-12,12>-Pitch
    return float_param ("Semitones", -12.0, 12.0, 0.0, 1.);
  }
  //--------------------------------------------------------------------------
  struct post_shift_tag {};

  void set (post_shift_tag, float v)
  {
    if (v == _post_shift) {
      return;
    }
    _post_shift = v;
    _main_shift.set_reader (_main_shift_read, v);
  }

  static constexpr auto get_parameter (post_shift_tag)
  {
    // Original slider line: slider13:0<-12,12>-Post Shift
    return float_param ("Semitones", -12.0, 12.0, 0.0, 1.);
  }
  //----------------------------------------------------------------------------
  using add_ducker::get_parameter;
  using add_ducker::set;
  using ducking_speed_tag     = add_ducker::ducking_speed_tag;
  using ducking_threshold_tag = add_ducker::ducking_threshold_tag;
  //----------------------------------------------------------------------------
  using parameters = mp_list<
    time_tag,
    feedback_tag,
    density_tag,
    cascade_tag,
    mod_rate_tag,
    mod_depth_tag,
    stereo_tag,
    smooth_tag,
    hipass_tag,
    lopass_tag,
    mix_tag,
    pre_shift_tag,
    post_shift_tag,
    ducking_speed_tag,
    ducking_threshold_tag>;
  //----------------------------------------------------------------------------
  void reset (plugin_context& pc)
  {
    _plugcontext = &pc;
    _t_spl       = 1.f / pc.get_sample_rate();
    reinitialize (true);
    _filters.reset_states<filter_lp>();
    _filters.reset_states<filter_hp>();
    _lfo.reset();
    _out_prev  = vec_set<1> (0.f);
    _pre_shift = _post_shift = -999999.;
    add_ducker::reset (pc.get_sample_rate());
    mp11::mp_for_each<parameters> ([=] (auto param) {
      set (param, get_parameter (param).defaultv);
    });
  }
  //----------------------------------------------------------------------------
  template <class T>
  void process (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    add_ducker::process (
      outs,
      ins,
      samples,
      [=] (xspan<T*> outs_fw, xspan<T const*> ins_fw, uint samples_fw) {
        this->process_intern (outs_fw, ins_fw, samples_fw);
      });
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  void process_intern (xspan<T*> outs, xspan<T const*> ins, uint samples)
  {
    assert (outs.size() >= (n_outputs * (uint) bus_type));
    assert (ins.size() >= (n_inputs * (uint) bus_type));

    for (uint i = 0; i < samples; ++i) {
      vec<float, 1> l {ins[0][i]};
      vec<float, 1> r;
      if (unlikely (_feedback >= 1.f)) {
        l = vec_set<1> (0.f);
      }
      l += (_out_prev * _feedback * 0.9f);
      auto old = l;

      l = _mod.get (_mod_delay_spls, _lfo.tick_sine()[0], _mod_depth_spls, 0);
      _mod.push (xspan {&old, 1});
      for (auto& ap : _ap) {
        l = ap.tick (l);
      }
      l = _filters.tick_cascade (l);
      _fb_shift.push (l);
      _out_prev = _fb_shift.read (_fb_shift_read);
      if (unlikely (_feedback >= 1.f)) {
        l = vec_set<1> (0.f);
      }
      r = _hhmatrix[0].tick (l);
      l = _hhmatrix[1].tick (l);
      vec<float, 2> out {l[0], r[0]};
      _main_shift.push (out);
      out        = _main_shift.read (_main_shift_read);
      outs[0][i] = ins[0][i] * (1. - _mix) + out[0] * _dlevel * _mix;
      outs[1][i] = ins[1][i] * (1. - _mix) + out[1] * _dlevel * _mix;
    }
  }
  //----------------------------------------------------------------------------
  void reset_mod()
  {
    float srate      = _plugcontext->get_sample_rate();
    _mod_delay_spls  = srate * (1. / 75.);
    float hz_correct = std::min (1.0 / (2.0 * _mod_hz), 1.);
    _mod_depth_spls  = hz_correct * _mod_depth * _mod_delay_spls * (0.5f);
    _lfo.set_freq (make_vec (_mod_hz), _t_spl);
  }
  //----------------------------------------------------------------------------
  void reinitialize (bool reallocate = false)
  {
    uint delay_samples = _plugcontext->get_sample_rate() / 2; // 500ms
    uint delay_size = 1u << (last_bit_set (delay_samples) - 1); // floor (log2)
    delay_size      = delay_size < delay_samples ? delay_size << 1 : delay_size;
    uint pitch_size = delay_size / 4; // 8192 at 44/48KHz

    constexpr uint n_delays = //
      1 + // _mod
      ap_size + //
      hh_matrix_size * 2;

    constexpr uint n_pitch_buffers = //
      1 + //_fb_shift
      2; // _main_shift

    if (reallocate) {
      _mem.clear();
      _mem.resize (n_delays * delay_size + n_pitch_buffers * pitch_size);
    }
    else {
      xspan_memset (xspan {_mem}, 0);
    }

    auto ptr = _mem.data();

    _mod.reset (xspan {ptr, delay_size}.cast (float_x1 {}), 1);
    ptr += delay_size;

    for (uint i = 0; i < _ap.size(); ++i) {
      _ap[i].reset (xspan {ptr, delay_size}.cast (float_x1 {}));
      ptr += delay_size;
    }

    for (uint i = 0; i < _hhmatrix.size(); ++i) {
      _hhmatrix[i].reset (
        xspan {ptr, delay_size * _hhmatrix[i].size}.cast (float_x1 {}));
      ptr += delay_size * _hhmatrix[i].size;
    }

    _fb_shift.reset (xspan {ptr, pitch_size}.cast (float_x1 {}));
    ptr += pitch_size;

    _main_shift.reset (xspan {ptr, pitch_size * 2}.cast (vec<float, 2> {}));
  }
  //----------------------------------------------------------------------------
  void set_delay_lengths (float length, float warp, float swidth)
  {
    float l1    = length;
    float l2    = length;
    float warp2 = warp;

    std::array<std::array<uint, hh_matrix_size>, 2> hh_lengths {};

    static_assert (ap_size == (2 * hh_matrix_size), "Rewrite the loop");

    for (uint i = 0; i < _ap.size(); ++i) {
      _ap[i].set_time (std::max ((uint) l1, 4u));

      if (i < hh_matrix_size) {
        bool hhidx            = swidth < 0.;
        uint lster            = (uint) (l2 * (1. - std::abs (swidth) * 0.005));
        hh_lengths[hhidx][i]  = std::max ((uint) l2, 4u);
        hh_lengths[!hhidx][i] = std::max ((uint) lster, 4u);
      }
      l1 *= warp;
      l2 *= warp2;
    }
    _hhmatrix[0].set_times (hh_lengths[0]);
    _hhmatrix[1].set_times (hh_lengths[1]);
  }
  //----------------------------------------------------------------------------
  static constexpr uint ap_size        = 8;
  static constexpr uint hh_matrix_size = 4;

  enum filter_type { filter_lp, filter_hp };
  using filter_types = mp_list<andy::svf_lowpass, andy::svf_highpass>;

  lfo<1>                                                     _lfo;
  modulable_delay_line<float_x1, linear_interp, false, true> _mod;

  std::array<allpass_with_params<float_x1>, ap_size>            _ap;
  part_classes<filter_types, float_x1>                          _filters;
  std::array<fdn_verb_householder<float_x1, hh_matrix_size>, 2> _hhmatrix;
  pitch_shift_sin<float_x1, false>                              _fb_shift;
  decltype (_fb_shift)::reader                                  _fb_shift_read;
  pitch_shift_sin<vec<float, 2>, true>                          _main_shift;
  decltype (_main_shift)::reader _main_shift_read;

  std::vector<float_x1> _mem;

  vec<float, 1> _out_prev;

  float _mod_delay_spls;
  float _mod_depth_spls;

  // sliders
  float           _rtime; //      slider1
  float           _feedback; //   slider2
  float           _dlevel;
  float           _density; //    slider3
  float           _cascade; //    slider4
  float           _mod_hz; //     slider5
  float           _mod_depth; //  slider6
  float           _stereo; //     slider7
  float           _smooth; //     slider8
  float           _hp; //         slider9
  float           _lp; //         slider10
  float           _mix; //        slider11
  float           _pre_shift; //  slider12
  float           _post_shift; // slider13
  plugin_context* _plugcontext;
  float           _t_spl;
};
//------------------------------------------------------------------------------

}} // namespace artv::shabtronic
