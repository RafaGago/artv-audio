#pragma once

#pragma once

#include <cassert>
#include <limits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/vec_math.hpp"
#include "artv-common/misc/xspan.hpp"

#include "artv-common/dsp/own/parts/misc/envelope.hpp"
#include "artv-common/dsp/own/parts/parts_to_class.hpp"

namespace artv {
// ducker, for e.g. usage on reverbs/echoes
//------------------------------------------------------------------------------
template <class V, enable_if_vec_of_float_point_t<V>* = nullptr>
class ducker {
public:
  //----------------------------------------------------------------------------
  void set_speed (V factor, vec_value_type_t<V> t_spl)
  {
    using T = vec_value_type_t<V>;
    // 1us to 600.1ms
    V t = (T) 0.0001 + (T) 0.6 * factor;
    _env.template reset_coeffs<ducker_idx> (t, (T) t_spl);
    t = (T) 0.0001 + (T) 0.025 * ((T) 1. - factor);
    _env.template reset_coeffs<smooth_idx> (t, (T) t_spl);
  }
  //----------------------------------------------------------------------------
  void set_threshold (V db) { _threshold_lin = vec_db_to_gain (db); }
  //----------------------------------------------------------------------------
  // returns a gain
  V tick (V in)
  {
    using T          = vec_value_type_t<V>;
    constexpr T tiny = -(std::numeric_limits<T>::min() * 1e6f);

    in *= in;
    auto env = _env.template tick<ducker_idx> (in);
    auto rms = vec_sqrt (vec_max (tiny, env));
    rms      = vec_max (tiny, rms);
    auto gr = rms < _threshold_lin ? vec_set<V> ((T) 1.) : _threshold_lin / rms;
    // TODO: use this for smoothing?
    // https://cytomic.com/files/dsp/DynamicSmoothing.pdf
    return _env.template tick<smooth_idx> (gr);
  }
  //----------------------------------------------------------------------------
private:
  enum { ducker_idx, smooth_idx };
  V                                                   _threshold_lin {};
  part_classes<mp_list<envelope, envelope>, V, false> _env;
};
//------------------------------------------------------------------------------
} // namespace artv
