#pragma once

// Implementation heavily "inspired" on:
// https://github.com/jatinchowdhury18/ADAA

#include <cmath>

#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv { namespace adaa {

//------------------------------------------------------------------------------
static constexpr double epsilon (double)
{
  return 0.0000000001;
}
static constexpr float epsilon (float)
{
  return 0.00001;
}
//------------------------------------------------------------------------------
template <class functions, uint order>
class waveshaper;
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 1> {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x1_int, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    T x_int   = functions::int_fn (x);
    T x1v     = st[x1];
    T x1v_int = st[x1_int];

    st[x1]     = x;
    st[x1_int] = x_int;

    T diff = x - x1v;
    if (abs (diff) >= epsilon (T {})) {
      return (x_int - x1v_int) / diff;
    }
    else {
      return (functions::fn ((x + x1v) / (T) 2.));
    }
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <class functions>
class waveshaper<functions, 2> {
public:
  enum coeffs { n_coeffs };
  enum state { x1, x2, x2_der, x1_int2, n_states };
  //----------------------------------------------------------------------------
  template <class T>
  static T tick (crange<const T>, crange<T> st, T x)
  {
    using fns = functions;

    T x2v     = st[x2];
    T x2v_der = st[x2_der];
    T x1v_der = get_derivative_x1 (st, x);
    T diff    = x - x2v;
    T ret;

    if (abs (diff) >= epsilon (T {})) {
      ret = ((T) 2. * (x1v_der - x2v_der)) / diff;
    }
    else {
      // fallback
      T a  = (x + x2v) * (T) 0.5;
      diff = a - x;

      if (abs (diff) >= epsilon (T {})) {
        ret = ((T) 2. / diff)
          * (fns::int_fn (a) + (fns::int2_fn (x) - fns::int2_fn (a)) / diff);
      }
      else {
        ret = fns::fn ((a + x) * (T) 0.5);
      }
    }

    st[x2]     = st[x1];
    st[x1]     = x;
    st[x2_der] = x1v_der;

    return ret;
  }
  //----------------------------------------------------------------------------
private:
  //----------------------------------------------------------------------------
  template <class T>
  static T get_derivative_x1 (crange<T> st, T x)
  {
    T x1v      = st[x1];
    T x1v_int2 = st[x1_int2];

    T x_int2 = functions::int2_fn (x);
    T diff   = x - x1v;

    st[x1_int2] = x_int2;

    if (abs (diff) >= epsilon (T {})) {
      return (x_int2 - x1v_int2) / diff;
    }
    else {
      return functions::int_fn ((x + x1v) * (T) 0.5);
    }
  }
  //----------------------------------------------------------------------------
};
}} // namespace artv::adaa
