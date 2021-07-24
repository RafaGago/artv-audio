#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "artv-common/dsp/own/blocks/waveshapers/adaa.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/util.hpp"

namespace artv {
// clang-format off
/*
There are many approximations out there suitable for a waveshaper, but is there
one that is integrable analitically? Can the integral be approximated instead?

https://varietyofsound.wordpress.com/2011/02/14/efficient-tanh-computation-using-lamberts-continued-fraction/
(x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x) / (28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135)

https://www.kvraudio.com/forum/viewtopic.php?f=33&t=388650&start=45

inline double vox_fasttanh2( const double x )
{
	const double ax = fabs( x );
	const double x2 = x * x;

	return( x * ( 2.45550750702956 + 2.45550750702956 * ax + ( 0.893229853513558 + 0.821226666969744 * ax ) * x2 ) / ( 2.44506634652299 + ( 2.44506634652299 + x2 ) * fabs( x + 0.814642734961073 * x * ax )));
}

Octave reminder:

x = sym('x') # symbolic
fn = x^7 + 378 * x^5 + 17325 * x^3 + 135135 *x
fd = 28 * x^6 + 3150 * x^4 + 62370* x^2 + 135135
f = fn / fd
ccode (horner (fn, x))
ccode (horner (fd, x))

At the end this was unused.


Tan based on exp?

tanh = exp(x) - exp(-x) / exp(x) + exp (-x)
cosh = 0.5 * (exp(x) + exp(-x))

Notice that:
exp(-x) = 1/exp(x)

So for ADAA it is needed:

-exp (x)
-exp (-x)
-exp (x1)
-exp (-x1)
-exp (x1*x*0.5).
-exp (-(x1*x*0.5)).

exp(x1) and exp(-x1) are stored from the previous round.

So 4 "exp" calls are needed, alternatively 2 exp calls and 2 divisions.
*/

// clang-format on
struct tanh_waveshaper_functions {
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T fn (T x)
  {
    return tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> fn (simd_batch<T, N> x)
  {
    return xsimd::tanh (x);
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int_fn (T x)
  {
#if 0
    return log (cosh (x));
#else
    static constexpr T log2 = (T) 0.69314718056;
    return log (exp (x) + exp (-x)) - log2;
#endif
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int_fn (simd_batch<T, N> x)
  {
    return xsimd::log (xsimd::cosh (x));
  }
  //----------------------------------------------------------------------------
  template <class T, std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
  static T int2_fn (T x)
  {
    // not solveable analitically in a practical/easy way.
    static_assert (!std::is_same_v<T, T>, "TBI?");
    return (T) 0.;
  }
  //----------------------------------------------------------------------------
  template <class T, size_t N>
  static simd_batch<T, N> int2_fn (simd_batch<T, N> x)
  {
    // not solveable analitically in a practical/easy way.
    static_assert (!std::is_same_v<T, T>, "TBI?");
    return {(T) 0.};
  }
  //----------------------------------------------------------------------------
};
//------------------------------------------------------------------------------
template <uint order>
using tanh_waveshaper_adaa = adaa::waveshaper<tanh_waveshaper_functions, order>;
//------------------------------------------------------------------------------
} // namespace artv
