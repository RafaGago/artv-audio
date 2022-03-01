#pragma once

#define register
#include <WDL/MersenneTwister.h>
#undef register

namespace artv { namespace jsfx_engine {
namespace impl {
extern MTRand jsfx_mersenne_twister;
}
// using just rand() or a float white noise generator wasn't giving the same
// results on some JSFX (e.g Liteon's "nonlinear"), so this uses WDL's mersene
// twister implementation. It is very likely that this is what JSFX uses under
// the hood.
inline double rand (double max = 1.)
{
  return impl::jsfx_mersenne_twister.rand53() * max;
}

}} // namespace artv::jsfx_engine
