#pragma once

#include "artv-common/misc/short_ints.hpp"

namespace artv {

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

} // namespace artv
