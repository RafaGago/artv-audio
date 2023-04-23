#pragma once

#include <type_traits>

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/short_ints.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

// not winning a competition with these, but easy code for the intended range.
//------------------------------------------------------------------------------
template <class T>
inline constexpr T primes_table_size_guess (T minv, T maxv)
{
  T n_guess = maxv - minv;
  n_guess   = div_ceil<T> (n_guess, 2); // even numbers
  n_guess   = div_ceil<T> (n_guess * 2, 3); // divisible by three
  return n_guess;
}
//------------------------------------------------------------------------------
template <class T>
inline constexpr xspan<T> make_primes_table (xspan<T> mem, T minv, T maxv)
{
  static_assert (std::is_unsigned_v<T>);

  auto it = mem;
  T    n  = 0;
  auto v  = minv;

  while (it.size() && v < maxv) {
    T prime = 1;
    T i     = 2;
    while (i <= (v / 2) && (prime == 1)) {
      prime = (v % i) != 0;
      i += 1;
    }
    if (prime == 1) {
      it[0] = v;
      it.cut_head (1);
    }
    v += 1;
  }
  return mem.cut_head ((v >= maxv) ? mem.size() - it.size() : 0);
}
//------------------------------------------------------------------------------
} // namespace artv
