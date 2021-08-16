#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <gtest/gtest.h>

#include "artv-common/misc/mp11.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {
//------------------------------------------------------------------------------
TEST (simd, sanity_check)
{
  vec<int, 4> a {0, 1, 2, 3};
  vec<int, 4> b {3, 2, 1, 0};
  auto        r = a + b;
  ASSERT_EQ (r[0], 3);
  ASSERT_EQ (r[1], 3);
  ASSERT_EQ (r[2], 3);
  ASSERT_EQ (r[3], 3);
}
//------------------------------------------------------------------------------
TEST (simd, broadcast)
{
  using T = vec<int, 4>;
  auto r  = vec_set<T> (3);
  ASSERT_EQ (r[0], 3);
  ASSERT_EQ (r[1], 3);
  ASSERT_EQ (r[2], 3);
  ASSERT_EQ (r[3], 3);
}
//------------------------------------------------------------------------------
TEST (simd, unaligned_load_store)
{
  using array_type = std::array<u64, 2>;
  using vec_type   = vec<u64, 2>;

  // enforce unalignment of the src and dst memory blocks.
  u8   pool[(sizeof (array_type) * 2) + sizeof (array_type::value_type) - 1];
  uint memstart = 0;

  constexpr uint mask = sizeof (array_type) - 1;
  while ((reinterpret_cast<uintptr_t> (&pool[memstart]) & mask)
         != sizeof (array_type::value_type)) {
    ++memstart;
  }
  array_type* src = new (&pool[memstart]) array_type {};
  array_type* dst = new (&pool[memstart + sizeof *src]) array_type {};

  for (auto& v : *src) {
    v = (array_type::value_type) rand();
  }
  auto a = vec_load_unaligned<vec_type> (*src);
  vec_store_unaligned (*dst, a);

  ASSERT_EQ (*src, *dst);
}
//------------------------------------------------------------------------------
TEST (simd, aligned_load_store)
{
  using array_type = std::array<u64, 2>;
  using vec_type   = vec<u64, 2>;

  alignas (vec_type) array_type src, dst;
  for (auto& v : src) {
    v = (array_type::value_type) rand();
  }
  auto a = vec_load<vec_type> (src);
  vec_store (dst, a);

  ASSERT_EQ (src, dst);
}
//------------------------------------------------------------------------------
TEST (simd, cast)
{
  vec<float, 4> src = {1.4, 2.4, 3.4, 4.4};
  auto          dst = vec_cast<int> (src);
  ASSERT_EQ (dst[0], 1);
  ASSERT_EQ (dst[1], 2);
  ASSERT_EQ (dst[2], 3);
  ASSERT_EQ (dst[3], 4);
}
//------------------------------------------------------------------------------
TEST (simd, shuffle)
{
  vec<int, 4> src  = {1, 2, 3, 4};
  vec<int, 4> src2 = {5, 6, 7, 8};

  auto dst = vec_shuffle (src, src2, 0, 1, 4, 5);

  ASSERT_EQ (dst[0], 1);
  ASSERT_EQ (dst[1], 2);
  ASSERT_EQ (dst[2], 5);
  ASSERT_EQ (dst[3], 6);
}
//------------------------------------------------------------------------------
// one of the functions having 1 parameter
TEST (simd, exp)
{
  mp11::mp_for_each<mp_list<float, double>> ([] (auto typeval) {
    using T       = decltype (typeval);
    using vectype = vec<T, 16 / sizeof (T)>;

    T    scalar = (T) 0.12345678901234567890;
    auto vect   = vec_set<vectype> (scalar);

    scalar = exp (scalar);
    vect   = vec_exp (vect);

    ASSERT_NEAR (scalar, vect[0], (T) 0.0000001);
  });
}
//------------------------------------------------------------------------------
// one of the functions having 2 parameters
TEST (simd, pow)
{
  mp11::mp_for_each<mp_list<float, double>> ([] (auto typeval) {
    using T       = decltype (typeval);
    using vectype = vec<T, 16 / sizeof (T)>;

    T    scalar = (T) 0.12345678901234567890;
    auto vect   = vec_set<vectype> (scalar);

    scalar = pow (scalar, scalar);
    vect   = vec_pow (vect, vect);

    ASSERT_NEAR (scalar, vect[0], (T) 0.0000001);
  });
}
//------------------------------------------------------------------------------
} // namespace artv
