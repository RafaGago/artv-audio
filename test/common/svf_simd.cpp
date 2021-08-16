#include <array>
#include <cstring>
#include <gtest/gtest.h>

#include "artv-common/dsp/own/blocks/filters/andy_svf.hpp"
#include "artv-common/dsp/own/misc.hpp"

namespace artv {

//------------------------------------------------------------------------------
// too lazy to look for a reminder on the GTEST documentation for fixtures.
class svf_test {
public:
  static constexpr uint vecsize = 4;
  using scalar_type             = float;
  using vector_type             = vec<scalar_type, vecsize>;

  static constexpr scalar_type samplerate = 44100;

  std::array<scalar_type, andy::svf::n_coeffs> scalar_coeffs;
  std::array<scalar_type, andy::svf::n_states> scalar_states;
  alignas (vector_type)
    std::array<scalar_type, andy::svf::n_coeffs * vecsize> vector_coeffs;
  alignas (vector_type)
    std::array<scalar_type, andy::svf::n_states * vecsize> vector_states;

  void operator()()
  {
    memset (scalar_states.data(), 0, sizeof scalar_states);
    memset (vector_states.data(), 0, sizeof vector_states);

    for (uint i = 0; i < 2048; ++i) {
      scalar_type in = _whitenoise (1.);
      scalar_type scalar
        = andy::svf::tick<scalar_type> (scalar_coeffs, scalar_states, in);
      auto vect = andy::svf::tick_simd (
        vector_coeffs, vector_states, vec_set<vector_type> (in));
      ASSERT_NEAR (scalar, vect[0], (scalar_type) 0.00001);
    }
  }

private:
  white_noise_generator _whitenoise;
};
//------------------------------------------------------------------------------
TEST (svf, lowpass)
{
  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  andy::svf::lowpass ({test.scalar_coeffs}, f, q, test.samplerate);
  andy::svf::lowpass_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (svf, allpass)
{
  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  andy::svf::allpass ({test.scalar_coeffs}, f, q, test.samplerate);
  andy::svf::allpass_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
TEST (svf, bell)
{
  svf_test              test;
  svf_test::scalar_type f  = 500.;
  svf_test::scalar_type q  = 2.;
  svf_test::scalar_type db = 12.;

  andy::svf::bell ({test.scalar_coeffs}, f, q, db, test.samplerate);
  andy::svf::bell_simd (
    test.vector_coeffs,
    vec_set<svf_test::vector_type> (f),
    vec_set<svf_test::vector_type> (q),
    vec_set<svf_test::vector_type> (db),
    test.samplerate);

  test();
}
//------------------------------------------------------------------------------
} // namespace artv
