#include <array>
#include <cstdio>
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

    std::array<scalar_type, 2048> res_scalar, res_vec;

    for (uint i = 0; i < res_scalar.size(); ++i) {
      scalar_type in = _whitenoise (1.);
      res_scalar[i]
        = andy::svf::tick<scalar_type> (scalar_coeffs, scalar_states, in);
      res_vec[i] = andy::svf::tick_simd (
        vector_coeffs, vector_states, vec_set<vector_type> (in))[0];
    }
    // TODO: Do a comparison that prints the failing index. GTEST suggested
    // GMock Matchers and macro stuff I didn't want to lose half an hour with
    // just now.
    for (uint i = 0; i < res_scalar.size(); ++i) {
      constexpr auto epsilon = 0.00001;
      if (abs (res_scalar[i] - res_vec[i]) > epsilon) {
        printf ("On index: %u\n", i);
        ASSERT_NEAR (res_scalar[i], res_vec[i], (scalar_type) epsilon);
      }
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
TEST (svf, highpass)
{
  svf_test              test;
  svf_test::scalar_type f = 500.;
  svf_test::scalar_type q = 2.;

  andy::svf::highpass ({test.scalar_coeffs}, f, q, test.samplerate);
  andy::svf::highpass_simd (
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
