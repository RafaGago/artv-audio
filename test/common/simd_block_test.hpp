#include <array>
#include <cstdio>
#include <cstring>
#include <gtest/gtest.h>

#include "artv-common/dsp/own/misc.hpp"
#include "artv-common/misc/simd.hpp"

namespace artv {

//------------------------------------------------------------------------------
// too lazy to look for a reminder on the GTEST documentation for fixtures.
template <class Dsp>
class simd_block_test {
public:
  static constexpr uint vecsize = 4;
  using scalar_type             = float;
  using vector_type             = vec<scalar_type, vecsize>;

  static constexpr scalar_type samplerate = 44100;

  std::array<scalar_type, Dsp::n_coeffs> scalar_coeffs;
  std::array<scalar_type, Dsp::n_states> scalar_states;
  alignas (
    vector_type) std::array<scalar_type, Dsp::n_coeffs * vecsize> vector_coeffs;
  alignas (
    vector_type) std::array<scalar_type, Dsp::n_states * vecsize> vector_states;

  void operator()()
  {
    memset (scalar_states.data(), 0, sizeof scalar_states);
    memset (vector_states.data(), 0, sizeof vector_states);

    std::array<scalar_type, 2048> res_scalar, res_vec;

    for (uint i = 0; i < res_scalar.size(); ++i) {
      scalar_type in = _whitenoise (1.);
      res_scalar[i]
        = Dsp::template tick<scalar_type> (scalar_coeffs, scalar_states, in);
      res_vec[i] = Dsp::tick_simd (
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

} // namespace artv
