#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

// clang-format off
#include <hiir/PolyphaseIir2Designer.h>
#include <hiir/PolyphaseIir2Designer.cpp>
// clang-format on

#include "artv-common/misc/util.hpp"

// NOTE: Not using <iostream> is deliberate. "fmt" is not available on C++17.
namespace artv {
//------------------------------------------------------------------------------
bool read_double (double& result, char const* str)
{
  char* endptr;
  result = strtod (str, &endptr);
  if (*endptr != 0) {
    fprintf (stderr, "Could not parse as floating point: %s.\n", str);
    return false;
  }
  return true;
}
//------------------------------------------------------------------------------
void hiir_print (
  crange<double> coeffs,
  double         att_db,
  double         transition_width_percent)
{
  printf ("Order: %lu\n", (coeffs.size() * 2) + 1);
  printf ("Passband Att (fB): %lf\n", att_db);
  printf ("Transition Width (%%): %lf\n", transition_width_percent);

  printf ("[] = { \n");
  for (double coeff : coeffs) {
    printf ("  %.17lg, \n", coeff);
  }
  printf ("}\n");

  double tp = 0.25 - (transition_width_percent / (2. * 2. * 100.));
  auto   freqs
    = make_array (0.0005, 0.005, 0.01, 0.05, .1, .15, .2, .24, .25, tp);
  for (auto freq : freqs) {
    double gd = hiir::PolyphaseIir2Designer::compute_group_delay (
      coeffs.data(), coeffs.size(), freq, false);
    float hz_44100 = 44100. * 2. * freq;
    float hz_48000 = 48000. * 2. * freq;

    printf (
      "group delay at (%lg x FS(2x), (%g Hz(44k)), (%g Hz(48k))): %lf\n",
      freq,
      hz_44100,
      hz_48000,
      gd);
  }
}
//------------------------------------------------------------------------------
constexpr uint hiir_max_coeffs
  = (hiir::PolyphaseIir2Designer::_max_order - 1) / 2;
//------------------------------------------------------------------------------
void hiir_calculate (
  double transition_width_percent,
  double order,
  double att_db)
{
  std::array<double, hiir_max_coeffs> coeffs;
  auto transition = transition_width_percent / (2. * 100.); // F / FS, max 0.5

  int n_coeffs;
  if (order != 0.) {
    n_coeffs = ((uint) order - 1) / 2;
    if ((((uint) order) % 2) != 1) {
      puts ("Warning: Orders are always odd. Reducing order");
    }
    att_db = hiir::PolyphaseIir2Designer::compute_atten_from_order_tbw (
      n_coeffs, transition);
    hiir::PolyphaseIir2Designer::compute_coefs_spec_order_tbw (
      coeffs.data(), n_coeffs, transition);
  }
  else {
    n_coeffs = hiir::PolyphaseIir2Designer::compute_coefs (
      coeffs.data(), att_db, transition);
  }

  hiir_print (make_crange (coeffs, n_coeffs), att_db, transition_width_percent);
}
//------------------------------------------------------------------------------
int oversample_util (int argnum, char const** args)
{
  if (argnum == 3 && strcmp (args[0], "hiir-calculate") == 0) {
    double att_db, transition;

    if (read_double (att_db, args[1]) && read_double (transition, args[2])) {
      hiir_calculate (transition, 0., att_db);
      return 0;
    }
    else {
      return 1;
    }
  }
  if (argnum == 3 && strcmp (args[0], "hiir-calculate-best") == 0) {
    double order, transition;

    if (read_double (order, args[1]) && read_double (transition, args[2])) {
      hiir_calculate (transition, order, 0.);
      return 0;
    }
    else {
      return 1;
    }
  }

  puts (R"END(Usage:
 )END");
  return strcmp (args[0], "help") == 0 ? 0 : 1;
}

} // namespace artv
//------------------------------------------------------------------------------
int main (int argc, char const* argv[])
{
  return artv::oversample_util (argc - 1, &argv[1]);
}
//------------------------------------------------------------------------------
