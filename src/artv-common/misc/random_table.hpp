#pragma once

#include "artv-common/misc/misc.hpp"
#include "artv-common/misc/xspan.hpp"

namespace artv {

// a table with values from /dev/random, for predictable "randomness".
namespace detail {

static constexpr auto rnd_table_raw = make_array<double> (
  0.9171434810105141,
  0.8569858412166442,
  0.5178699413011407,
  0.8658419727056448,
  0.09615608560228828,
  0.8657091878698523,
  0.8569333970393207,
  0.3780605117952399,
  0.26031208092491054,
  0.5635124119976632,
  0.9790658438505838,
  0.8562823856318246,
  0.21556298702180277,
  0.8600632971753791,
  0.662714633786504,
  0.2225621933588111,
  0.6457530747930535,
  0.7827105700278855,
  0.6705869303441022,
  0.5154710337106151,
  0.815454332575039,
  0.6179902227520485,
  0.7115313466684177,
  0.9378033055153567,
  0.21433529585823752,
  0.8701474992411431,
  0.7086038807361402,
  0.30052303721084295,
  0.28393219786694557,
  0.5983530311667046,
  0.20020536916058207,
  0.6392286472751323,
  0.37143886775293566,
  0.6898805855917455,
  0.1884387811019529,
  0.5686068227042015,
  0.9620012698662417,
  0.4707056753390745,
  0.5257648252025556,
  0.6742146878570825,
  0.7082473720416732,
  0.13154771079490413,
  0.881639016990223,
  0.5319574884475743,
  0.37221870621745656,
  0.29767888275867493,
  0.7901841537170252,
  0.9446750496773592,
  0.9225501410767506,
  0.9805160330674125,
  0.37215404486327974,
  0.8896940430361793,
  0.4460397289641458,
  0.5596925008309813,
  0.972865691753777,
  0.09152757909635534,
  0.8255157060110575,
  0.5475708401593774,
  0.8558832930841987,
  0.944978975726182,
  0.265799720765052,
  0.7421827049142025,
  0.7365250938202301,
  0.9102159605707534);
}

// xspan has assertions on member access.
static constexpr auto rnd_table = xspan {detail::rnd_table_raw};

} // namespace artv
