#pragma once

// some OSS library lifting the DSP as a GIT submodule, CMake subproject, etc.
// might not need to pull MP11.
#ifndef ARTV_NO_MP11

#include <boost/mp11.hpp>

namespace artv {

namespace mp11 = boost::mp11;

template <class... T>
using mp_list = mp11::mp_list<T...>;

//------------------------------------------------------------------------------
template <class T, class Idx>
struct type_with_index {
  using type  = T;
  using index = Idx;
};

template <class T>
using type_with_index_type = typename T::type;

template <class T>
using type_with_index_index = typename T::index;
//------------------------------------------------------------------------------
namespace detail {

template <class T>
class to_typelist_with_indexes;

template <template <class...> class T, class... Ts>
class to_typelist_with_indexes<T<Ts...>> {
private:
  using index_seq       = mp11::mp_iota_c<sizeof...(Ts)>;
  using index_seq_typed = mp11::mp_rename<index_seq, T>;

public:
  using type = mp11::mp_transform<type_with_index, T<Ts...>, index_seq_typed>;
};
} // namespace detail

template <class Tl>
using to_typelist_with_indexes =
  typename detail::to_typelist_with_indexes<Tl>::type;
//------------------------------------------------------------------------------
template <template <class...> class L, class... Ts, class F>
void mp_foreach_idx (L<Ts...>, F&& f)
{
  using with_idx = to_typelist_with_indexes<L<Ts...>>;

  mp11::mp_for_each<with_idx> ([&] (auto type) {
    f (typename decltype (type)::index {}, typename decltype (type)::type {});
  });
}
//------------------------------------------------------------------------------

#else

namespace artv {

template <class... T>
struct mp_list {}
};
}

#endif
} // namespace artv
