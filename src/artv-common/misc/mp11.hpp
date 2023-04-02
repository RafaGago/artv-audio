#pragma once

// some OSS library lifting the DSP as a GIT submodule, CMake subproject, etc.
// might not need to pull MP11.
#ifndef ARTV_NO_MP11

#include "artv-common/misc/short_ints.hpp"
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
template <uint End, class F>
void mp_foreach_idx (F&& f)
{
  mp11::mp_for_each<
    mp11::mp_from_sequence<mp11::make_integer_sequence<uint, End>>> (f);
}
//------------------------------------------------------------------------------
template <template <class...> class L, class... Ts, class F>
void mp_foreach_idx_reverse (L<Ts...>, F&& f)
{
  using with_idx = to_typelist_with_indexes<L<Ts...>>;

  mp11::mp_for_each<mp11::mp_reverse<with_idx>> ([&] (auto type) {
    f (typename decltype (type)::index {}, typename decltype (type)::type {});
  });
}
//------------------------------------------------------------------------------
template <uint Start, class F>
void mp_foreach_idx_reverse (F&& f)
{
  mp11::mp_for_each<mp11::mp_reverse<
    mp11::mp_from_sequence<mp11::make_integer_sequence<uint, Start + 1>>>> (f);
}
//------------------------------------------------------------------------------
// clang-format off
template <class L, class V>
using mp_not_in_list =
  mp11::mp_bool<mp11::mp_find<L, V>::value < mp11::mp_size<L>::value>;
// clang-format on

// "mp_not_in_list" as a quoted metafunction for a given list.
template <class L>
struct mp_not_in_list_qm {
  template <class V>
  using fn = mp_not_in_list<L, V>;
};
//------------------------------------------------------------------------------
// Removes the types on the list "L_rm" from list "L"
template <class L, class L_rm>
using mp_remove_all = mp11::mp_remove_if_q<L, mp_not_in_list_qm<L_rm>>;
//------------------------------------------------------------------------------
template <uint Offset, class T>
struct add_offset;

template <uint Offset, uint... Idx>
struct add_offset<Offset, std::index_sequence<Idx...>> {
  using type = std::index_sequence<(Offset + Idx)...>;
};

template <uint Offset, uint... Idx>
struct add_offset<Offset, mp11::index_sequence<Idx...>> {
  using type = mp11::index_sequence<(Offset + Idx)...>;
};

template <uint Offset, class T>
using add_offset_t = typename add_offset<Offset, T>::type;
//------------------------------------------------------------------------------

} // namespace artv

#else

namespace artv {

template <class... T>
struct mp_list {}
};
}

#endif
