#pragma once

#if 1

#include "artv-common/misc/mp11.hpp"

namespace artv {

template <class... Ts>
using typelist = mp_list<Ts...>;

}

#else
// This is very naîve Better to use mp11, they benchmark compiler performance,
// ifdef based on compiler features, etc...

#include <type_traits>

namespace artv {

//------------------------------------------------------------------------------
template <int...>
struct intlist {};
//------------------------------------------------------------------------------
template <class T>
struct intlist_sizeof;

template <int... Vs>
struct intlist_sizeof<intlist<Vs...>> {
  static constexpr unsigned value = sizeof...(Vs);
};

template <int... Vs1, int... Vs2>
constexpr auto operator+ (intlist<Vs1...>, intlist<Vs2...>)
{
  return intlist<Vs1..., Vs2...> {};
}
//------------------------------------------------------------------------------
template <class...>
struct typelist {};
//------------------------------------------------------------------------------
template <class T>
struct typelist_first {};

template <class T, class... Ts>
struct typelist_first<typelist<T, Ts...>> {
  using type = T;
};
//------------------------------------------------------------------------------
template <class... Ts>
constexpr typelist<Ts...> make_typelist (Ts... args)
{
  return typelist<Ts...> {};
};
//------------------------------------------------------------------------------
template <class... Ts>
constexpr size_t typelist_sizeof (typelist<Ts...>)
{
  return sizeof...(Ts);
};
//------------------------------------------------------------------------------
template <class... Ts1, class... Ts2>
constexpr auto operator+ (typelist<Ts1...>, typelist<Ts2...>)
{
  return typelist<Ts1..., Ts2...> {};
}
#if 0
template <class... Ts1, class... Ts2>
constexpr auto operator+ (typelist<Ts1...>, typelist<typelist<Ts2...>>)
{
  return typelist<Ts1..., Ts2...> {};
}
#endif
//------------------------------------------------------------------------------
namespace detail {

template <class Td, class Ts>
struct typelist_flatten;

template <class... Td>
struct typelist_flatten<typelist<Td...>, typelist<>> {
  using type = typelist<Td...>;
};

template <class... Td, class T, class... Ts>
struct typelist_flatten<typelist<Td...>, typelist<T, Ts...>>
  : public typelist_flatten<typelist<Td..., T>, typelist<Ts...>> {};

template <class... Td, class... Tn, class... Ts>
struct typelist_flatten<typelist<Td...>, typelist<typelist<Tn...>, Ts...>>
  : public typelist_flatten<typelist<Td...>, typelist<Tn..., Ts...>> {};

} // namespace detail

template <class T>
struct typelist_flatten : public detail::typelist_flatten<typelist<>, T> {};

#if 1 /*test/sanity check/documentation*/

using tflatten_data = typelist<
  int,
  typelist<typelist<float>>,
  typelist<typelist<typelist<char>>>>;

static_assert (
  std::is_same<
    typelist_flatten<tflatten_data>::type,
    typelist<int, float, char>>::value,
  "fail");
#endif

//------------------------------------------------------------------------------
struct typelist_has_args {};
struct typelist_no_more_args {};
struct typelist_null {};
//------------------------------------------------------------------------------
template <class T>
struct typelist_fwit;

template <>
struct typelist_fwit<typelist<>> {
  using head               = typelist_no_more_args;
  using tail               = typelist<>;
  using has_more_args_type = typelist_no_more_args;
};

template <class T, class... Ts>
struct typelist_fwit<typelist<T, Ts...>> {
  using head               = T;
  using tail               = typelist<Ts...>;
  using has_more_args_type = typelist_has_args;
};
//------------------------------------------------------------------------------
namespace detail {
template <class T, int N, class Tl>
struct typelist_type_index
  : public typelist_type_index<T, N + 1, typelist_fwit<typename Tl::tail>> {};

// found
template <int N, class Tl>
struct typelist_type_index<typename Tl::head, N, Tl> {
  static constexpr int value = N;
};

// not found
template <class T, int N>
struct typelist_type_index<T, N, typelist_fwit<typelist<>>> {
  static constexpr int value = -1;
};
} // namespace detail

// definition
template <class T, class Tl>
struct typelist_type_index;

// specialization
template <class T, class... Ts>
struct typelist_type_index<T, typelist<Ts...>>
  : public detail::typelist_type_index<T, 0, typelist_fwit<typelist<Ts...>>> {};

#if 0 /*test/sanity check/documentation*/
constexpr int v = typelist_type_index<int, typelist<float, float, int>>::value;
static_assert (v == 2, "fail");
constexpr int v2
  = typelist_type_index<char, typelist<float, float, int>>::value;
static_assert (v2 == -1, "fail");
#endif
//------------------------------------------------------------------------------
namespace detail {
template <int N, int I, class T_fwd, class T_back>
struct typelist_split;

template <int N, int I, class T, class... T_fwd, class... T_back>
struct typelist_split<N, I, typelist<T, T_fwd...>, typelist<T_back...>>
  : public typelist_split<
      N,
      I + 1,
      typelist<T_fwd...>,
      typelist<T_back..., T>> {};

template <int N, class T, class... T_fwd, class... T_back>
struct typelist_split<N, N, typelist<T, T_fwd...>, typelist<T_back...>> {
  using type      = T;
  using backwards = typelist<T_back...>;
  using forwards  = typelist<T_fwd...>;
};

template <int N, int I, class... T_back>
struct typelist_split<N, I, typelist<>, typelist<T_back...>> {
  using type      = typelist_no_more_args;
  using backwards = typelist<T_back...>;
  using forwards  = typelist<>;
};

} // namespace detail

template <int N, class T>
struct typelist_split;

template <int N, class... Ts>
struct typelist_split<N, typelist<Ts...>>
  : public detail::typelist_split<N, 0, typelist<Ts...>, typelist<>> {
  static_assert (N >= 0, "invalid range");
};

#if 0 /*test/sanity check/documentation*/

using split_0 = typelist_split<0, typelist<int, float, char>>;
static_assert (std::is_same<split_0::type, int>::value, "fail");
static_assert (
  std::is_same<split_0::forwards, typelist<float, char>>::value,
  "fail");
static_assert (std::is_same<split_0::backwards, typelist<>>::value, "fail");

using split_1 = typelist_split<1, typelist<int, float, char>>;
static_assert (std::is_same<split_1::type, float>::value, "fail");
static_assert (std::is_same<split_1::forwards, typelist<char>>::value, "fail");
static_assert (std::is_same<split_1::backwards, typelist<int>>::value, "fail");

using split_2 = typelist_split<2, typelist<int, float, char>>;
static_assert (std::is_same<split_2::type, char>::value, "fail");
static_assert (std::is_same<split_2::forwards, typelist<>>::value, "fail");
static_assert (
  std::is_same<split_2::backwards, typelist<int, float>>::value,
  "fail");

using split_3 = typelist_split<3, typelist<int, float, char>>;
static_assert (
  std::is_same<split_3::type, typelist_no_more_args>::value,
  "fail");
static_assert (std::is_same<split_3::forwards, typelist<>>::value, "fail");
static_assert (
  std::is_same<split_3::backwards, typelist<int, float, char>>::value,
  "fail");

#endif
//------------------------------------------------------------------------------
namespace detail {

template <
  int N,
  template <class, int>
  class Pred,
  class Unvisited_tl,
  class Match_tl,
  class Match_idxs_tl,
  class Reject_tl,
  class Reject_idxs_tl>
struct typelist_filter;

template <
  int N,
  template <class, int>
  class Pred,
  class T,
  class... Ts,
  class... Match,
  int... Match_idxs,
  class... Reject,
  int... Reject_idxs>
struct typelist_filter<
  N,
  Pred,
  typelist<T, Ts...>,
  typelist<Match...>,
  intlist<Match_idxs...>,
  typelist<Reject...>,
  intlist<Reject_idxs...>>
  : public typelist_filter<
      N + 1,
      Pred,
      typelist<Ts...>,
      typename std::conditional<
        Pred<T, N>::value,
        typelist<Match..., T>,
        typelist<Match...>>::type,
      typename std::conditional<
        Pred<T, N>::value,
        intlist<Match_idxs..., N>,
        intlist<Match_idxs...>>::type,
      typename std::conditional<
        !Pred<T, N>::value,
        typelist<Reject..., T>,
        typelist<Reject...>>::type,
      typename std::conditional<
        !Pred<T, N>::value,
        intlist<Reject_idxs..., N>,
        intlist<Reject_idxs...>>::type> {};

template <
  int N,
  template <class, int>
  class Pred,
  class... Match,
  int... Match_idxs,
  class... Reject,
  int... Reject_idxs>
struct typelist_filter<
  N,
  Pred,
  typelist<>,
  typelist<Match...>,
  intlist<Match_idxs...>,
  typelist<Reject...>,
  intlist<Reject_idxs...>> {
  using matched          = typelist<Match...>;
  using matched_indexes  = intlist<Match_idxs...>;
  using rejected         = typelist<Reject...>;
  using rejected_indexes = intlist<Reject_idxs...>;
  using type             = matched;
};

} // namespace detail

// filters a typelist based on a predicate/filter of this form
//
// template <class T, size_t index>
// struct predicate;
//
// The current type and index of the current iterated element would be passed to
// the predicate. It has to match the interface of a bool
// "std::integral_constant"

template <template <class, int> class Pred, class T>
struct typelist_filter;

template <template <class, int> class Pred, class... types>
struct typelist_filter<Pred, typelist<types...>>
  : public detail::typelist_filter<
      0,
      Pred,
      typelist<types...>,
      typelist<>,
      intlist<>,
      typelist<>,
      intlist<>> {};

#if 0 /*test/sanity check/documentation*/
template <class T, int>
using filter = std::is_same<T, float>;

using filtered = typelist_filter<filter, typelist<int, float, char>>;
static_assert (std::is_same<filtered::matched, typelist<float>>::value, "fail");
static_assert (
  std::is_same<filtered::matched_indexes, intlist<1>>::value,
  "fail");
static_assert (
  std::is_same<filtered::rejected, typelist<int, char>>::value,
  "fail");
static_assert (
  std::is_same<filtered::rejected_indexes, intlist<0, 2>>::value,
  "fail");
#endif
/*----------------------------------------------------------------------------*/
namespace detail {
template <class Count, class F>
void typelist_foreach (F& func, Count c)
{}

template <class Count, class F, class T, class... Ts>
void typelist_foreach (F& func, Count c, T, Ts...)
{
  func (c, T {});
  auto next = std::integral_constant<int, Count::value + 1> {};
  detail::typelist_foreach (func, next, Ts {}...);
}
} // namespace detail

template <class... Ts, class F>
void typelist_foreach (typelist<Ts...>, F func)
{
  detail::typelist_foreach (func, std::integral_constant<int, 0> {}, Ts {}...);
}
/*----------------------------------------------------------------------------*/
namespace detail {

template <class Count, class F>
void typelist_foreach (F& func, Count, typelist<>, typelist<>)
{}

template <class Count, class F, class T, class... Ts>
void typelist_foreach (F& func, Count, typelist<T, Ts...>, typelist<>)
{}

template <class Count, class F, class T, class... Ts>
void typelist_foreach (F& func, Count, typelist<>, typelist<T, Ts...>)
{}

template <class Count, class F, class T1, class T2, class... Ts1, class... Ts2>
void typelist_foreach (
  F&    func,
  Count c,
  typelist<T1, Ts1...>,
  typelist<T2, Ts2...>)
{
  func (c, T1 {}, T2 {});
  auto next = std::integral_constant<int, Count::value + 1> {};
  detail::typelist_foreach (
    func, next, typelist<Ts1...> {}, typelist<Ts2...> {});
}
} // namespace detail

// variant that accepts two typelists, stops when one of them is empty.
template <class... Ts1, class... Ts2, class F>
void typelist_foreach (typelist<Ts1...>, typelist<Ts2...>, F func)
{
  detail::typelist_foreach (
    func,
    std::integral_constant<int, 0> {},
    typelist<Ts1...> {},
    typelist<Ts2...> {});
}
/*----------------------------------------------------------------------------*/
namespace detail {
template <class F>
void typelist_invoke_on_index (F&, unsigned, unsigned)
{}

template <class F, class T, class... Ts>
void typelist_invoke_on_index (F& func, unsigned index, unsigned i, T, Ts...)
{
  if (i == index) {
    func (T {});
    return;
  }
  detail::typelist_invoke_on_index (func, index, ++i, Ts {}...);
}
} // namespace detail

template <class Functor, class... Ts>
void typelist_invoke_on_index (unsigned index, typelist<Ts...>, Functor func)
{
  detail::typelist_invoke_on_index (func, index, 0, Ts {}...);
}
/*----------------------------------------------------------------------------*/
namespace detail {
template <int Freq, int Start_idx>
struct get_select_every_n_filter {

  template <int Index>
  static constexpr bool match()
  {
    return (Start_idx <= Index) && ((Index - Start_idx) % Freq) == 0;
  }

  template <class T, int Index>
  struct type : public std::integral_constant<bool, match<Index>()> {};
};

} // namespace detail

template <class T, int Freq, int Start_idx = 0>
struct typelist_select_multiples
  : public typelist_filter<
      detail::get_select_every_n_filter<Freq, Start_idx>::template type,
      T> {};

#if 0 /*test/sanity check/documentation*/
using tf_test_data = typelist<int, float, double, int, float, double, int>;

template <class T>
struct print;

print<typelist_select_every_n<tf_test_data, 3, 2>::type> p;

static_assert (
  std::is_same<
    typelist_select_every_n<tf_test_data, 3, 2>::type,
    typelist<double, double>>::value,
  "fail");
#endif
/*----------------------------------------------------------------------------*/
namespace detail {
template <class Count, class F, class... Tfs>
auto typelist_transform (F& func, Count, typelist<Tfs...>)
{
  return typelist<Tfs...> {};
}

template <class Count, class F, class T, class... Tfs, class... Ts>
auto typelist_transform (F& func, Count c, typelist<Tfs...>, T, Ts...)
{
  auto transfmed = func (c, T {});
  auto next      = std::integral_constant<int, Count::value + 1> {};
  return detail::typelist_transform (
    func, next, typelist<Tfs..., decltype (transfmed)> {}, Ts {}...);
}
} // namespace detail

template <class... Ts, class F>
auto typelist_transform (typelist<Ts...>, F func)
{
  return detail::typelist_transform (
    func, std::integral_constant<int, 0> {}, typelist<> {}, Ts {}...);
}
/*----------------------------------------------------------------------------*/
} // namespace artv

#endif // #if 0
