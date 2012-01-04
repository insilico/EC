/**
 * \file StringUtils.h
 *
 * \brief Various string-related utilities.
 *
 * This is originally from Nate Barney circa Moore Lab days 2003-2007.
 * His function naming follows lowercase with underscores style,
 * while my additions are camelCase.
 *
 * \author Bill White, Nate Barney
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 10/7/04
 */

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <string>
#include <cctype>
#include <vector>
#include <clocale>
#include <functional>
#include <algorithm>
#include <iterator>
#include <climits>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace insilico
{
  // Predicate class to test if a character is of a given
  // class.  For example: is_classified<std::ctype_base::space>()
  // tests to see if a character is whitespace for the given
  // locale.  There are similar types for each of the is* functions
  // in <cctype>.
  template <std::ctype_base::mask Type, class charT = char>
  class is_classified : public std::unary_function<charT, bool>
  {
  public:
    // ctor from a ctype
    is_classified(std::ctype<charT> &ct) : m_ctype(ct) {
    }
    // ctor from a locale (for convenience)
    is_classified(const std::locale &loc = std::locale())
    : m_ctype(std::use_facet<std::ctype<charT> >(loc)) {
    }
    bool operator()(charT c) const {
      return m_ctype.is(Type, c);
    }
  private:
    std::ctype<charT> const & m_ctype;
  };

  // Unary functor to call toupper under the specified locale
  template <class charT = char>
  class do_to_upper : public std::unary_function<charT, charT>
  {
  public:
    // ctor from a ctype
    do_to_upper(std::ctype<charT> &ct) : m_ctype(ct) {
    }
    // ctor from a locale (for convenience)
    do_to_upper(const std::locale &loc = std::locale())
    : m_ctype(std::use_facet<std::ctype<charT> >(loc)) {
    }
    charT operator()(charT c) const {
      return m_ctype.toupper(c);
    }
  private:
    std::ctype<charT> const & m_ctype;
  };

  // Unary functor to call tolower under the specified locale

  template <class charT = char>
  class do_to_lower : public std::unary_function<charT, charT>
  {
  public:
    // ctor from a ctype
    do_to_lower(std::ctype<charT> &ct) : m_ctype(ct) {
    }
    // ctor from a locale (for convenience)
    do_to_lower(const std::locale &loc = std::locale())
    : m_ctype(std::use_facet<std::ctype<charT> >(loc)) {
    }
    charT operator()(charT c) const {
      return m_ctype.tolower(c);
    }
  private:
    std::ctype<charT> const & m_ctype;
  };

  // For all string functions, the stringT type must be a template instance
  // of std::basic_string.

  // remove leading whitespace
  template <typename stringT>
  stringT trim_left(const stringT &s, const std::locale &loc = std::locale()) {
    // find first non-whitespace character
    typename stringT::const_iterator it = std::find_if(s.begin(), s.end(),
                                                       std::not1(is_classified < std::ctype_base::space,
                                                                 typename stringT::value_type > (loc)));
    // return appropriate substring
    return stringT(it, s.end());
  }

  // remove trailing whitespace
  template <typename stringT>
  stringT trim_right(const stringT &s, const std::locale &loc = std::locale()) {
    // find last non-whitespace character
    typename stringT::const_reverse_iterator it = std::find_if(s.rbegin(),
                                                               s.rend(), std::not1(is_classified < std::ctype_base::space,
                                                                                   typename stringT::value_type > (loc)));
    // return appropriate substring
    return stringT(s.begin(), it.base());
  }

  // remove leading and trailing whitespace
  template <typename stringT>
  stringT trim(const stringT &s,
               const std::locale &loc = std::locale()) {
    typename stringT::const_iterator b;
    typename stringT::const_reverse_iterator e;
    // find first non-whitespace character
    b = std::find_if(s.begin(), s.end(),
                     std::not1(is_classified < std::ctype_base::space,
                               typename stringT::value_type > (loc)));
    // if none, return empty string
    if(b == s.end())
      return stringT();
    // find last non-whitespace character
    e = std::find_if(s.rbegin(), s.rend(),
                     std::not1(is_classified < std::ctype_base::space,
                               typename stringT::value_type > (loc)));
    // return appropriate substribg
    return stringT(b, e.base());
  }

  // For all split functions, the Container type must contain stringT's
  // and it must have a push_back(const stringT &) member to add items
  // to the end of the sequence.

  // perl-like split on whitespace
  template <typename Container, typename stringT>
  inline void split(Container &cont, const stringT &s,
                    const std::locale &loc = std::locale()) {
    // split on any whitespace character
    split_if(cont, s, is_classified < std::ctype_base::space,
             typename stringT::value_type > (loc));
  }
  // perl-like split with delimiter
  template <typename Container, typename stringT>
  void split(Container &cont, const stringT &s, const stringT &delim) {
    typename stringT::size_type i, j;
    // loop through positions where delimiter is found
    for(i = 0; i < s.size(); i = j + delim.size()) {
      // get next delimiter pos
      j = s.find(delim, i);
      // if not found, return remainder of string
      if(j == std::string::npos) {
        cont.push_back(s.substr(i));
        break;
      }
      // store current field, if nonempty
      if(j != i)
        cont.push_back(s.substr(i, j - i));
    }
  }

  // perl-like split with predicate
  // (pred must be a unary predicate that takes a stringT::value_type
  // as input and returns a boolean.  It is used to determine whether
  // a given character is a delimiter.)
  template <typename Container, typename stringT, typename Pred>
  void split_if(Container &cont, const stringT &s, const Pred &pred) {
    typename stringT::const_iterator i, j;
    // loop through positions where delimiter is found
    for(i = s.begin(); i != s.end(); i = j + 1) {
      // get next delimiter pos
      j = find_if(i, s.end(), pred);
      // if not found, return remainder of string
      if(j == s.end()) {
        cont.push_back(s.substr(i - s.begin()));
        break;
      }
      // store current field, if nonempty
      if(j != i)
        cont.push_back(s.substr(i - s.begin(), j - i));
    }
  }

  // perl-like join
  // It type must be an iterator that points to stringT's.
  template <typename It, typename stringT>
  stringT join(const It &begin, const It &end, const stringT &delim) {
    stringT ret;
    // loop through input fields
    for(It i = begin; i != end; ++i) {
      // add delimiter
      if(i != begin)
        ret += delim;
      // add field
      ret += *i;
    }
    // return result
    return ret;
  }

  // return uppercased copy of string
  template <typename stringT>
  stringT to_upper(const stringT &str,
                   const std::locale &loc = std::locale()) {
    stringT s = str;
    std::transform(s.begin(), s.end(), s.begin(),
                   do_to_upper<typename stringT::value_type > (loc));
    return s;
  }

  // return lowercased copy of string
  template <typename stringT>
  stringT to_lower(const stringT &str,
                   const std::locale &loc = std::locale()) {
    stringT s = str;
    std::transform(s.begin(), s.end(), s.begin(),
                   do_to_lower<typename stringT::value_type > (loc));
    return s;
  }

  // char and wchar_t versions
  inline std::string trim_left(const char *s,
                               const std::locale &loc = std::locale()) {
    return trim_left(std::string(s), loc);
  }
  inline std::wstring trim_left(const wchar_t *s,
                                const std::locale &loc = std::locale()) {
    return trim_left(std::wstring(s), loc);
  }
  inline std::string trim_right(const char *s,
                                const std::locale &loc = std::locale()) {
    return trim_right(std::string(s), loc);
  }
  inline std::wstring trim_right(const wchar_t *s,
                                 const std::locale &loc = std::locale()) {
    return trim_right(std::wstring(s), loc);
  }
  inline std::string trim(const char *s,
                          const std::locale &loc = std::locale()) {
    return trim(std::string(s), loc);
  }
  inline std::wstring trim(const wchar_t *s,
                           const std::locale &loc = std::locale()) {
    return trim(std::wstring(s), loc);
  }
  template <typename Container>
  inline void split(Container &cont, const char *s,
                    const std::locale &loc = std::locale()) {
    split(cont, std::string(s), loc);
  }
  template <typename Container>
  inline void split(Container &cont, const wchar_t *s,
                    const std::locale &loc = std::locale()) {
    split(cont, std::wstring(s), loc);
  }
  template <typename Container>
  inline void split(Container &cont, const std::string &s, const char *delim) {
    split(cont, s, std::string(delim));
  }
  template <typename Container>
  inline void split(Container &cont, const char *s, const std::string &delim) {
    split(cont, std::string(s), delim);
  }
  template <typename Container>
  inline void split(Container &cont, const char *s, const char *delim) {
    split(cont, std::string(s), std::string(delim));
  }
  template <typename Container>
  inline void split(Container &cont, const std::wstring &s, const wchar_t *delim) {
    split(cont, s, std::wstring(delim));
  }
  template <typename Container>
  inline void split(Container &cont, const wchar_t *s, const std::wstring &delim) {
    split(cont, std::wstring(s), delim);
  }
  template <typename Container>
  inline void split(Container &cont, const wchar_t *s, const wchar_t *delim) {
    split(cont, std::wstring(s), std::wstring(delim));
  }
  template <typename Container, typename Pred>
  inline void split_if(Container &cont, const char *s, const Pred &pred) {
    split_if(cont, std::string(s), pred);
  }
  template <typename Container, typename Pred>
  inline void split_if(Container &cont, const wchar_t *s, const Pred &pred) {
    split_if(cont, std::wstring(s), pred);
  }
  template <typename It>
  inline std::string join(const It &begin, const It &end, const char *delim) {
    return join(begin, end, std::string(delim));
  }
  template <typename It>
  inline std::wstring join(const It &begin, const It &end, const wchar_t *delim) {
    return join(begin, end, std::wstring(delim));
  }
  inline std::string to_upper(const char *s,
                              const std::locale &loc = std::locale()) {
    return to_upper(std::string(s), loc);
  }
  inline std::wstring to_upper(const wchar_t *s,
                               const std::locale &loc = std::locale()) {
    return to_upper(std::wstring(s), loc);
  }
  inline std::string to_lower(const char *s,
                              const std::locale &loc = std::locale()) {
    return to_lower(std::string(s), loc);
  }
  inline std::wstring to_lower(const wchar_t *s,
                               const std::locale &loc = std::locale()) {
    return to_lower(std::wstring(s), loc);
  }
  template<typename T>
  std::string get_bits(T value) {
    int size = sizeof(value) * CHAR_BIT;
    std::string ret;
    ret.reserve(size);
    for(int i = size - 1; i >= 0; --i)
      ret += (value & (1 << i)) == 0 ? '0' : '1';
    return ret;
  }
  template<typename T>
  std::string zeroPadNumber(T num, int padSize) {
    std::ostringstream ss;
    ss << std::setw(padSize) << std::setfill('0') << num;
    return ss.str();
  }
} // namespace insilico

#endif // STRINGUTIL_H
