// $Id: statistics.hpp 933 2011-12-12 13:28:02Z veralh $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NRLIB_STATISTICS_HPP
#define NRLIB_STATISTICS_HPP

#include <algorithm>
#include <iterator>

#include "../exception/exception.hpp"

namespace NRLib {
  class NoDataError : public Exception
  {
  public:
    explicit NoDataError(const std::string& msg = "")
      : Exception(msg) {}

    virtual ~NoDataError() throw() {}
  };

  class InvalidDataError : public Exception
  {
    public:
    explicit InvalidDataError(const std::string& msg = "")
      : Exception(msg) {}

    virtual ~InvalidDataError() throw() {}
  };

  template<class FI>
  int Count(FI begin, FI end);

  template<class FI, class M>
  int Count(FI begin, FI end, M IsMissing);

  template<class FI>
  FI MaxElement(FI begin, FI end);

  /// Returns end if all elements are missing.
  template<class FI, class M>
  FI MaxElement(FI begin, FI end, M isMissing);

  template<class FI>
  FI MinElement(FI begin, FI end);

  /// Returns end if all elements are missing.
  template<class FI, class M>
  FI MinElement(FI begin, FI end, M isMissing);

  template<class FI>
  double Sum(FI begin, FI end);

  template<class FI, class M>
  double Sum(FI begin, FI end, M isMissing);

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class FI>
  double Mean(FI begin, FI end, bool mayThrow = true);

  /// If there is no data or all data is missing an exception is thrown if
  /// mayThrow is true, else 0 is returned.
  template<class FI, class M>
  double Mean(FI begin, FI end, M isMissing, bool mayThrow = true);

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class A>
  double Mean(A collection, bool mayThrow = true)
  { return Mean (collection.begin(), collection.end(), mayThrow); }

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class A, class M>
  double Mean(A collection, M isMissing, bool mayThrow = true)
  { return Mean (collection.begin(), collection.end(),
                 isMissing, mayThrow);
  }

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class FI>
  double Var(FI begin, FI end, bool mayThrow = true);

  /// If there is no data or all data is missing an exception is thrown if
  /// mayThrow is true, else 0 is returned.
  template<class FI, class M>
  double Var(FI begin, FI end, M isMissing, bool mayThrow = true);

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class A>
  double Var(A collection, bool mayThrow = true)
  { return Var(collection.begin(), collection.end(), mayThrow); }

  /// If there is no data an exception is thrown if mayThrow is true, else
  /// 0 is returned.
  template<class A, class M>
  double Var(A collection, M isMissing, bool mayThrow = true)
  { return Var(collection.begin(), collection.end(), isMissing, mayThrow); }


  /// If there is no data or all data is missing an exception is thrown if
  /// mayThrow is true, else 0 is returned.
  template<class FIA, class FIB>
  double Cov(FIA beginA, FIA endA, FIB beginB, FIB endB,
             bool mayThrow = true);

  /// If there is no data or all data is missing an exception is thrown if
  /// mayThrow is true, else 0 is returned.
  template<class A, class B>
  double Cov(A collection1, B collection2, bool mayThrow = true)
  { return Cov(collection1.begin(), collection1.end(),
               collection2.begin(), collection2.end(), mayThrow); }
} // namespace NRLib


// ------ IMPLEMENTATION ------

template <class FI>
inline int NRLib::Count(FI begin, FI end)
{
  return static_cast<int>(std::distance(begin, end));
}


template <class FI, class M>
inline int NRLib::Count(FI begin, FI end, M IsMissing)
{
  int n = 0;
  for (FI it = begin; it != end; ++it) {
    if (!IsMissing(*it)) {
      ++n;
    }
  }
  return n;
}


template <class FI>
inline FI NRLib::MaxElement(FI begin, FI end)
{
  return std::max_element(begin, end);
}


template <class FI, class M>
inline FI NRLib::MaxElement(FI begin, FI end, M IsMissing)
{
  FI it = begin;
  while (it != end && IsMissing(*it)) {
    ++it;
  }
  if (it == end) {
    return it;
  }
  FI max = it;
  for (++it; it != end; ++it) {
    if (!IsMissing(*it) && *it > *max) {
      max = it;
    }
  }
  return max;
}


template <class FI>
inline FI NRLib::MinElement(FI begin, FI end)
{
  return std::min_element(begin, end);
}


template <class FI, class M>
inline FI NRLib::MinElement(FI begin, FI end, M IsMissing)
{
  FI it = begin;
  while (it != end && IsMissing(*it)) {
    ++it;
  }
  if (it == end) {
    return it;
  }
  FI min = it;
  for (++it; it != end; ++it) {
    if (!IsMissing(*it) && *it < *min) {
      min = it;
    }
  }
  return min;
}


template <class FI>
inline double NRLib::Sum(FI begin, FI end)
{
  double sum = 0.0;
  for (FI it = begin; it != end; ++it) {
    sum += *it;
  }
  return sum;
}


template <class FI, class M>
inline double NRLib::Sum(FI begin, FI end, M IsMissing)
{
  double sum = 0.0;
  for (FI it = begin; it != end; ++it) {
    if (!IsMissing(*it)) {
      sum += *it;
    }
  }
  return sum;
}


template <class FI>
inline double NRLib::Mean(FI begin, FI end, bool mayThrow)
{
  int c = NRLib::Count(begin, end);
  if (c == 0) {
    if (mayThrow) {
      throw NoDataError("No input data.");
    }
    else {
      return 0.0;
    }
  }
  return NRLib::Sum(begin, end) / c;
}

template <class FI, class M>
double NRLib::Mean(FI begin, FI end, M IsMissing, bool mayThrow)
{
  int c = NRLib::Count(begin, end, IsMissing);
  if (c == 0) {
    if (mayThrow) {
      throw NoDataError("No input data, or all data is missing.");
    }
    else {
      return 0.0;
    }
  }
  return NRLib::Sum(begin, end, IsMissing) / c;
}


template <class FI>
double NRLib::Var(FI begin, FI end, bool mayThrow)
{
  int c = NRLib::Count(begin, end);
  if (c == 0) {
    if (mayThrow) {
      throw NoDataError("No input data");
    }
    else {
      return 0.0;
    }
  }

  if (c == 1) {
    return 0.0;
  }

  double sum = 0.0;
  double mean = NRLib::Mean(begin, end, mayThrow);
  for (FI it = begin; it != end; ++it) {
    sum += (*it - mean) * (*it - mean);
  }
  return sum/(c-1);
}


template <class FI, class M>
double NRLib::Var(FI begin, FI end, M IsMissing, bool mayThrow)
{
  int c = NRLib::Count(begin, end, IsMissing);
  if (c == 0) {
    if (mayThrow) {
      throw NoDataError("No input data, or all data is missing.");
    }
    else {
      return 0.0;
    }
  }

  if (c == 1) {
    return 0.0;
  }

  double sum = 0.0;
  double mean = NRLib::Mean(begin, end, IsMissing, mayThrow);
  for (FI it = begin; it != end; ++it) {
    if (!IsMissing(*it)) {
      sum += (*it - mean) * (*it - mean);
    }
  }
  return sum/(c-1);
}

template <class FIA, class FIB>
double NRLib::Cov(FIA beginA, FIA endA, FIB beginB, FIB endB,
                   bool mayThrow)
{
  int cA = NRLib::Count(beginA, endA);
  int cB = NRLib::Count(beginB, endB);

  if (cA != cB) {
    if (mayThrow) {
      throw InvalidDataError("The two datasets have different lengths.");
    }
    else {
      return 0.0;
    }
  }
  if (cA == 0) {
    if (mayThrow) {
      throw NoDataError("No input data.");
    }
    else {
      return 0.0;
    }
  }

  double sum = 0.0;
  double meanA = NRLib::Mean(beginA, endA, mayThrow);
  double meanB = NRLib::Mean(beginB, endB, mayThrow);

  FIA itA = beginA;
  FIB itB = beginB;
  for ( ; itA != endA; ++itA, ++itB) {
    sum += (*itA - meanA) * (*itB - meanB);
  }
  return sum/(cA-1);
}


#endif // NRLIB_STATISTICS_HPP
