// $Id: functions.hpp 883 2011-09-26 09:17:05Z perroe $

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

#ifndef NRLIB_RANDOM_FUNCTIONS_HPP
#define NRLIB_RANDOM_FUNCTIONS_HPP

/// @file Functions for use in statistical distributions.

#include <cstddef>

namespace NRLib {

// We use the "standard" names, as given in C99 and boost.

/// Error function.
double erf(double z);

/// Complementary error function, 1 - erf(z).
double erfc(double z);

/// Inverse error function.
double erf_inv(double z);

/// Complementary inverse error function, 1 - erf_inv(z).
double erfc_inv(double z);

// Evaluation of polynomials, Based on boost version.
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Compile time sized polynomials, just inline forwarders to the
// implementations above:
//
template <std::size_t N, class T, class V>
inline V evaluate_polynomial(const T(&a)[N], const V& val)
{
   V sum = static_cast<V>(a[N - 1]);
   for(int i = static_cast<int>(N) - 2; i >= 0; --i)
   {
      sum *= val;
      sum += static_cast<V>(a[i]);
   }
   return sum;
}


} // namespace NRLib

#endif // NRLIB_RANDOM_FUNCTIONS_HPP
