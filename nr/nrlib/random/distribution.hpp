// $Id: distribution.hpp 1169 2013-05-22 12:54:23Z anner $

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

#ifndef NRLIB_RANDOM_DISTRIBUTION_HPP
#define NRLIB_RANDOM_DISTRIBUTION_HPP

#include "random.hpp"
#include "randomgenerator.hpp"

namespace NRLib {

/// Abstract distribution class.
template<class A>
class Distribution {
public:
  virtual ~Distribution();

  /// Make a deep copy of the correct type.
  virtual Distribution<A>* Clone() const = 0;

  /// Probability density function.
  virtual double Pdf(A x) const = 0;

  /// Cumulative density function
  virtual double Cdf(A x) const = 0;

  /// Quantile - inverse of Cdf.
  /// @param z Quantile in the [0, 1] range.
  virtual A Quantile(double q) const = 0;

  /// Draw a number from the distribution.
  /// The random generator should already be initialized.
  virtual A Draw() const {
    // Default implementation.
    double q = NRLib::Random::Unif01();
    return Quantile(q);
  }

  virtual A Draw(RandomGenerator & g) const {
    // Default implementation.
    double q = g.Unif01();
    return Quantile(q);
  }

};

template<class A>
Distribution<A>::~Distribution() {}
}

#endif // NRLIB_RANDOM_DISTRIBUTION_HPP
