// $Id: chisquared.hpp 1081 2012-09-28 09:11:07Z perroe $

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

#ifndef NRLIB_RANDOM_CHI_SQUARED_HPP
#define NRLIB_RANDOM_CHI_SQUARED_HPP

#include "distribution.hpp"
#include "gamma.hpp"
#include "../exception/exception.hpp"

#include <cassert>
#include <cmath>

namespace NRLib {

class ChiSquared : public Distribution<double>
{
public:
  ChiSquared(double degrees_of_freedom);

  ~ChiSquared();

  Distribution<double>* Clone() const;

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double q) const;

private:
  double df_;
  double fac_;
};


// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS ------------

double ChiSquared::Pdf(double x) const
{
  if (x < 0)
    return 0.0;
  else
    return std::exp(-0.5*(x-(df_-2)*std::log(x)) - fac_);
}


double ChiSquared::Cdf(double x) const
{
  if (x < 0)
    return 0.0;
  else
    return gammp(0.5*df_, 0.5*x);
}


double ChiSquared::Quantile(double q) const
{
  assert(q >= 0.0 && q <= 1.0);

  return 2*invgammp(q, 0.5*df_);
}


} // namespace NRLib

#endif // NRLIB_RANDOM_CHI_SQUARED_HPP
