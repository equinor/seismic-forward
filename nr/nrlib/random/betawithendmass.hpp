// $Id: betawithendmass.hpp 1150 2013-04-09 09:00:00Z veralh $

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

// Class copied from PCube

#ifndef NRLIB_RANDOM_BETAWITHENDMASS_HPP
#define NRLIB_RANDOM_BETAWITHENDMASS_HPP

#include "distribution.hpp"
#include "beta.hpp"
#include "../iotools/logkit.hpp"

#include <limits>

namespace NRLib {

class BetaWithEndMass : public Distribution<double> {
public:
  BetaWithEndMass(double min_val, double max_val, double p_min, double p_max, double a, double b, int nx = 200);
  BetaWithEndMass();
  ~BetaWithEndMass();

  Distribution<double>* Clone() const;

  inline double Cdf(double x) const;
  inline double Pdf(double x) const;
  inline double Quantile(double p) const;

  double GetMinVal()     const { return beta_.GetMinVal();   }
  double GetMaxVal()     const { return beta_.GetMaxVal();   }
  double GetInterval()   const { return beta_.GetInterval(); }

private:

   NRLib::Beta    beta_;
   int     nx_;
   double  p_min_;
   double  p_max_;
};

// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS ------------

double
BetaWithEndMass::Pdf(double x) const
{
  double min_val = GetMinVal();
  double max_val = GetMaxVal();

  double d = 0.0;
  if(x == min_val)
    d = std::numeric_limits<double>::infinity(); // consistent with delta distribution
  if(x > min_val && x < max_val)
  {
    d = beta_.Pdf(x);
    // Rescaled by a factor (1-p_min_-p_max_)
    d = d*(1-p_min_-p_max_);
  }
  if(x == max_val)
    d = std::numeric_limits<double>::infinity();

  return d;
}

double
BetaWithEndMass::Cdf(double x) const
{
  double min_val = GetMinVal();
  double max_val = GetMaxVal();

  double p = 0.0;
  if(x < min_val)
    p = 0.0;
  else
    if(x >= max_val)
      p = 1.0;
    else
    {
      // Rescale by (1-p_min_-p_max_) and shift up p_min_.
      p = beta_.Cdf(x)*(1 - p_min_ - p_max_) + p_min_;
    }
  return p;
}

double
BetaWithEndMass::Quantile(double p) const
{
  // Inverse of CDF
  double min_val = GetMinVal();
  double max_val = GetMaxVal();

  double x = min_val;
  if(p <= p_min_)
    x = min_val;
  else
    if(p >= 1.0 - p_max_)
      x = max_val;
    else
    {
      // Rescale the probability value to use in the beta_.Quantile
      double p2 = (p - p_min_)/(1 - p_min_ - p_max_);
      x = beta_.Quantile(p2);
    }
  return x;
}

} //namespace NRLib

#endif // NRLIB_RANDOM_BETAWITHENDMASS_HPP
