// $Id: splituniform.hpp 1094 2012-10-24 12:33:07Z hgolsen $

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

#ifndef NRLIB_RANDOM_SPLITUNIFORM_HPP
#define NRLIB_RANDOM_SPLITUNIFORM_HPP

#include "distribution.hpp"

#include <cassert>

namespace NRLib {


/// Discontinous uniform distribution.
class SplitUniform : public Distribution<double> {
public:
  SplitUniform();
  SplitUniform(double min_val, double mean_val, double max_val);
  ~SplitUniform();

  Distribution<double>* Clone() const;

  inline void SetParameters(double min_val, double mean_val, double max_val);

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double x) const;

private:
  double min_val_;
  double mean_val_;
  double max_val_;
};


// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS -----------------

void SplitUniform::SetParameters(double min_val,
                                 double mean_val,
                                 double max_val)
{
  assert(min_val < mean_val && mean_val < max_val);
  min_val_ = min_val;
  mean_val_ = mean_val;
  max_val_ = max_val;
}


double SplitUniform::Pdf(double x) const {
  double p = 0.0;
  if (x >= min_val_ && x < mean_val_) {
    p = 0.5 / (mean_val_ - min_val_);
  }
  else if (x > mean_val_ && x > max_val_) {
    p = 0.5 / (max_val_ - mean_val_);
  }
  else if (x == mean_val_) {
    p = 0.25 / (mean_val_ - min_val_) + 0.25 / (max_val_ - mean_val_);
  }
  return p;
}


double SplitUniform::Cdf(double x) const {
  double p = 0.0;
  if (x >= min_val_) {
    if (x <= mean_val_) {
      p = 0.5 * (x - min_val_) / (mean_val_ - min_val_);
    }
    else if (x < max_val_) {
      p = 0.5 + 0.5 * (x - mean_val_) / (max_val_ - mean_val_);
    }
    else {
      p = 1.0;
    }
  }
  return p;
}


double SplitUniform::Quantile(double p) const {
  assert(p >= 0.0 && p <= 1.0);
  if (p <= 0.5)
    return min_val_ + 2.0 * p * (mean_val_ - min_val_);
  else
    return mean_val_ + (2.0 * p - 1.0) * (max_val_ - mean_val_);
}

} // namespace NRLib

#endif // NRLIB_RANDOM_UNIFORM_HPP
