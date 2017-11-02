// $Id: triangular.hpp 1081 2012-09-28 09:11:07Z perroe $

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

#ifndef NRLIB_RANDOM_TRIANGULAR_HPP
#define NRLIB_RANDOM_TRIANGULAR_HPP

#include "distribution.hpp"

#include <cassert>
#include <cmath>

namespace NRLib {

class Triangular : public Distribution<double> {
public:
  Triangular();
  Triangular(double min_val, double mode, double max_val);
  ~Triangular();

  Distribution<double>* Clone() const;

  inline void SetParameters(double min_val, double mode, double max_val);

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double x) const;

private:
  double min_val_;
  double mode_;
  double max_val_;
};


// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS -----------------

void Triangular::SetParameters(double min_val, double mode, double max_val)
{
  min_val_ = min_val;
  mode_ = mode;
  max_val_ = max_val;
}


double Triangular::Pdf(double x) const {
  if (x < min_val_ || x > max_val_) {
    return 0.0;
  }
  else if (x == min_val_) {
    return (mode_ == min_val_) ? 2.0 / (max_val_ - min_val_) : 0.0;
  }
  else if (x == max_val_) {
    return (mode_ == max_val_) ? 2.0 / (max_val_ - min_val_) : 0.0;
  }
  else if (x <= mode_) {
    return 2.0 * (x - min_val_) / ((mode_ - min_val_)*(max_val_ - min_val_));
  }
  else { // mode_ < x < max_val
    return 2.0 * (max_val_ - x) / ((max_val_ - mode_)*(max_val_ - min_val_));
  }
}


double Triangular::Cdf(double x) const {
  double p = 0.0;
  if (x > min_val_) {
    if (x <= mode_) {
      p = (x - min_val_) * (x - min_val_) / ((mode_ - min_val_)*(max_val_ - min_val_));
    }
    else if (x < max_val_) {
      p = 1.0 - (max_val_ - x) * (max_val_ - x) / ((max_val_ - mode_)*(max_val_ - min_val_));
    }
    else {
      p = 1.0;
    }
  }
  return p;
}


double Triangular::Quantile(double p) const {
  assert (p >= 0.0 && p <= 1.0);

  if (p == 0.0) {
    return min_val_;
  }
  if (p == 1.0) {
    return max_val_;
  }
  double p0 = (mode_ - min_val_) / (max_val_ - min_val_);
  if (p < p0) {
    return std::sqrt((max_val_ - min_val_) * (mode_ - min_val_) * p) + min_val_;
  }
  else if (p == p0) {
    return mode_;
  }
  else {
    return max_val_ - std::sqrt((max_val_ - min_val_) * (max_val_ - mode_) * (1.0 - p));
  }
}

} // namespace NRLib

#endif // NRLIB_RANDOM_TRIANGULAR_HPP
