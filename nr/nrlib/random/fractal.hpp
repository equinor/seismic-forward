// $Id: fractal.hpp 1228 2014-01-03 12:30:11Z gudmundh $

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

#ifndef NRLIB_RANDOM_FRACTAL_HPP
#define NRLIB_RANDOM_FRACTAL_HPP

#include "distribution.hpp"
#include <cassert>
#include <cmath>

namespace NRLib {

class Fractal : public Distribution<double> {
public:
  Fractal(double min, double max, double dim);
  ~Fractal();

  Distribution<double>* Clone() const;

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double p) const;

  double        DrawFromUpperQuantile(double qu, bool has_limit_info = false, double limit = 0.0) const;
  double        DrawFromLowerQuantile(double ql) const;
  double        DrawFromMiddleQuantile(double ql, double qu) const;
  double        GetMin() const { return min_; }

private:
  double min_;
  double max_;
  double dim_;

  double min_fractal_;
  double max_fractal_;

};

// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS -----------------
double Fractal::Pdf(double x) const {
  double pdf(0.0);
  if (x > min_ && x < max_) {
    double constant = dim_ / (min_fractal_ - max_fractal_);
    pdf = constant * std::pow(x, -dim_ - 1.0);
  }
  return (pdf);
}


double Fractal::Cdf(double x) const {
  double cdf;
  if (x <= min_)
    cdf = 0.0;
  else if (x < max_) {
    double d = std::pow(x, -dim_);
    cdf = (min_fractal_ - d) / (min_fractal_ - max_fractal_);
  }
  else
    cdf = 1.0;
  return(cdf);
}


double Fractal::Quantile(double p) const {
  assert (p >= 0.0 && p <= 1.0);
  double d = (1.0 - p) * min_fractal_ + p * max_fractal_;
  double q = std::pow(d, - 1.0/dim_);
  return(q);
}


} // namespace NRLib

#endif // NRLIB_RANDOM_NORMAL_HPP
