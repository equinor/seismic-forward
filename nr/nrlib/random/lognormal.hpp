// $Id: lognormal.hpp 1228 2014-01-03 12:30:11Z gudmundh $

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

#ifndef NRLIB_RANDOM_LOGNORMAL_HPP
#define NRLIB_RANDOM_LOGNORMAL_HPP

#include "distribution.hpp"
#include "functions.hpp"
#include "normal.hpp"


#include "../math/constants.hpp"
#include "../exception/exception.hpp"
#include <cassert>
#include <cmath>

namespace NRLib {

class LogNormal : public Distribution<double> {
public:
  LogNormal();
  LogNormal(double mu, double sigma);
  ~LogNormal();

  Distribution<double>* Clone() const;

  inline void   SetParameters(double mu, double sigma);
  inline double Draw() const;
  inline double Draw(RandomGenerator & g) const;

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double x) const;

private:
  Normal normal_;
  double mu_; // mean in normal distr.
  double sigma_; //stdev in normal distr.
  double exp_; // Expectation in lognormal distribution
  double std_dev_; // Standard deviation in lognormal distribution
};

void LogNormal::SetParameters(double mu, double sigma)
{
  // If mean and std. dev. are in the log normal distribution, we calculate the corresponding
  // mu and sigma in the normal distribution
  normal_.SetParameters(mu_, sigma_);
  mu_ = mu;
  sigma_ = sigma;
  exp_ = std::exp(mu_ + 0.5 * sigma_ * sigma_);
  std_dev_ = std::sqrt((std::exp(sigma_*sigma_)-1.0) * std::exp(2.0 * mu_ + sigma_ * sigma_));
}

double LogNormal::Draw() const {
  double normal = normal_.Draw();
  double value =  std::exp(normal);
  return value;
}

double LogNormal::Draw(RandomGenerator & g) const {
  double normal = normal_.Draw(g);
  double value =  std::exp(normal);
  return value;
}

double LogNormal::Pdf(double x) const {
  double var = sigma_*sigma_;
  return std::exp(-(log(x)-mu_)*(log(x)-mu_)/(2.0*var)) / (sqrt(2.0*NRLib::Pi)*sigma_*x);
}

double LogNormal::Cdf(double x) const {
  double z = (log(x) - mu_)/sigma_;
  return erfc(-z * NRLib::Sqrt1_2) / 2.0;
}

double LogNormal::Quantile(double p) const {
  assert (p >= 0.0 && p <= 1.0);
  return exp(mu_ + sigma_ * NRLib::Sqrt2 * (-erfc_inv(2*p))); // is this correct??
}



}

#endif
