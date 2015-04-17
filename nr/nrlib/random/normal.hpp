// $Id: normal.hpp 1228 2014-01-03 12:30:11Z gudmundh $

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

#ifndef NRLIB_RANDOM_NORMAL_HPP
#define NRLIB_RANDOM_NORMAL_HPP

#include "distribution.hpp"
#include "functions.hpp"
#include "../math/constants.hpp"
#include "../exception/exception.hpp"
#include <cassert>
#include <cmath>

namespace NRLib {

class Normal : public Distribution<double> {
public:
  Normal();
  Normal(double mean, double std_dev);
  ~Normal();

  Distribution<double>* Clone() const;

  inline void   SetParameters(double mean, double std_dev);

  inline double Pdf(double x) const;
  inline double Cdf(double x) const;
  inline double Quantile(double p) const;

  inline double Draw() const;
  inline double Draw(RandomGenerator & g) const;
  double        DrawLowerTruncated(double lower) const; //Not tested
  double        DrawLowerTruncated(double lower, RandomGenerator &g) const; //Not tested
  double        DrawUpperTruncated(double upper) const; //Not tested
  double        DrawUpperTruncated(double upper, RandomGenerator &g) const; //Not tested
  double        DrawTwoSidedTruncated(double min, double max) const; // Not tested
  double        DrawTwoSidedTruncated(double min, double max, RandomGenerator &g) const; // Not tested

  double        PhiInverse(double y) const;  //Not tested
  inline double Potential(double x) const; // -Log of pdf
  inline double PotentialTwoSidedTruncated(double x, double min, double max) const;
  inline double PotentialZeroTruncated(double x) const;

  /// Quantile for Normal(0, 1) distribution.
  inline static double Quantile01(double p);

private:
  double mean_;
  double std_dev_;

  void   PHIMinusyAndphi(double x,
                         double y,
                         double & Phi_minus_y,
                         double & phi) const;

  double RootPHI(double x1,
                 double x2,
                 double xacc,
                 int MaxIt,
                 double y) const; // not tested

  double DrawLowerTruncatedAcceptReject(double mu,
                                        double sigma,
                                        double lower,
                                        RandomGenerator *g = NULL) const;

  double DrawTwoSidedTruncatedAcceptReject(double mu,
                                           double sigma,
                                           double lower,
                                           double upper,
                                           RandomGenerator *g = NULL) const;
};

// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS -----------------

void Normal::SetParameters(double mean, double std_dev) {
  mean_ = mean;
  std_dev_ = std_dev;
}


double Normal::Pdf(double x) const {
  double var = std_dev_*std_dev_;
  return std::exp(-(x-mean_)*(x-mean_)/(2.0*var)) / sqrt(2.0*NRLib::Pi*var);
}


double Normal::Cdf(double x) const {
  double z = (x - mean_)/std_dev_;
  return erfc(-z * NRLib::Sqrt1_2) / 2.0;
}


double Normal::Quantile(double p) const {
  assert (p >= 0.0 && p <= 1.0);
  return mean_ + std_dev_ * NRLib::Sqrt2 * (-erfc_inv(2*p));
}


double Normal::Draw() const {
  double z = NRLib::Random::Norm01();
  return mean_ + std_dev_ * z;
}


double Normal::Draw(RandomGenerator & g) const {
    // Default implementation.
    double z = g.Norm01();
    return mean_ + std_dev_ * z;
}


double Normal::Potential(double x) const {
  double var = std_dev_*std_dev_;
  double y = x - mean_;
  double potential = 0.5*(log(2.0*NRLib::Pi*var) + ((y*y)/var));
  return potential;
}


double Normal::PotentialTwoSidedTruncated(double x, double min, double max) const {
  assert(max > min);
  double scale = Cdf(max)-Cdf(min);
  double potential;
  if(scale !=0.0)
   potential = Potential(x)+log(scale);
  else
   potential = Potential(x)+log(max-min);
  return(potential);
}


double Normal::PotentialZeroTruncated(double x) const
{
  assert (x > 0.0);
  double scale = 1.0 - Cdf(0.0);
  double potential = Potential(x);
  if (scale > 0.0)
    potential += log(scale);
  return(potential);
}

double Normal::Quantile01(double p)
{
  assert (p >= 0.0 && p <= 1.0);
  return NRLib::Sqrt2 * (-erfc_inv(2*p));
}


} // namespace NRLib

#endif // NRLIB_RANDOM_NORMAL_HPP
