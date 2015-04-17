// $Id: poisson.cpp 1235 2014-02-05 06:45:53Z anner $

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

#include "poisson.hpp"
#include "gamma.hpp"

#include <cassert>

namespace NRLib {

Poisson::Poisson(double lambda)
  : lambda_(lambda)
{
  exp_mean_ = std::exp(-lambda);
  smu_ = sqrt(lambda);
  b_ = 0.931 + 2.53 * smu_;
  a_ = -0.059 + 0.02483 * b_;
  inv_alpha_ = 1.1239 + 1.1328 / (b_ - 3.4);
  v_r_ = 0.9277 - 3.6224 / (b_ - 2);
  poisson_table_.resize(10);
  poisson_table_[0] = 0.0;
  poisson_table_[1] = 0.0;
  poisson_table_[2] = 0.69314718055994529;
  poisson_table_[3] = 1.7917594692280550;
  poisson_table_[4] = 3.1780538303479458;
  poisson_table_[5] = 4.7874917427820458;
  poisson_table_[6] = 6.5792512120101012;
  poisson_table_[7] = 8.5251613610654147;
  poisson_table_[8] = 10.604602902745251;
  poisson_table_[9] = 12.801827480081469;

}


Poisson::~Poisson()
{}


Distribution<int>*
Poisson::Clone() const
{
  return new Poisson(*this);
}


double Poisson::Pdf(int x) const
{
  int x_fac = 1;
  int i = 1;
  while (i < x) {
    i++;
    x_fac = x_fac * i;
  }
  double pdf = std::exp(x*std::log(lambda_)-lambda_)/x_fac;
  return pdf;
}


double Poisson::Cdf(int x) const
{
  int x_fac = 1;
  int i = 1;
  while (i < x) {
    i++;
    x_fac = x_fac * i;
  }
  double x1 = x + 1.0;
  double cdf = NRLib::gammp(x1, lambda_)/x_fac;
  return cdf;

}


int Poisson::Quantile(double q) const
{
   assert(q >= 0.0 && q <= 1.0);
   int i = 0;
   while (Cdf(i) < q)
     i++;
   return i;
}


// from boost
int Poisson::Draw() const
{

  if (lambda_ < 10) {
    double  p = exp_mean_;
    int x = 0;
    double u = NRLib::Random::Unif01();
    while(u > p) {
      u = u - p;
      ++x;
      p = lambda_ * p / x;
    }
    return x;

  }
  else {
    while(true) {
      double u;
      double v = NRLib::Random::Unif01();
      if(v <= 0.86 * v_r_) {
        u = v / v_r_ - 0.43;
        return static_cast<int>(floor(
          (2*a_/(0.5-abs(u)) + b_)*u + lambda_ + 0.445));
      }

      if(v >= v_r_) {
        u = NRLib::Random::Unif01() - 0.5;
      } else {
        u = v/v_r_ - 0.93;
        u = ((u < 0)? -0.5 : 0.5) - u;
        v = NRLib::Random::Unif01() * v_r_;
      }

      double us = 0.5 - abs(u);
      if(us < 0.013 && v > us) {
        continue;
      }

      double k = floor((2*a_/us + b_)*u+lambda_+0.445);
      v = v*inv_alpha_/(a_/(us*us) + b_);

      double log_sqrt_2pi = 0.91893853320467267;

      if(k >= 10) {
        if(log(v*smu_) <= (k + 0.5)*log(lambda_/k)
          - lambda_
          - log_sqrt_2pi
          + k
          - (1/12. - (1/360. - 1/(1260.*k*k))/(k*k))/k) {
            return static_cast<int>(k);
        }
      } else if(k >= 0) {
        if(log(v) <= k*log(lambda_)
          - lambda_
          - poisson_table_[static_cast<int>(k)]) {
            return static_cast<int>(k);
        }
      }
    }


  }

}


// from boost
int Poisson::Draw(RandomGenerator &g) const
{

  if (lambda_ < 10) {
     double  p = exp_mean_;
     int x = 0;
     double u = g.Unif01();
        while(u > p) {
            u = u - p;
            ++x;
            p = lambda_ * p / x;
        }
        return x;

  }
  else {
    while(true) {
            double u;
            double v = g.Unif01();
            if(v <= 0.86 * v_r_) {
                u = v / v_r_ - 0.43;
                return static_cast<int>(floor(
                    (2*a_/(0.5-abs(u)) + b_)*u + lambda_ + 0.445));
            }

            if(v >= v_r_) {
                u = g.Unif01() - 0.5;
            } else {
                u = v/v_r_ - 0.93;
                u = ((u < 0)? -0.5 : 0.5) - u;
                v = g.Unif01() * v_r_;
            }

            double us = 0.5 - abs(u);
            if(us < 0.013 && v > us) {
                continue;
            }

            double k = floor((2*a_/us + b_)*u+lambda_+0.445);
            v = v*inv_alpha_/(a_/(us*us) + b_);

            double log_sqrt_2pi = 0.91893853320467267;

            if(k >= 10) {
                if(log(v*smu_) <= (k + 0.5)*log(lambda_/k)
                               - lambda_
                               - log_sqrt_2pi
                               + k
                               - (1/12. - (1/360. - 1/(1260.*k*k))/(k*k))/k) {
                    return static_cast<int>(k);
                }
            } else if(k >= 0) {
                if(log(v) <= k*log(lambda_)
                           - lambda_
                           - poisson_table_[static_cast<int>(k)]) {
                    return static_cast<int>(k);
                }
            }
    }


  }

}


} // namespace NRLib
