// $Id: binomial.cpp 1081 2012-09-28 09:11:07Z perroe $

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

#include "binomial.hpp"
#include "gamma.hpp"

#include <cassert>

namespace NRLib {

Binomial::Binomial(int n, double p)
  : n_(n), p_(p)
{
  assert(p_ >= 0.0 && p_ <= 1.0);

  // Compute cdf and store in table
  cdf_.resize(n_ + 1);
  cdf_[0] = Pdf(0);
  for(int i = 1; i < n_+1; i++)
    cdf_[i] = cdf_[i-1]+ Pdf(i);

}


Binomial::~Binomial()
{}


Distribution<int>*
Binomial::Clone() const
{
  return new Binomial(*this);
}


double Binomial::Pdf(int x) const
{
  if (x < 0 || x > n_)
    return 0.0;
  else if(p_ == 1.0 && x == n_)
    return 1.0;
  else if(p_ == 0.0 && x == 0)
    return 1.0;
  else if(p_ == 1.0 && x < n_)
    return 0.0;
  else if(p_ == 0.0 && x > 0)
    return 0.0;
  else
  {
    double pdf = bico(n_,x)*std::exp((x*std::log(p_))+((n_-x)*std::log(1.0-p_)));
    return pdf;
  }
}


double Binomial::Cdf(int x) const
{
  if (x < 0)
    return 0.0;
  else if(x > n_)
    return 1.0;
  else
    return cdf_[x];
}


int Binomial::Quantile(double q) const
{
  assert (q >= 0.0 && q <= 1.0);

  int high, low;
  low = 0;
  high = n_;
  int middle;
  while (high - low > 1)
  {
    middle = (high + low) >> 1;
    if (cdf_[middle] <=q)
      low = middle;
    else
      high = middle;
    //    middle = int(0.5*(high+low));
  }

  if (high == low)
    return low;
  else
    if (cdf_[high] > q && cdf_[low] <= q)
      return high;
    else
      return low;
}

} // namespace NRLib
