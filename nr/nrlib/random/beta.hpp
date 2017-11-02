// $Id: beta.hpp 1143 2013-04-03 09:01:08Z veralh $

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

#ifndef NRLIB_RANDOM_BETA_HPP
#define NRLIB_RANDOM_BETA_HPP

#include "distribution.hpp"
#include "../iotools/logkit.hpp"

namespace NRLib {

class Beta : public Distribution<double> {
public:
  Beta(double min_val, double max_val, double a, double b, int nx = 200);
  Beta();
  ~Beta();

  Distribution<double>* Clone() const;

  inline double Cdf(double x) const;
  inline double Pdf(double x) const;
  inline double Quantile(double p) const;

  double GetMinVal()    const { return min_val_; }
  double GetMaxVal()    const { return max_val_; }
  double GetInterval()  const { return max_val_ - min_val_; }

private:
   void    ComputePdfAndCdf();
   void    ComputeQuantile();

   double  min_val_;
   double  max_val_;
   double  a_;
   double  b_;
   int     nx_;        // number of bins in intervall [minVal_, maxVal_]
   double  dx_;        // width of one bin
   std::vector<double> pdf_;       // pdf at centerpoint of bin
   std::vector<double> cdf_;       // cdf at left side of bin
   std::vector<double> quantile_;  // cdf at left side of bin

};

// ---------------- IMPLEMENTATION OF INLINE FUNCTIONS ------------

double
Beta::Pdf(double x) const
{
  // round x to cell-center
  int ix = int( (x-min_val_)/(max_val_-min_val_) * nx_-0.5 );
  double d = 0.0;
  if(x >= min_val_ && x <= max_val_)
   d = pdf_[ix];

  return d;
}

double
Beta::Cdf(double x) const
{
  // linear interpolation between cell bounaries
  int ix = int( (x-min_val_)/dx_);

  double lowP  = 0.0;
  double p     = 0.0;
  double w     = 0.0;
  double highP = 1.0;


  if(x <= min_val_)
   p=0.0;
  else
    if(x >= max_val_)
      p=1.0;
    else
    {
      if(ix==0)
        lowP=0.0;
      else
        lowP = cdf_[ix-1];
      highP = cdf_[ix];

      w = ((x-min_val_)/dx_-ix);
      p = (1-w)*lowP+ w*highP;
    }
  return p;
}

double
Beta::Quantile(double p) const
{
  // Inverse of CDF

  // linear interpolation between cell bounaries
  int ip = int( p*nx_);

  double lowX  = min_val_;
  double x     = min_val_;
  double w     = 0.0;
  double highX = max_val_;


  if(p <= 0)
   x=min_val_;
  else
    if(p >= 1.0)
      x = max_val_;
    else
    {
      if(ip==0)
        lowX=min_val_;
      else
        lowX = quantile_[ip-1];
      highX = quantile_[ip];

      w = (p*nx_-ip);
      x = (1-w)*lowX+ w*highX;
    }
  return x;
}

} //namespace NRLib

#endif // NRLIB_RANDOM_BETA_HPP
