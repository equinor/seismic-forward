// $Id: beta.cpp 1131 2012-12-20 17:25:05Z hauge $

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

#include "beta.hpp"
#include <cmath>
#include <cstdio>
#include <cassert>

namespace NRLib {

Beta::Beta(double min_val, double max_val, double a, double b, int nx)
{
  assert(a > 0.0);
  assert(b > 0.0);
  assert(max_val > min_val);
  min_val_  = min_val;
  max_val_  = max_val;
  a_        = a;
  b_        = b;
  nx_       = nx;
  dx_       = (max_val-min_val)/nx;

  pdf_.resize(nx);
  cdf_.resize(nx);
  quantile_.resize(nx);

  ComputePdfAndCdf();
  ComputeQuantile();
}

Beta::Beta()
{
}

Beta::~Beta()
{
}

Distribution<double>*
Beta::Clone() const
{
  return new Beta(*this);
}

void
Beta::ComputePdfAndCdf()
{
  double x,invCdfNormF,invPdfNormF, max_log_pdf;
  std::vector<double> log_pdf(nx_);
  double normF = 0.0;
  int i;

  //To avoid numerical problems, we first compute loglikelihoods.
  max_log_pdf = -1e+30;
  for(i=0;i<nx_;i++)
  {
    x          = (i+0.5)/double(nx_);
    log_pdf[i] = (a_-1.0)*log(x)+(b_-1.0)*log(1-x);
    if(log_pdf[i] > max_log_pdf)
      max_log_pdf = log_pdf[i];
  }

  //We are going to renormalise later, so an arbitrary scale will do here. Set maximum element to 1.
  for(i=0;i<nx_;i++) {
    pdf_[i] = exp(log_pdf[i] - max_log_pdf);
    normF += pdf_[i];
  }

  invCdfNormF = 1.0/normF;
  invPdfNormF = 1.0/(normF*dx_);

  cdf_[0]  = pdf_[0] * invCdfNormF;//NB note order important
  pdf_[0] *= invPdfNormF;

  for(i=1;i<nx_;i++)
  {
   cdf_[i]   = cdf_[i-1]+pdf_[i] * invCdfNormF; // NB note order important
   pdf_[i]  *= invPdfNormF;
  }
}


void
Beta::ComputeQuantile()
{
  int i;
  int cdfCounter = 0;
  double w=0.0;
  double dp=1.0/nx_;
  double lowX = min_val_;
  double highX= min_val_;
  double lowP = 0.0;;
  double highP= 0.0;;

  for(i=0;i<nx_;i++)
  {
    while( (cdf_[cdfCounter] < dp*(i+1) ) &&  ( cdfCounter < nx_-1)  )
      cdfCounter++;
    if(cdfCounter > 0)
    {
      lowP = cdf_[cdfCounter-1];
      lowX = min_val_+(cdfCounter)*dx_;
    }
    if(cdfCounter >= nx_)
    {
      cdfCounter = nx_-1;
    }
    highP = cdf_[cdfCounter];
    highX = min_val_+(cdfCounter+1)*dx_;
    if(highP <= lowP)
      w = 1.0;
    else
      w= (dp*(i+1)-lowP)/(highP-lowP);

    quantile_[i] = (1-w)*lowX + w*highX;
  }
}

} // namespace NRLib
