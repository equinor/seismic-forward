// $Id: normal.cpp 1220 2013-11-27 09:54:01Z hauge $

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


#include "normal.hpp"

#include "../iotools/stringtools.hpp"
#include "../math/mathutility.hpp"

namespace NRLib {

Normal::Normal()
  : mean_(0.0), std_dev_(1.0)
{}

Normal::Normal(double mean, double std_dev)
  : mean_(mean), std_dev_(std_dev)
{}

Normal::~Normal()
{}

Distribution<double>*
Normal::Clone() const
{
  return new Normal(*this);
}

double
Normal::PhiInverse(double y) const
{
  double lower_limit, upper_limit;

  if (y < 0.0 || y > 1.0)
    throw Exception("y in PhiInverse must be in [0,1]\n");

  if (y == 0.0) return - exp(1000.0);
  if (y == 1.0) return exp(1000.0);

  if (y > 0.5) {
    double step = std_dev_;
    lower_limit = mean_;
    upper_limit = mean_ + step;
    while (Cdf(upper_limit) < y) {
      lower_limit  = upper_limit;
      step        *= 10.0;
      upper_limit  = mean_ + step;
    }
  }
  else if (y < 0.5) {
    double step = std_dev_;
    upper_limit = mean_;
    lower_limit = mean_ - step;
    while (Cdf(lower_limit) > y) {
      upper_limit = lower_limit;
      step        *= 10.0;
      lower_limit  = mean_ - step;
    }
  }
  else
    return mean_;

  if (IsEqual(Cdf(upper_limit),y))
    return upper_limit;

  if (IsEqual(Cdf(lower_limit),y))
    return lower_limit;

  return RootPHI(lower_limit, upper_limit, 1.0e-10, 1000, y);
}

/* PHIMinusyAndphi ************************************************************

DESCRIPTION: computes PHI(x) - y and phi(x) with y being a global variable
             and PHI() and phi() being cummulative
             and density of standard Gaussian distribution


SIDE EFFECTS: none

RETURN VALUE: none

******************************************************************************/
void
Normal::PHIMinusyAndphi(double x, double y, double & Phi_minus_y, double & phi) const
{
  Phi_minus_y = Cdf(x) - y;
  phi         = Pdf(x);

  return;
}

/* rtsafe *********************************************************************

DESCRIPTION: find roots of specified function by Newton-Raphson algorithm.
             Function is taken from numerical recipes.

SIDE EFFECTS: none

RETURN VALUE: root

******************************************************************************/
double Normal::RootPHI(double lower,double upper,double xacc,int MaxIt, double y) const
{
  int j;
  double df,fh,fl;
  void nrerror();

  PHIMinusyAndphi(lower, y, fl, df);
  PHIMinusyAndphi(upper, y, fh, df);

  if (IsZero(fl))
    return lower;
  else if (IsZero(fh))
    return upper;

  if (fl*fh >= 0.0)
    throw Exception("Root must be bracketed in NRLib::Normal::RootPHI\n");

  double xh, xl;
  if (fl < 0.0) {
    xl = lower;
    xh = upper;
  }
  else {
    xh = lower;
    xl = upper;
  }
  double rts=0.5*(lower + upper);
  double dxold=fabs(upper - lower);
  double dx=dxold;

  double f;
  PHIMinusyAndphi(rts, y, f, df);

  for (j=1;j<=MaxIt;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
        || (fabs(2.0*f) > fabs(dxold*df))) {
      dxold = dx;
      dx = 0.5*(xh-xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else {
      dxold = dx;
      dx = f/df;
      double temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (IsZero(dx, xacc))
      return rts;

    PHIMinusyAndphi(rts, y, f, df);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  throw Exception("Maximum number of iterations exceeded in NRLib::Normal::RootPHI");
}


double Normal::DrawLowerTruncated(double lower) const
{
  if(IsZero(std_dev_)) {
    if(mean_ >= lower)
      return mean_;
    else
      throw Exception("Cannot draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }
  double out = DrawLowerTruncatedAcceptReject(mean_, std_dev_, lower);

  return(out);
}


double Normal::DrawLowerTruncated(double lower, RandomGenerator &g) const
{
  if(IsZero(std_dev_)) {
    if(mean_ >= lower)
      return mean_;
    else
      throw Exception("Cannot draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }

  double out = DrawLowerTruncatedAcceptReject(mean_, std_dev_, lower, &g);

  return(out);
}


double Normal::DrawUpperTruncated(double upper) const
{
  if(IsZero(std_dev_)) {
    if(mean_ <= upper)
      return mean_;
    else
      throw Exception("Cannot draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }
  upper = 2*mean_-upper;
  double out = DrawLowerTruncatedAcceptReject(mean_, std_dev_, upper);

  out = 2*mean_ - out;

  return(out);
}


double Normal::DrawUpperTruncated(double upper, RandomGenerator &g) const
{
  if(IsZero(std_dev_)) {
    if(mean_ <= upper)
      return mean_;
    else
      throw Exception("Cannot draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }
  upper = 2*mean_-upper;
  double out = DrawLowerTruncatedAcceptReject(mean_, std_dev_, upper, &g);

  out = 2*mean_ - out;

  return(out);
}


double Normal::DrawTwoSidedTruncated(double lower, double upper) const
{
  if(lower > upper) {
    std::string err_txt = "Wrong input to truncated normal, should have min <= max). (min, max) = ("
      + NRLib::ToString(lower) + ", " + NRLib::ToString(upper) + ").";
    throw Exception(err_txt);
  }

  if(IsZero(std_dev_)) {
    if(lower <= mean_ && mean_ <= upper)
      return mean_;
    else
      throw Exception("Can not draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }

  double out = DrawTwoSidedTruncatedAcceptReject(mean_, std_dev_, lower, upper);

  return(out);
}


double Normal::DrawTwoSidedTruncated(double lower, double upper, RandomGenerator &g) const
{
  if(lower > upper) {
    std::string err_txt = "Wrong input to truncated normal, should have min <= max). (min, max) = ("
      + NRLib::ToString(lower) + ", " + NRLib::ToString(upper) + ").";
    throw Exception(err_txt);
  }

  if(IsZero(std_dev_)) {
    if(lower <= mean_ && mean_ <= upper)
      return mean_;
    else
      throw Exception("Can not draw from a truncated distribution with zero variance, when the mean is outside the truncation values");
  }

  double out = DrawTwoSidedTruncatedAcceptReject(mean_, std_dev_, lower, upper, &g);

  return(out);
}


double Normal::DrawLowerTruncatedAcceptReject(double mu,
                                              double sigma,
                                              double lower,
                                              RandomGenerator *g) const
{
  double xL = (lower - mu)/sigma;
  double out(0.0);
  Normal norm01(0,1);
  double phi_lower = norm01.Cdf(xL);
  double phi_upper = 1.0;
  const double max_attempts(100000);

  int attempts(1);
  double ru01;
  if(g == NULL)
    ru01 = NRLib::Random::Unif01Open();
  else
    ru01 = g->Unif01Open();
  double ordinate = phi_lower + ru01 * (phi_upper - phi_lower);
  bool ok(false);

  if(phi_upper == phi_lower)
    out = lower;
  else {
    while(!ok && attempts < max_attempts) {
        out = norm01.PhiInverse(ordinate);
        out *= sigma;
        out += mu;
      if(out < lower) {
        attempts++;
        if(g == NULL)
          ru01 = NRLib::Random::Unif01Open();
        else
          ru01 = g->Unif01Open();
        ordinate = phi_lower + ru01 * (phi_upper - phi_lower);
      }
      else
        ok = true;
    }
    if(!ok) { //have tried many times without success
      double tm = mu + sigma * norm01.Pdf(xL)/(phi_upper - phi_lower);
      if(g == NULL)
        ru01 = NRLib::Random::Unif01Open();
      else
        ru01 = g->Unif01Open();
      //returns a uniform sample in (lower_, 2*tm-lower_)
      //where tm is the mean of the truncated normal
      out = lower + ru01*2*(tm - lower);
    }
  }

  return(out);
}

double Normal::DrawTwoSidedTruncatedAcceptReject(double mu,
                                                 double sigma,
                                                 double lower,
                                                 double upper,
                                                 RandomGenerator *g) const
{
  double xL = (lower - mu)/sigma;
  double xU =  (upper - mu)/sigma;
  Normal norm01(0,1);
  double phi_lower = norm01.Cdf(xL);
  double phi_upper = norm01.Cdf(xU);
  const double max_attempts(100000);

  int attempts(1);
  double ru01;
  if(g == NULL)
    ru01 = NRLib::Random::Unif01Open();
  else
    ru01 = g->Unif01Open();
  double ordinate = phi_lower + ru01 * (phi_upper - phi_lower);
  bool ok(false);

  double out(0.0);
  if(phi_upper == phi_lower)
    out = (upper + lower)/2.0;
  else {
    while(!ok && attempts < max_attempts) {
      if (ordinate == 1.0)
        out = upper;
      else {
        out = norm01.PhiInverse(ordinate);
        out *= sigma;
        out += mu;
      }
      if(out < lower || out > upper) {
        attempts++;
        if(g == NULL)
          ru01 = NRLib::Random::Unif01Open();
        else
          ru01 = g->Unif01Open();
        ordinate = phi_lower + ru01 * (phi_upper - phi_lower);
      }
      else {
        ok = true;
      }
    }
  }

  if (attempts == max_attempts)
    out = (upper + lower)/2.0;

  return out;
}

} // namespace NRLib
