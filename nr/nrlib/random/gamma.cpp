// $Id: gamma.cpp 883 2011-09-26 09:17:05Z perroe $

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

#include "gamma.hpp"
#include "../exception/exception.hpp"
#include <limits>
#include <cmath>

// Functions follow Numerical Recipes by Press, Teukolsky, Vetterling and Flannery


namespace NRLib {
  static const int    ASWITCH = 100;
  static const int    ngau    = 18;
  static const double EPS     = std::numeric_limits<double>::epsilon();
  static const double FPMIN   = std::numeric_limits<double>::min()/EPS;
  double              gln;

  const double y[18] = {0.0021695375159141994,
  0.011413521097787704,0.027972308950302116,0.051727015600492421,
  0.082502225484340941, 0.12007019910960293,0.16415283300752470,
  0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
  0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
  0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
  0.87126389619061517, 0.95698180152629142};

  static const double w[18] = {0.0055657196642445571,
  0.012915947284065419,0.020181515297735382,0.027298621498568734,
  0.034213810770299537,0.040875750923643261,0.047235083490265582,
  0.053244713977759692,0.058860144245324798,0.064039797355015485,
  0.068745323835736408,0.072941885005653087,0.076598410645870640,
  0.079687828912071670,0.082187266704339706,0.084078218979661945,
  0.085346685739338721,0.085983275670394821};

  double gammln(const double xx)
  {
    // Returns the value ln(Gamma(xx)) for xx>0
    int    j;
    double x;
    double tmp;
    double y;
    double ser;
    static const double cof[14]={57.1562356658629235,-59.5979603554754912,
      14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
      .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
      -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
      .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

    if (xx <= 0)
      throw Exception("Bad argument in gammln");

    y=x=xx;
    tmp = x+5.24218750000000000;
    tmp =(x+0.5)*std::log(tmp)-tmp;
    ser = 0.999999999999997092;
    for (j=0; j<14; j++)
    {
      y++;
      ser += cof[j]/y;
    }
    return tmp+log(2.5066282746310005*ser/x);
  }

  double gammp(const double a, const double x)
  {
    // Returns incomplete gamma function P(a,x)
    if (x < 0 || a <= 0)
      throw Exception("Bad arguments");
    if (x==0)
      return 0;
    else if (a >= ASWITCH)
      return gammpapprox(a,x,1);
    else if (x < a+1)
      return gser(a,x);
    else
      return 1-gcf(a,x);
  }

  double gser(const double a, const double x)
  {
    // Returns the incomplete gamma function P(a,x) evaluated by its series representation. Also sets ln gamma(a) as gln.
    double sum;
    double del;
    double ap;

    gln = gammln(a);
    ap  = a;
    del = 1/a;
    sum = del;
    for ( ; ; )
    {
      ap++;
      del *= x/ap;
      sum += del;
      if (std::abs(del) < std::abs(sum)*EPS)
        return sum*std::exp(-x+a*std::log(x)-gln);
    }
  }

  double gcf(const double a, const double x)
  {
    //Returns the incomplete gamma function Q(a,x) evaluated by its continued fraction representation. Also sets ln Gamma(a) as gln.
    int    i;
    double an,b,c,d,del,h;

    gln = gammln(a);
    b   = x+1-a;
    c   = 1/FPMIN;
    d   = 1/b;
    h   = d;

    for (i=1;;i++)
    {
      an = -i*(i-a);
      b += 2;
      d  = an*d+b;

      if (std::abs(d) < FPMIN )
        d = FPMIN;

      c = b+an/c;

      if (std::abs(c) < FPMIN)
        c = FPMIN;

      d   = 1/d;
      del = d*c;
      h  *= del;

      if (std::abs(del-1) <= EPS)
        break;
    }
    return std::exp(-x+a*std::log(x)-gln)*h;
  }

  double gammpapprox(double a, double x, int psig)
  {
    // Incomplete gamma by quadrature. Returns P(a,x) or Q(a,x) when psig is 1 or 0, respectively
    int    j;
    double xu,t,sum,ans;

    double a1     = a-1;
    double lna1   = std::log(a1);
    double sqrta1 = std::sqrt(a1);

    gln = gammln(a);

    if (x > a1)
      xu = std::max(a1+11.5*sqrta1, x+6*sqrta1);
    else
    {
      double min = std::min(a1-7.5*sqrta1, x-5*sqrta1);
      xu = std::max(0.0, min);
    }

    sum = 0;

    for (j=0; j<ngau; j++)
    {
      t    = x+(xu-x)*y[j];
      sum += w[j]*std::exp(-(t-a1)+a1*(std::log(t)-lna1));
    }
    ans = sum*(xu-x)*std::exp(a1*(lna1-1)-gln);

    if( psig == 1 )
    {
      if (ans > 0 )
        return 1-ans;
      else
        return -ans;
    }
    else
    {
      if (ans >= 0)
        return ans;
      else
        return 1+ans;
    }
  }

  double invgammp(double p, double a)
  {
    int    j;
    double x,err,t,u,pp;
    double lna1      = 0.0;
    double afac      = 0.0;
    double a1        = a-1;
    const double EPS = 1.e-8;

    gln = gammln(a);
    if (a <= 0.0)
      throw Exception("a must be positive in invgammap");
    if (p >= 1.0)
      return std::max(100.0,a + 100.0*std::sqrt(a));
    if (p <= 0.0)
      return 0.0;
    if (a > 1.0)
    {
      lna1 = std::log(a1);
      afac = std::exp(a1*(lna1-1.0)-gln);
      if (p < 0.5)
        pp = p;
      else
        pp = 1.0 - p;
      t = std::sqrt(-2.0*std::log(pp));
      x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
      if (p < 0.5)
        x = -x;
      x = std::max(1.e-3,a*std::pow(1.0-1.0/(9.0*a)-x/(3.0*sqrt(a)),3));
    }
    else
    {
      t = 1.0 - a*(0.253+a*0.12);
      if (p < t)
        x = std::pow(p/t,1.0/a);
      else
        x = 1.0-std::log(1.0-(p-t)/(1.0-t));
    }
    for (j=0; j<12; j++)
    {
      if (x <= 0.0)
        return 0.0;
      err = gammp(a,x) - p;
      if (a > 1.0)
        t = afac*std::exp(-(x-a1)+a1*(std::log(x)-lna1));
      else
        t = std::exp(-x+a1*std::log(x)-gln);
      u = err/t;
      x -= (t=u/(1.0-0.5*std::min(1.0,u*((a-1.0)/x - 1))));
      if (x <= 0.0)
        x = 0.5*(x + t);
      if (std::abs(t) < EPS*x )
        break;
    }
    return x;
  }

  // Logarithm of binomial coefficient
  double bico(int n, int k)
  {
    if(n<0 || k<0 || k>n)
      throw Exception("Bad args in bico");
    if(n==k || k==0)
      return 1.0;
    else
    {
      double result = floor(0.5+std::exp(gammln(double(n+1.))-gammln(double(k+1.))-gammln(double(n-k+1.0))));
      return(result);
    }

  }

} // namespace NRLib
