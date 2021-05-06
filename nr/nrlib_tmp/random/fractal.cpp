// $Id: fractal.cpp 1083 2012-10-01 10:55:34Z georgsen $

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


#include "fractal.hpp"

//#include "../iotools/stringtools.hpp"
//#include "../math/mathutility.hpp"

namespace NRLib {

Fractal::Fractal(double min, double max, double dim)
  : min_(min),
    max_(max),
    dim_(dim)
{
  min_fractal_ = std::pow(min_, -dim_);
  max_fractal_ = std::pow(max_, -dim_);
}

Fractal::~Fractal()
{}

Distribution<double>*
Fractal::Clone() const
{
  return new Fractal(*this);
}

double Fractal::DrawFromUpperQuantile(double qu,
                                      bool has_limit_info,
                                      double limit) const
{
  assert(0.0 <= qu && qu <= 1.0);
  double unif01 = NRLib::Random::Unif01();
  double d_max(max_fractal_);
  if (max_ < 0.0)
    d_max = 0.0;
  double d_min(min_fractal_);
  if (has_limit_info) {
    if (limit >= min_)
      d_min = std::pow(limit, -dim_);
  }

  double d = std::pow(d_max + unif01 * (d_min - d_max) * (1.0 - qu),
                      -1.0/dim_);
  return(d);
}


double Fractal::DrawFromLowerQuantile(double ql) const
{
  assert(0.0 <= ql && ql <= 1.0);
  double unif01 = NRLib::Random::Unif01();
  double d_max(max_fractal_);
  if (max_ < 0.0)
    d_max = 0.0;
  double d_min(min_fractal_);
  double d = std::pow(d_min - unif01 * (d_min - d_max) * ql,
                      -1.0/dim_);
  return(d);
}

double Fractal::DrawFromMiddleQuantile(double ql,
                                       double qu) const
{
  assert(0.0 <= ql && ql < qu  && qu <= 1.0);
  double unif01 = NRLib::Random::Unif01();
  double d_max(max_fractal_);
  if (max_ < 0.0)
    d_max = 0.0;
  double d_min(min_fractal_);
  double d = std::pow(d_min + (ql + unif01 * (qu - ql)) * (d_max-d_min),
                      -1.0/dim_);
  return(d);
}

} // namespace NRLib
