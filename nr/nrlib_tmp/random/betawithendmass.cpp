// $Id: betawithendmass.cpp 1146 2013-04-08 10:33:13Z veralh $

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

#include "betawithendmass.hpp"
#include "beta.hpp"
#include <cmath>
#include <cstdio>
#include <cassert>

namespace NRLib {

BetaWithEndMass::BetaWithEndMass(double min_val, double max_val, double p_min, double p_max, double a, double b, int nx)
{
  assert(a > 0.0);
  assert(b > 0.0);
  assert(max_val > min_val);
  assert(p_min >= 0 && p_min <= 1);
  assert(p_max >= 0 && p_max <= 1);
  p_min_ = p_min;
  p_max_ = p_max;
  nx_    = nx;

  beta_ = NRLib::Beta(min_val, max_val, a, b, nx);
  //Implicitt call of beta_.ComputePdfAndCdf() and beta_.ComputeQuantile();
}

BetaWithEndMass::BetaWithEndMass()
{
}

BetaWithEndMass::~BetaWithEndMass()
{
}

Distribution<double>*
BetaWithEndMass::Clone() const
{
  return new BetaWithEndMass(*this);
}

} // namespace NRLib
