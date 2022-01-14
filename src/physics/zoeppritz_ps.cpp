// $Id: zoeppritz_ps.cpp 35 2013-10-22 10:57:21Z anner $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �    Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �    Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include <stdio.h>
#include "zoeppritz_ps.hpp"
#include <math.h>

ZoeppritzPS::ZoeppritzPS(void) {

}

ZoeppritzPS::~ZoeppritzPS() {

}

void ZoeppritzPS::ComputeConstants(double theta)
{
  a1_ = 0.0;
  sin2theta_ = sin(theta) * sin(theta);
  sin_theta_ = sin(theta);
  cos_theta_ = cos(theta);
}

double ZoeppritzPS::GetReflection(double diffvp, double meanvp, double diffrho, double meanrho, double diffvs, double meanvs)
{
  double cos_phi = cos(asin(meanvs * sin_theta_ / meanvp));
  a2_ = 2 * sin_theta_ * (meanvs * meanvs * sin2theta_ / (cos_phi * meanvp * meanvp) - meanvs * cos_theta_ / meanvp);
  a3_ = sin_theta_ * (-0.5 + meanvs * meanvs * sin2theta_ / (meanvp * meanvp) - meanvs * cos_theta_ * cos_phi / meanvp) / cos_phi;
  double refl = a2_ * diffvs / meanvs + a3_ * diffrho / meanrho;
  return refl;

}

