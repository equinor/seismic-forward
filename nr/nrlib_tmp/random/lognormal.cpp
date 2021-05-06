// $Id: lognormal.cpp 1083 2012-10-01 10:55:34Z georgsen $

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


#include "lognormal.hpp"

#include "../iotools/stringtools.hpp"
#include "../math/mathutility.hpp"

namespace NRLib {


LogNormal::LogNormal()
  : normal_(Normal()),
    mu_(0.0),
    sigma_(0.0)
{
  exp_ = std::exp(0.5);
  std_dev_ = std::sqrt((std::exp(1.0) - 1.0) * std::exp(1.0));
}

LogNormal::LogNormal(double mu, double sigma)
  : normal_(Normal(mu, sigma)),
    mu_(mu),
    sigma_(sigma)
{
  exp_ = std::exp(mu_ + 0.5 * sigma_ * sigma_);
  std_dev_ = std::sqrt((std::exp(sigma_*sigma_)-1.0) * std::exp(2.0 * mu_ + sigma_ * sigma_));
}

LogNormal::~LogNormal()
{
}

Distribution<double>*
LogNormal::Clone() const
{
  return new LogNormal(*this);
}

}
