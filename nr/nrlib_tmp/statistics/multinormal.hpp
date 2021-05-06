// $Id: multinormal.hpp 1241 2014-02-18 12:17:26Z anner $

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


#ifndef NRLIB_STATISTICS_MULTINORMAL_HPP
#define NRLIB_STATISTICS_MULTINORMAL_HPP

#include "../math/constants.hpp"
#include "../exception/exception.hpp"
#include "../random/randomgenerator.hpp"
#include <cstdlib>
#include <vector>
#include "../flens/nrlib_flens.hpp"


namespace NRLib {

  class MultiNormal {

  public:
    MultiNormal();
    MultiNormal(const Vector    & exp,
                SymmetricMatrix   cov);
    ~MultiNormal();
    inline void SetParameters(const Vector & exp,
                              SymmetricMatrix cov);
    Vector Draw() const;
    Vector Draw(RandomGenerator & g) const;
    double Potential(const Vector & x) const; // -Log of pdf

  private:
    Vector exp_;
    SymmetricMatrix chol_fac_cov_;
    SymmetricMatrix cov_;
  };


void MultiNormal::SetParameters(const Vector & exp,
                              SymmetricMatrix cov)
{
  exp_ = exp;
  cov_ = cov;

}


}

#endif
