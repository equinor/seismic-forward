// $Id: eclipsetransmult.hpp 1465 2017-04-07 13:28:48Z perroe $

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

#ifndef NRLIB_ECLIPSEGRID_ECLIPSETRANSMULT_HPP
#define NRLIB_ECLIPSEGRID_ECLIPSETRANSMULT_HPP

#include "eclipsetools.hpp"

#include <fstream>
#include <list>

namespace NRLib {

class EclipseTransMult{
public:
  EclipseTransMult(int ix1, int ix2, int jy1, int jy2, int kz1, int  kz2,
                   double value, NRLib::Face & face);
  ~EclipseTransMult() {};

  static void WriteMultiply(std::ofstream                     & out_file,
                            const std::list<EclipseTransMult> & trans_mult);

private:
  int ix1_, ix2_, jy1_, jy2_, kz1_, kz2_;
  double multiplier_;

  NRLib::Face face_type_;
};

}

#endif
