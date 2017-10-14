// $Id: eclipsetransmult.cpp 1472 2017-04-27 11:41:04Z eyaker $

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

#include "eclipsetransmult.hpp"

#include "../iotools/stringtools.hpp"

namespace NRLib{

  EclipseTransMult::EclipseTransMult(int ix1, int ix2, int jy1, int jy2, int kz1, int kz2,
                                     double value, NRLib::Face & face)
    :ix1_(ix1),
    ix2_(ix2),
    jy1_(jy1),
    jy2_(jy2),
    kz1_(kz1),
    kz2_(kz2),
    multiplier_(value),
    face_type_(face)
  {}

void EclipseTransMult::WriteMultiply(std::ofstream                     & out_file,
                                     const std::list<EclipseTransMult> & trans_mult)
{
  if (trans_mult.empty()) {
    return;
  }

  out_file << "MULTIPLY\n";
  out_file << std::setprecision(10);

  std::list<EclipseTransMult>::const_iterator it;
  for (it = trans_mult.begin(); it != trans_mult.end(); ++it){
    int ix = it->ix1_;   // same as ix2_
    int jy = it->jy1_;   // same as jy2_
    int kz = it->kz1_;   // same as kz2_

    if (it->face_type_ == PosX){
      out_file << "'TRANX'" << std::setw(15);
    }
    else if (it->face_type_ == NegX){
      out_file << "'TRANX'" << std::setw(15);
      ix = ix - 1;
    }
    if (it->face_type_ == PosY){
      out_file << "'TRANY'" << std::setw(15);
    }
    else if (it->face_type_ == NegY){
      out_file << "'TRANY'" << std::setw(15);
      jy = jy - 1;
    }
    if (it->face_type_ == PosZ){
      out_file << "'TRANZ'" << std::setw(15);
    }
    else if(it->face_type_ == NegZ){
      out_file << "'TRANZ'" << std::setw(15);
      kz = kz - 1;
    }
    out_file << it->multiplier_ << std::setw(7);
    out_file << ix << std::setw(7) << ix << std::setw(7)
             << jy << std::setw(7) << jy << std::setw(7)
             << kz << std::setw(7) << kz << std::setw(5);

    out_file << "/\n";
  }
  out_file << "/\n\n";
}



} // namespace NRLib

