// $Id: eclipseeditnnc.cpp 1472 2017-04-27 11:41:04Z eyaker $

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

#include "eclipseeditnnc.hpp"

#include "../iotools/stringtools.hpp"

namespace NRLib{

EclipseEditNNC::EclipseEditNNC(int ix, int iy, int iz, int jx, int jy, int jz, double value)
  : ix_(ix),
    iy_(iy),
    iz_(iz),
    jx_(jx),
    jy_(jy),
    jz_(jz),
    tran_(value)
{}


void EclipseEditNNC::WriteEditNNC(std::ofstream                   & out_file,
                                  const std::list<EclipseEditNNC> & editNNC)
{
    if (editNNC.empty())
    return;

  out_file << "EDITNNC\n";
  out_file << std::setprecision(10);

  out_file << "--" << std::setw(4)
    << "IX" << std::setw(7) << "IY" << std::setw(7)
    << "IZ" << std::setw(7) << "JX" << std::setw(7)
    << "JY" << std::setw(7) << "JZ" << std::setw(15)
    << "TRANM" << "\n\n";

  std::list<EclipseEditNNC>::const_iterator it;
  for (it = editNNC.begin(); it != editNNC.end(); ++it) {
    out_file << std::setw(5);
    out_file << it->ix_ << std::setw(7);
    out_file << it->iy_ << std::setw(7);
    out_file << it->iz_ << std::setw(7);
    out_file << it->jx_ << std::setw(7);
    out_file << it->jy_ << std::setw(7);
    out_file << it->jz_ << std::setw(18);
    out_file << it->tran_ << std::setw(10) << "/\n";
  }

  out_file << "/\n\n";
}

}
