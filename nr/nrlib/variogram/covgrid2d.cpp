
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

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdarg.h>
#include <ostream>

#include "nrlib/exception/exception.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "covgrid2d.hpp"
#include "variogram.hpp"

using namespace NRLib;

CovGrid2D::CovGrid2D(const Variogram & variogram,
                     int      nx,
                     int      ny,
                     double   dx,
                     double   dy)
  : cov_(2*nx*ny),
    nx_(nx),
    ny_(ny),
    dx_(dx),
    dy_(dy)
{
  for(int i=0;i<nx;i++) {
    for(int j=-ny;j<ny;j++) {
      int index = i*2*ny + ny + j;
      cov_[index] = static_cast<float>(variogram.GetCorr(static_cast<double>(i*dx), static_cast<double>(j*dy)));
    }
  }
}

float
CovGrid2D::GetCov(int deltai, int deltaj) const
{
  if(deltai<0) {
    deltai = -deltai;
    deltaj = -deltaj;
  }
  int index = deltai*2*ny_ + deltaj + ny_;
  return cov_[index];
}

double
CovGrid2D::GetRangeX() const
{
  //
  // Find range in X direction: max distance where cov >= 0.05
  //
  double rx = 0;
  for(int i=0;i<nx_;i++) {
    for(int j=-ny_;j<ny_;j++) {
      int index = i*2*ny_ + ny_ + j;
      if (cov_[index] >= 0.05) {
        if (i*dx_ > rx)
          rx = i*dx_;
      }
    }
  }
  return rx;
}

double
CovGrid2D::GetRangeY() const
{
  //
  // Find range in Y direction: max distance where cov >= 0.05
  //
  double ry = 0;
  for(int i=0;i<nx_;i++) {
    for(int j=-ny_;j<ny_;j++) {
      int index = i*2*ny_ + ny_ + j;
      if (cov_[index] >= 0.05) {
        if (j*dy_ > ry)
          ry = j*dy_;
      }
    }
  }
  return ry;
}

void
CovGrid2D::WriteToFile(const std::string & name) const
{
  //
  // Write grid using an ASCII Irap Classic surface format
  //
  std::ofstream file;
  NRLib::OpenWrite(file, name);
  file.precision(14);
  file << 2*nx_-1  << " "
       << 2*ny_-1  << " "
       << dx_      << " "
       << dy_      << "\n"
       << -dx_*nx_ << " "
       <<  dx_*nx_ << " "
       << -dy_*ny_ << " "
       <<  dy_*ny_ << "\n";

  for(int j=-ny_+1;j<ny_;j++)
    for(int i=-nx_+1;i<nx_;i++)
      file << GetCov(i,j) << " ";
  file << std::endl;
  file.close();
}

