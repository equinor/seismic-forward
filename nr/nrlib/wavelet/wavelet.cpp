// $Id: wavelet.cpp 1380 2016-10-31 13:13:39Z perroe $

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

#include "wavelet.hpp"

#include <cmath>
#include <string>

#include "waveletio.hpp"
#include "../iotools/stringtools.hpp"
#include "../math/constants.hpp"
#include "../math/mathutility.hpp"

namespace NRLib {

Wavelet::Wavelet()
  : wavelet_(0),
    dt_(1.0)
{}


Wavelet::Wavelet(size_t nw,
                 double dt)
  : wavelet_(nw, 0.0),
    dt_(dt)
{}


Wavelet::Wavelet(double peak_f,
                 size_t nw,
                 double dt,
                 double scaling_factor)
 : wavelet_(nw),
   dt_(dt)
{
  for (size_t i = 0; i < nw; ++i) {
    double t = GetTimeForIndex(i);
    double c = NRLib::Pi * NRLib::Pi * peak_f * peak_f * t * t * 1e-6;
    wavelet_[i] = (1 - 2*c) * std::exp(-c);
  }

  if (scaling_factor != 1.0)
    Scale(scaling_factor);
}


Wavelet::Wavelet(const std::string& filename,
                 WaveletFileType    file_type,
                 double             scaling_factor)
{
  if (file_type == UnknownWavelet) {
    file_type = FindWaveletFileType(filename);
  }
  if (file_type == UnknownWavelet) {
    throw FileFormatError(filename + " has unknown wavelet format.");
  }
  ReadFromFile(filename, file_type);

  if (scaling_factor != 1.0)
    Scale(scaling_factor);
}


void Wavelet::Resample(double new_dt)
{
  if (NRLib::IsEqual(new_dt, dt_))
    return;

  if (new_dt < dt_ || !NRLib::IsZero(std::fmod(new_dt, dt_)))
    throw NRLib::Exception("Resampling only supported if new dt is a multiple of original dt.");

  int subsampling = static_cast<int>(std::floor(new_dt / dt_ + 0.5));

  int old_center = static_cast<int>(GetCenterIndex());
  int old_nw = static_cast<int>(GetNSamples());
  int nw = 0;
  if (old_nw > 0)
    nw = 1 + old_center / subsampling + (old_nw - old_center - 1) / subsampling;
  std::vector<double> new_wavelet(nw, 0.0);

  int new_center = nw/2;

  for (int i = 0; i < nw; ++i) {
    int old_i = (i - new_center)*subsampling + old_center;
    new_wavelet[i] = wavelet_[old_i];
  }

  wavelet_.swap(new_wavelet);
  dt_ = new_dt;
}


void Wavelet::Resize(size_t nw)
{
  if (nw == wavelet_.size())
    return;

  std::vector<double> new_wavelet(nw, 0.0);
  size_t old_center = wavelet_.size() / 2;
  size_t new_center = nw/2;
  size_t start = 0;  // Start in old wavelet;
  size_t end   = wavelet_.size();
  size_t shift = 0;  // Start in new wavelet;
  if (old_center > new_center) {
    // Wavelet has become smaller.
    start = old_center - new_center;
    end   = nw + start;
  }
  else {
    shift = new_center - old_center;
  }

  for (size_t i = start; i < end; ++i) {
    new_wavelet[i + shift] = wavelet_[i - start];
  }

  wavelet_.swap(new_wavelet);
}


void Wavelet::ReadFromFile(const std::string& filename, WaveletFileType type)
{
  switch (type) {
  case JasonAscii:
    ReadJasonAsciiFile(filename, *this);
    break;
  case SyntoolOW:
    ReadSyntoolOWFile(filename, *this);
    break;
  case Strata:
    ReadStrataFile(filename, *this);
    break;
  case LandmarkAscii:
    ReadLandmarkAsciiFile(filename, *this);
    break;
  default:
    throw Exception("Unknown wavelet file format: " + ToString(type));
  }
}


void Wavelet::WriteToFile(const std::string              & filename,
                          WaveletFileType                  type,
                          const std::vector<std::string> & text_header)
{
  switch (type) {
  case JasonAscii:
    WriteJasonAsciiFile(filename, *this, text_header);
    break;
  default:
    throw Exception("Format not supported for wavelet output: " + ToString(type));
  }
}


void Wavelet::Scale(double scaling_factor)
{
  for (size_t i = 0; i < wavelet_.size(); ++i) {
    wavelet_[i] *= scaling_factor;
  }
}


} // namespace NRLib
