// $Id: wavelet.hpp 1380 2016-10-31 13:13:39Z perroe $

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

#ifndef NRLIB_WAVELET_WAVELET_HPP
#define NRLIB_WAVELET_WAVELET_HPP

#include <cassert>
#include <string>
#include <vector>

#include "waveletio.hpp"

namespace NRLib {

class Wavelet {
public:
  typedef std::vector<double>::iterator        iterator;
  typedef std::vector<double>::const_iterator  const_iterator;

  /// Default constructor - constructs empty wavelet.
  Wavelet();

  /// Empty wavelet of given length
  Wavelet(size_t nw, double dt);

  /// Ricker wavelet.
  Wavelet(double peak_frequency, size_t nw, double dt, double scaling_factor = 1.0);

  /// Wavelet from file.
  /// Type of file is automatically deduced if type is UnknownWavelet.
  Wavelet(const std::string & filename,
          WaveletFileType     type           = UnknownWavelet,
          double              scaling_factor = 1.0);

  double& operator[](size_t i)
  { assert(i < wavelet_.size()); return wavelet_[i]; }

  double  operator[](size_t i) const
  { assert(i < wavelet_.size()); return wavelet_[i]; }

  iterator       begin()       { return wavelet_.begin(); }
  iterator       end()         { return wavelet_.end(); }

  const_iterator begin() const { return wavelet_.begin(); }
  const_iterator end()   const { return wavelet_.end(); }

  size_t GetNSamples() const
  { return wavelet_.size(); }

  /// Time for wavelet sample. Time for center index is 0.
  double GetTimeForIndex(size_t i) const
  { assert(i < wavelet_.size()); return (static_cast<int>(i) - static_cast<int>(GetCenterIndex())) * dt_; }

  size_t GetCenterIndex() const                  { return wavelet_.size()/2; }

  double GetSampleInterval() const               { return dt_; }
  void   SetSampleInterval(double dt)            { dt_ = dt; }

  /// Resample to new sampling interval.
  void Resample(double dt);

  /// Resizes. Keeps the data. The center-index is kept at nw/2.
  /// We don't do any tapering if wavelets are shrunk.
  void   Resize(size_t nw);

  void ReadFromFile(const std::string& filename, WaveletFileType type);

  void WriteToFile(const std::string              & filename,
                   WaveletFileType                  type,
                   const std::vector<std::string> & text_header = std::vector<std::string>(0));

  void Scale(double scaling_factor);

private:
  /// Wavelet data
  /// The index for t=0 lies at wavelet_.size()/2.
  std::vector<double> wavelet_;

  /// Sample interval.
  double dt_;
};

}

#endif // NRLIB_WAVELET_WAVELET_HPP
