// $Id: fft.hpp 1006 2012-03-22 15:07:24Z hauge $

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

#ifndef NRLIB_FFT_HPP
#define NRLIB_FFT_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iterator>
#include <vector>

#include "fftw3.h"

// The "public" declarations
namespace NRLib {
  /// If minPadSize < 0 : sets the minimum padding size to distance between begin and end.
  /// If minPadSize == 0: No padding is added.
  /// *begin and *end are of type double
  /// container could be a std::vector or std::list of doubles
  template<class FI>
  void ComputeFFT1D(FI begin, FI end, std::vector<std::complex<double> >& vOut, bool scale_forward, int minPadSize = -1);
  template<class container>
  void ComputeFFT1D(const container& vIn, std::vector<std::complex<double> >& vOut, bool scale_forward, int minPadSize = -1);

  /// Padding is assumed to be : Padding = vIn.size() - distance(begin, end).
  /// *begin and *end are of type double
  /// container could be a std::vector or std::list of doubles
  template<class FI>
  void ComputeFFTInv1D(const std::vector<std::complex<double> >& vIn, FI begin, FI end, bool scale_forward);
  template<class container>
  void ComputeFFTInv1D(const std::vector<std::complex<double> >& vIn, container& vOut, bool scale_forward);

  size_t FindNewSizeWithPadding(size_t minSize, bool must_be_even = false);
}

// The "private" declarations
namespace NRLib {
  namespace NRLibPrivate {
    template <typename T> void ComputeFFT1D(size_t size, T *in, std::complex<T> * out);
    template <typename T> void ComputeFFTInv1D(size_t size, std::complex<T> * in, T *out);
     // Declare specialization - to ensure Linux compilation.
    template <> void ComputeFFT1D(size_t size, double *in, std::complex<double> * out);
    template <> void ComputeFFTInv1D(size_t size, std::complex<double> * in, double *out);
    template<class FI> void AddPadding(FI begin, FI end, double *v, size_t padSize, size_t sizeOrig, double scale_factor);
  }
}

/// Implementation "public"
template<class container>
inline void NRLib::ComputeFFT1D(const container& vIn, std::vector<std::complex<double> >& vOut, bool scale_forward, int minPadSize) {
  return ComputeFFT1D(vIn.begin(), vIn.end(), vOut, scale_forward, minPadSize);
}

template<class container>
inline void NRLib::ComputeFFTInv1D(const std::vector<std::complex<double> >& vIn, container& vOut, bool scale_forward) {
  return ComputeFFTInv1D(vIn, vOut.begin(), vOut.end(), scale_forward);
}

template<class FI>
void NRLib::ComputeFFT1D(FI begin, FI end, std::vector<std::complex<double> >& vOut, bool scale_forward, int minPadSize) {
  if (begin == end)
    return; // No input data, nothing to do.

  //Copy values in vector to the correct internal data structures
  size_t orig_size = distance(begin, end);
  size_t pad_size = (minPadSize < 0 ? orig_size : static_cast<size_t>(minPadSize));
  if (pad_size)
    pad_size = FindNewSizeWithPadding(pad_size + orig_size) - orig_size;

  size_t tot_size = orig_size + pad_size;
  double * real_data = reinterpret_cast<double *> (fftw_malloc(tot_size * sizeof(double)));

  FI itLoop = begin;
  double scale_factor(1.0);
  if (scale_forward)
    scale_factor = 1.0 / sqrt(1.0 * tot_size);
  for (size_t i = 0; i < orig_size; i++)
    real_data[i] = scale_factor  * (*itLoop++);

  size_t complex_size = (tot_size / 2) + 1;
  std::complex<double> * complex_data = reinterpret_cast<std::complex<double> *> (fftw_malloc(complex_size * sizeof(std::complex<double>)));

  NRLibPrivate::AddPadding(begin, end, real_data, pad_size, orig_size, scale_factor);
  NRLibPrivate::ComputeFFT1D(tot_size, real_data, complex_data);

  ////copy result into vector out
  vOut.resize(tot_size);
  vOut[0] = complex_data[0];

  for (size_t i = 1; i < complex_size; i++) {
    vOut[i] = complex_data[i];;
    vOut[tot_size - i] = complex_data[i];
  }
  fftw_free(real_data);
  fftw_free(complex_data);
}

template<class FI>
void
NRLib::ComputeFFTInv1D(const std::vector<std::complex<double> >& vIn, FI begin, FI end, bool scale_forward) {
  using namespace std;
  if (vIn.empty())
    return; // No data, nothing to do.

  size_t tot_size = vIn.size();
  size_t data_size = distance(begin, end);

  //Copy values in vector to the correct internal data structures
  std::complex<double> * complex_data = reinterpret_cast<std::complex<double> *> (fftw_malloc(tot_size * sizeof(std::complex<double>)));
  for (size_t i = 0; i < tot_size; i++)
    complex_data[i] = vIn[i];

  double * real_data = reinterpret_cast<double *> (fftw_malloc(tot_size * sizeof(double)));

  NRLibPrivate::ComputeFFTInv1D(tot_size, complex_data, real_data);

  double scale_factor(1.0 / tot_size);
  if (scale_forward)
    scale_factor = 1.0 / sqrt(1.0*tot_size);
  //copy result into vector out
  for (size_t i = 0; i < data_size; i++) {
    *begin++ = real_data[i] * scale_factor;
  }
  fftw_free(real_data);
  fftw_free(complex_data);
}

/// Implementation "private"
//
template<class FI>
void
NRLib::NRLibPrivate::AddPadding(FI begin, FI end, double *v,
                                size_t padSize, size_t sizeOrig,
                                double scale_factor) {
  double endv = *begin;
  if (padSize == 1) {
    v[sizeOrig] = scale_factor*endv;
  }
  else {
    FI endMi = --end;
    size_t i;
    for (i = 0; i < padSize; i++) {
      size_t j = sizeOrig + i;
      double t = static_cast<double>(i)/(padSize - 1);
      v[j] = scale_factor*(endv*t + (1 - t)*(*endMi));
    }
  }

}

#endif
