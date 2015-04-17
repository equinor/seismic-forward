// $Id: fftgrid2d.hpp 994 2012-02-29 15:15:04Z perroe $

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

#ifndef NRLIB_FFT_FFTGRID2D_HPP
#define NRLIB_FFT_FFTGRID2D_HPP

#include <cassert>
#include <complex>

// Must set MKL's /include/fftw or fftw's include directory as additional include directory.
#include "fftw3.h"

#include "../grid/grid2d.hpp"
#include "fft.hpp"

#ifdef FFTW_DEBUG
// Debugging
#include "../surface/regularsurface.hpp"
#include "../surface/surfaceio.hpp"
#include "../iotools/stringtools.hpp"
#endif

namespace NRLib {

template <typename T>
class FFTGrid2D
{
public:
  /// \param padding_ni    Suggested padding in i dimension.
  /// \param padding_nj    Suggested padding in j dimension.
  /// \param scale_forward If true scale both in forward and inverse transform.
  FFTGrid2D(size_t ni, size_t nj,  size_t padding_ni, size_t padding_nj, bool scale_forward);

  // Don't need this one yet. We'll have to worry about tapering later....
  //FFTGrid2D(const Grid2D<T>& grid, size_t padding_ni, size_t padding_nj, bool scale_forward);

  virtual ~FFTGrid2D();


  /// Initialize the real grid.
  /// \param values Input values. The size of the grid must be the same as the size of the
  ///               real value grid.
  void                     Initialize(const Grid2D<T> &values);

  /// Convolution with filter in FFT domain. Performs FFT for both this grid and filter.
  void                     Convolve(FFTGrid2D<T> & filter);

  /// Convolution with covariance filter. (Takes the square root).
  void                     ConvolveCovariance(FFTGrid2D<T> & filter);

  // Accessors in real and complex mode.
  // Uses the local, unpadded grid.
  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  T&                       Real(size_t i, size_t j)
  { return real_data_[GetRealIndex(i, j)]; }

  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  const T&                 Real(size_t i, size_t j) const
  { return real_data_[GetRealIndex(i, j)]; }

  // Uses the total grid.
  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  std::complex<T>&         Complex(size_t i, size_t j)
  { return complex_data_[GetComplexIndex(i, j)]; }

  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  const std::complex<T>&   Complex(size_t i, size_t j) const
  { return complex_data_[GetComplexIndex(i, j)]; }

  /// Get copy of grid. Only the original grid, no padding.
  Grid2D<T>                GetRealGrid() const;
  /// Get copy of the complete grid in Fourier domain.
  Grid2D<std::complex<T> > GetComplexGrid() const;

  size_t GetRealNI()    const { return ni_; }
  size_t GetRealNJ()    const { return nj_; }
  size_t GetComplexNI() const { return ni_tot_/2 + 1; }
  size_t GetComplexNJ() const { return nj_tot_; }
  size_t GetNItot()     const { return ni_tot_;}
  size_t GetNJtot()     const { return nj_tot_;}

  void DoFFT();
  void DoInverseFFT();

  // Data accessors. Do not use, except with in FFTW interface and similar.
  T*                     RealData()           { return real_data_; }
  const T*               RealData()    const  { return real_data_; }
  std::complex<T>*       ComplexData()        { return complex_data_; }
  const std::complex<T>* ComplexData() const  { return complex_data_; }

#ifdef FFTW_DEBUG
  // File output - for debugging
  void WriteRealToFile(const std::string& filename);
  void WriteComplexToFile(const std::string& filename);
#endif

private:
  /// True if we shall scale both forward and inverse transform. If false only the
  /// inverse transform is scaled.
  bool scale_forward_;

  // Grid sizes, main grid and total grid in real domain.
  size_t ni_;
  size_t nj_;
  size_t ni_tot_;
  size_t nj_tot_;

  /// Real data. Column-major ordering.
  T*               real_data_;
  /// Complex data. Column major ordering.
  std::complex<T>* complex_data_;

// Make copying illegal.
  FFTGrid2D(FFTGrid2D<T>& rhs);
  FFTGrid2D<T>& operator=(FFTGrid2D<T>& rhs);


  /// Find array index in the main real grid.
  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  inline size_t GetRealIndex(size_t i, size_t j) const;


  /// Find array index in the main real grid.
  /// \throw IndexOutOfRange when compiled in with DEBUG if i or j is out of range.
  inline size_t GetComplexIndex(size_t i, size_t j) const;
};

namespace NRLibPrivate {
  /// Computes 2D FFT transform. No scaling. Assumes column major storage.
  template <typename T>
  void ComputeFFT2D(size_t ni, size_t nj, T* in, std::complex<T>* out);

  /// Computes inverse 2D FFT transform. No scaling. Assumes column major storage.
  template <typename T>
  void ComputeFFT2DInverse(size_t ni, size_t nj, std::complex<T>* in, T* out);

  // Declare specialization - to ensure Linux compilation.
  template <>
  void ComputeFFT2D<float>(size_t ni, size_t nj, float* in, std::complex<float>* out);
  template <>
  void ComputeFFT2D<double>(size_t ni, size_t nj, double* in, std::complex<double>* out);
  template <>
  void ComputeFFT2DInverse<float>(size_t ni, size_t nj,
                                  std::complex<float>* in, float* out);
  template <>
  void ComputeFFT2DInverse<double>(size_t ni, size_t nj,
                                   std::complex<double>* in, double* out);

}


// -----------------  TEMPLATE IMPLEMENTATION ---------------------

template <typename T>
FFTGrid2D<T>::FFTGrid2D(size_t ni, size_t nj,  size_t padding_ni,
                        size_t padding_nj, bool scale_forward)
  : scale_forward_(scale_forward), ni_(ni), nj_(nj)
{
  // Find total sizes.
  ni_tot_ = FindNewSizeWithPadding(ni + padding_ni, true);
  nj_tot_ = FindNewSizeWithPadding(nj + padding_nj);

  // Allocate aligned data for efficiency.
  // Add padding for use in inplace transform.
  real_data_ = reinterpret_cast<T*>(fftw_malloc(ni_tot_ * nj_tot_ * sizeof(T)));

  size_t n_complex = GetComplexNI() * GetComplexNJ();
  size_t complex_datalen = n_complex * sizeof(std::complex<T>);
  complex_data_ = reinterpret_cast<std::complex<T>*>(fftw_malloc(complex_datalen));

  for (size_t i = 0; i < ni_tot_*nj_tot_; ++i) {
    real_data_[i] = 0;
  }

  for (size_t i = 0; i < n_complex; ++i) {
    complex_data_[i] = 0;
  }
}


template <typename T>
FFTGrid2D<T>::~FFTGrid2D()
{
  fftw_free(real_data_);
  fftw_free(complex_data_);
}


template <typename T>
void FFTGrid2D<T>::DoFFT() {
  if (scale_forward_) {
    size_t n = ni_tot_ * nj_tot_;
    double scale = 1.0 / sqrt(1.0*n);
    for (size_t i = 0; i < n; ++i) {
      real_data_[i] *= static_cast<T>(scale);
    }
  }
  NRLibPrivate::ComputeFFT2D(ni_tot_, nj_tot_, real_data_, complex_data_);
}

template <typename T>
void FFTGrid2D<T>::DoInverseFFT() {
  NRLibPrivate::ComputeFFT2DInverse(ni_tot_, nj_tot_, complex_data_, real_data_);
  double scale;
  size_t n = ni_tot_ * nj_tot_;
  if (!scale_forward_)
    scale=1.0 / n;
  else
    scale=1.0 / sqrt(1.0*n);

  for (size_t i = 0; i < n; ++i) {
    real_data_[i] *= static_cast<T>(scale);
  }
}


template <typename T>
size_t FFTGrid2D<T>::GetRealIndex(size_t i, size_t j) const
{
  assert(i < GetRealNI() && j < GetRealNJ());
  return i + j*ni_tot_;
}


template <typename T>
size_t FFTGrid2D<T>::GetComplexIndex(size_t i, size_t j) const
{
  assert(i < GetComplexNI() && j < GetComplexNJ());
  return i + j*GetComplexNI();
}


template <typename T>
void FFTGrid2D<T>::Initialize(const Grid2D<T> &values)
{
  assert(ni_tot_ == values.GetNI() && nj_tot_ == values.GetNJ());

  size_t k = 0;
  for (size_t j = 0; j < values.GetNJ(); j++)
    for (size_t i = 0; i < values.GetNI(); i++) {
      real_data_[k] = values(i,j);
      k++;
    }
}


template <typename T>
void FFTGrid2D<T>::Convolve(FFTGrid2D<T> &filter)
{
  DoFFT();
  filter.DoFFT();
  for(size_t i = 0; i < GetComplexNI()*GetComplexNJ(); i++)
  {
    complex_data_[i] *= filter.complex_data_[i];
  }
  DoInverseFFT();
}


template <typename T>
void FFTGrid2D<T>::ConvolveCovariance(FFTGrid2D<T> &filter)
{
  DoFFT();
  filter.DoFFT();
  for(size_t i = 0; i < GetComplexNI()*GetComplexNJ(); i++)
  {
    complex_data_[i] *= std::sqrt(filter.complex_data_[i]);
  }
  DoInverseFFT();
}


template <typename T>
Grid2D<T>  FFTGrid2D<T>::GetRealGrid() const
{
  size_t i,j;
  Grid2D<T> output(ni_,nj_);
  for(i=0;i<ni_;i++)
    for(j=0;j<nj_;j++)
      output(i,j) = real_data_[i+ni_tot_*j];

  return output;

}

#ifdef FFTW_DEBUG
template <typename T>
void FFTGrid2D<T>::WriteRealToFile(const std::string& filename)
{
  RegularSurface<double> real_surf(0.0, 0.0, ni_tot_, nj_tot_, ni_tot_, nj_tot_);
  for (size_t i = 0; i < ni_tot_; ++i) {
    for (size_t j = 0; j < nj_tot_; ++j) {
      real_surf(i, j) = static_cast<double>(real_data_[i + j*ni_tot_]);
    }
  }
  WriteIrapClassicAsciiSurf(real_surf, filename);
}


template <typename T>
void FFTGrid2D<T>::WriteComplexToFile(const std::string& filename)
{
  RegularSurface<double> real_surf(0.0, 0.0, ni_tot_/2 + 1, nj_tot_,
                                   ni_tot_/2 + 1, nj_tot_);
  for (size_t i = 0; i < GetComplexNI(); ++i) {
    for (size_t j = 0; j < nj_tot_; ++j) {
      real_surf(i, j) = static_cast<double>(complex_data_[i + j*GetComplexNI()].real());
    }
  }
  WriteIrapClassicAsciiSurf(real_surf, "real_" + filename);
  for (size_t i = 0; i < GetComplexNI(); ++i) {
    for (size_t j = 0; j < nj_tot_; ++j) {
      real_surf(i, j) = static_cast<double>(complex_data_[i + j*GetComplexNI()].imag());
    }
  }
  WriteIrapClassicAsciiSurf(real_surf, "imag_" + filename);
}
#endif


// ----------------  FFTW interface functions ---------------------

// Non-specialized versions, should give error.

template <typename T>
void NRLibPrivate::ComputeFFT2D(size_t ni, size_t nj,
                                T* in, std::complex<T>* out)
{
  // 2D FFT computations are only implemented for float and double.
  assert(false);
}

template <typename T>
void NRLibPrivate::ComputeFFT2DInverse(size_t ni, size_t nj,
                                       std::complex<T>* in, T* out)
{
  // 2D FFT computations are only implemented for float and double.
  assert(false);
}


} // namespace NRLib


#endif // NRLIB_FFT_FFT2D_HPP
