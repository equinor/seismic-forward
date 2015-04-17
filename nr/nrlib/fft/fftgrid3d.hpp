// $Id: fftgrid3d.hpp 1192 2013-08-16 14:48:34Z ok $

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

#ifndef NRLIB_FFT_FFTGRID3D_HPP
#define NRLIB_FFT_FFTGRID3D_HPP

#include <cassert>
#include <complex>

// Must set MKL's /include/fftw or fftw's include directory as additional include directory.
#include "fftw3.h"


// #define  FFTW_DEBUG
#ifdef FFTW_DEBUG
#include <string>
#include "../iotools/fileio.hpp"
#include <iostream>
#include <fstream>
#endif

#include "../grid/grid.hpp"
#include "fft.hpp"

namespace NRLib {

template <typename T>
class FFTGrid3D
{
public:
  /// \param padding_ni    Suggested padding in i dimension.
  /// \param padding_nj    Suggested padding in j dimension.
  /// \param padding_nk    Suggested padding in k dimension.
  /// \param scale_forward If true scale both in forward and inverse transform.
  FFTGrid3D(size_t ni, size_t nj, size_t nk, size_t padding_ni, size_t padding_nj, size_t padding_nk, bool scale_forward);

  virtual ~FFTGrid3D();


  /// Initialize the real grid.
  /// \param values Input values. The size of the grid must be the same as the size of the
  ///               real value grid.
  void                     Initialize(const Grid<T>& values);

  /// Initialize the active part of the grid with the input grid. The rest is padded with zeros.
  /// \param values Input values. The size of the grid must be the same as the size of the
  ///               real value grid.
  void                     InitializePadZero(const Grid<T> &values);

  /// Initialize all cells the real grid to a constant value.
  /// \param value Input value.
  void                     InitializeConstant(const T& value);

  /// Convolution with filter in FFT domain. Performs FFT for both this grid and filter.
  void                     Convolve(FFTGrid3D<T> & filter);

  /// Convolution with covariance filter. (Takes the square root).
  void                     ConvolveCovariance(FFTGrid3D<T> & filter);

  /// Convolution with filter in FFT domain. Does not take FFT, so only multiplies in
  /// the Fourier domain.
  void                     ConvolveNoFFT(const FFTGrid3D<T> & filter);

  // Accessors in real and complex mode.
  // Uses the local, unpadded grid.
  T&                       Real(size_t i, size_t j, size_t k)
  { return real_data_[GetRealIndex(i, j, k)]; }

  const T&                 Real(size_t i, size_t j, size_t k) const
  { return real_data_[GetRealIndex(i, j, k)]; }

  // Cyclic accessors. Negative i means nt_tot - i, same for j and k.
  T&                       RealCyclic(int i, int j, int k)
  { return real_data_[GetRealIndexCyclic(i, j, k)]; }

  const T&                 RealCyclic(int i, int j, int k) const
  { return real_data_[GetRealIndexCyclic(i, j, k)]; }

  // Uses the total grid.
  std::complex<T>&         Complex(size_t i, size_t j, size_t k)
  { return complex_data_[GetComplexIndex(i, j, k)]; }

  const std::complex<T>&   Complex(size_t i, size_t j, size_t k) const
  { return complex_data_[GetComplexIndex(i, j, k)]; }

  /// Get copy of grid. Only the original grid, no padding.
  Grid<T>                GetRealGrid() const;
  /// Get copy of the complete grid in Fourier domain.
  Grid<std::complex<T> > GetComplexGrid() const;

  size_t GetRealNI()    const { return ni_; }
  size_t GetRealNJ()    const { return nj_; }
  size_t GetRealNK()    const { return nk_; }
  size_t GetComplexNI() const { return ni_tot_/2 + 1; }
  size_t GetComplexNJ() const { return nj_tot_; }
  size_t GetComplexNK() const { return nk_tot_; }
  size_t GetNItot()     const { return ni_tot_; }
  size_t GetNJtot()     const { return nj_tot_; }
  size_t GetNKtot()     const { return nk_tot_; }

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
  size_t nk_;
  size_t ni_tot_;
  size_t nj_tot_;
  size_t nk_tot_;

  /// Real data. Column-major ordering.
  T*               real_data_;
  /// Complex data. Column major ordering.
  std::complex<T>* complex_data_;

  // Make copying illegal.
  FFTGrid3D(FFTGrid3D<T>& rhs);
  FFTGrid3D<T>& operator=(FFTGrid3D<T>& rhs);


  /// Find array index in the main real grid.
  inline size_t GetRealIndex(size_t i, size_t j, size_t k) const;

  /// Find array index in the main real grid.
  inline size_t GetRealIndexCyclic(int i, int j, int k) const;

  /// Find array index in the main real grid.
  inline size_t GetComplexIndex(size_t i, size_t j, size_t k) const;
};

namespace NRLibPrivate {
  /// Computes 3D FFT transform. No scaling. Assumes column major storage.
  template <typename T>
  void ComputeFFT3D(size_t ni, size_t nj, size_t nk, T* in, std::complex<T>* out);

  /// Computes inverse 3D FFT transform. No scaling. Assumes column major storage.
  template <typename T>
  void ComputeFFT3DInverse(size_t ni, size_t nj, size_t nk, std::complex<T>* in, T* out);

  // Declare specialization - to ensure Linux compilation.
  template <>
  void ComputeFFT3D<float>(size_t ni, size_t nj, size_t nk, float* in, std::complex<float>* out);
  template <>
  void ComputeFFT3D<double>(size_t ni, size_t nj, size_t nk, double* in, std::complex<double>* out);
  template <>
  void ComputeFFT3DInverse<float>(size_t ni, size_t nj, size_t nk,
                                  std::complex<float>* in, float* out);
  template <>
  void ComputeFFT3DInverse<double>(size_t ni, size_t nj, size_t nk,
                                   std::complex<double>* in, double* out);

}


// -----------------  TEMPLATE IMPLEMENTATION ---------------------

template <typename T>
FFTGrid3D<T>::FFTGrid3D(size_t ni, size_t nj,  size_t nk, size_t padding_ni,
                        size_t padding_nj, size_t padding_nk, bool scale_forward)
  : scale_forward_(scale_forward), ni_(ni), nj_(nj), nk_(nk)
{
  // Find total sizes.
  ni_tot_ = FindNewSizeWithPadding(ni + padding_ni, true);
  nj_tot_ = FindNewSizeWithPadding(nj + padding_nj);
  nk_tot_ = FindNewSizeWithPadding(nk + padding_nk);

  // Allocate aligned data for efficiency.
  // Add padding for use in inplace transform.
  real_data_ = reinterpret_cast<T*>(fftw_malloc(ni_tot_ * nj_tot_ * nk_tot_ * sizeof(T)));

  size_t n_complex = GetComplexNI() * GetComplexNJ() * GetComplexNK();
  size_t complex_datalen = n_complex * sizeof(std::complex<T>);
  complex_data_ = reinterpret_cast<std::complex<T>*>(fftw_malloc(complex_datalen));

  for (size_t i = 0; i < ni_tot_*nj_tot_*nk_tot_; ++i)
    real_data_[i] = 0;

  for (size_t i = 0; i < n_complex; ++i)
    complex_data_[i] = 0;
}


template <typename T>
FFTGrid3D<T>::~FFTGrid3D()
{
  fftw_free(real_data_);
  fftw_free(complex_data_);
}


template <typename T>
void FFTGrid3D<T>::DoFFT() {
  if (scale_forward_) {
    size_t n = ni_tot_ * nj_tot_ * nk_tot_;
    double scale = 1.0 / sqrt(1.0*n);
    for (size_t i = 0; i < n; ++i)
      real_data_[i] *= static_cast<T>(scale);
  }
  NRLibPrivate::ComputeFFT3D(ni_tot_, nj_tot_, nk_tot_, real_data_, complex_data_);
}

template <typename T>
void FFTGrid3D<T>::DoInverseFFT()
{
  NRLibPrivate::ComputeFFT3DInverse(ni_tot_, nj_tot_, nk_tot_, complex_data_, real_data_);
  double scale;
  size_t n = ni_tot_ * nj_tot_ * nk_tot_;
  if (!scale_forward_)
    scale=1.0 / n;
  else
    scale=1.0 / sqrt(1.0*n);

  for (size_t i = 0; i < n; ++i)
    real_data_[i] *= static_cast<T>(scale);
}


template <typename T>
size_t FFTGrid3D<T>::GetRealIndex(size_t i, size_t j, size_t k) const
{
  assert (i < GetRealNI() && j < GetRealNJ() && k < GetRealNK());
  return i + j*ni_tot_ + k*ni_tot_*nj_tot_;
}


template <typename T>
size_t FFTGrid3D<T>::GetRealIndexCyclic(int i, int j, int k) const
{
  int ni = static_cast<int>(GetNItot());
  int nj = static_cast<int>(GetNJtot());
  int nk = static_cast<int>(GetNKtot());

#ifndef NDEBUG
  assert (i <  ni && j <  nj && k <  nk);
  assert (i > -ni && j > -nj && k > -nk);
#endif

  if (i < 0)
    i = ni + i;
  if (j < 0)
    j = nj + j;
  if (k < 0)
    k = nk + k;

  return static_cast<size_t>(i + j*ni + k*ni*nj);
}


template <typename T>
size_t FFTGrid3D<T>::GetComplexIndex(size_t i, size_t j, size_t k) const
{
  assert (i < GetComplexNI() && j < GetComplexNJ() && k < GetComplexNK());
  return i + j*GetComplexNI() + k*GetComplexNI()*GetComplexNJ();
}


template <typename T>
void FFTGrid3D<T>::Initialize(const Grid<T> &values)
{
  assert(ni_tot_ == values.GetNI() && nj_tot_ == values.GetNJ() && nk_tot_ == values.GetNK());
  size_t index(0);

  for (size_t k = 0; k < values.GetNK(); k++)
    for (size_t j = 0; j < values.GetNJ(); j++)
      for (size_t i = 0; i < values.GetNI(); i++) {
        real_data_[index] = values(i,j,k);
        index++;
      }
}


template <typename T>
void FFTGrid3D<T>::InitializePadZero(const Grid<T> &values)
{
  assert(ni_ == values.GetNI() && nj_ == values.GetNJ() && nk_ == values.GetNK());
  size_t index(0);

  for (size_t k = 0; k < nk_tot_; k++)
    for (size_t j = 0; j < nj_tot_; j++)
      for (size_t i = 0; i < ni_tot_; i++) {
        if (i < ni_ && j < nj_ && k < nk_)
          real_data_[index] = values(i,j,k);
        else
          real_data_[index] = static_cast<T>(0.0);
        index++;
      }
}


template <typename T>
void FFTGrid3D<T>::InitializeConstant(const T& value)
{
  for (size_t ind = 0; ind < ni_tot_*nj_tot_*nk_tot_; ++ind)
    real_data_[ind] = value;
}


template <typename T>
void FFTGrid3D<T>::Convolve(FFTGrid3D<T> &filter)
{
  DoFFT();
  filter.DoFFT();
  for(size_t i = 0; i < GetComplexNI()*GetComplexNJ()*GetComplexNK(); i++)
  {
    complex_data_[i] *= filter.complex_data_[i];
  }
  DoInverseFFT();
}


template <typename T>
void FFTGrid3D<T>::ConvolveCovariance(FFTGrid3D<T> &filter)
{
  DoFFT();
  filter.DoFFT();
  for(size_t i = 0; i < GetComplexNI()*GetComplexNJ()*GetComplexNK(); i++)
  {
    complex_data_[i] *= std::sqrt(filter.complex_data_[i]);
  }
  DoInverseFFT();
}


template <typename T>
void FFTGrid3D<T>::ConvolveNoFFT(const FFTGrid3D<T> &filter)
{
  for(size_t i = 0; i < GetComplexNI()*GetComplexNJ()*GetComplexNK(); i++)
  {
    complex_data_[i] *= filter.complex_data_[i];
  }
}


template <typename T>
Grid<T>  FFTGrid3D<T>::GetRealGrid() const
{
  size_t i,j,k;
  Grid<T> output(ni_,nj_,nk_);
  for(i=0;i<ni_;i++)
    for(j=0;j<nj_;j++)
      for(k=0;k<nk_;k++)
        output(i,j,k) = real_data_[i+ni_tot_*j+ni_tot_*nj_tot_*k];

  return output;

}


#ifdef FFTW_DEBUG
template <typename T>
void FFTGrid3D<T>::WriteRealToFile(const std::string& filename)
{
  std::ofstream out;
  NRLib::OpenWrite(out, filename);

  for (size_t k = 0; k < nk_tot_; ++k) {
    for (size_t j = 0; j < nj_tot_; ++j) {
      for (size_t i = 0; i < ni_tot_; ++i) {
        out << real_data_[i+ni_tot_*j+ni_tot_*nj_tot_*k]<< " ";
      }
    }
    out << "\n";
  }
}


template <typename T>
void FFTGrid3D<T>::WriteComplexToFile(const std::string& filename)
{
  ofstream outreal;
  NRLib::OpenWrite(outreal, "real_"+filename);

  for (size_t k = 0; k < nk_tot_; ++k) {
    for (size_t j = 0; j < nj_tot_; ++j) {
      for (size_t i = 0; i < ni_tot_/2+1; ++i) {
        outreal << complex_data_[i+ni_tot_*j+ni_tot_*nj_tot_*k].real() << " ";
      }
    }
    outreal << "\n";
  }

  ofstream outimag;
  NRLib::OpenWrite(outimag, "imag_"+filename);

  for (size_t k = 0; k < GetComplexNK(); ++k) {
    for (size_t j = 0; j < GetComplexNJ(); ++j) {
      for (size_t i = 0; i < GetComplexNI(); ++i) {
        outimag << complex_data_[i+ni_tot_*j+ni_tot_*nj_tot_*k].imag()<< " ";
      }
    }
   outimag << "\n";
  }
}
#endif


// ----------------  FFTW interface functions ---------------------

// Non-specialized versions, should give error.

template <typename T>
void NRLibPrivate::ComputeFFT3D(size_t ni, size_t nj, size_t nk,
                                T* in, std::complex<T>* out)
{
  // 3D FFT computations are only implemented for float and double.
  assert(0);
}

template <typename T>
void NRLibPrivate::ComputeFFT3DInverse(size_t ni, size_t nj, size_t nk,
                                       std::complex<T>* in, T* out)
{
  // 3D FFT computations are only implemented for float and double.
  assert(0);
}


} // namespace NRLib


#endif // NRLIB_FFT_FFT3D_HPP
