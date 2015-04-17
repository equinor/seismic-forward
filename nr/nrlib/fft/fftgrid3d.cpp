// $Id: fftgrid3d.cpp 994 2012-02-29 15:15:04Z perroe $

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

#include "fftgrid3d.hpp"

namespace NRLib {

template <>
void NRLibPrivate::ComputeFFT3D<double>(size_t ni, size_t nj, size_t nk,
                                        double* in, std::complex<double>* out)
{
  fftw_complex* out_data = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_r2c_3d(static_cast<int>(nk), static_cast<int>(nj), static_cast<int>(ni),
                                     in, out_data, FFTW_ESTIMATE);
  assert (p != 0);

  fftw_execute(p);
  fftw_destroy_plan(p);
}


template <>
void NRLibPrivate::ComputeFFT3DInverse<double>(size_t ni, size_t nj, size_t nk,
                                               std::complex<double>* in, double* out)
{
  fftw_complex* in_data = reinterpret_cast<fftw_complex*>(in);
  fftw_plan p = fftw_plan_dft_c2r_3d(static_cast<int>(nk), static_cast<int>(nj), static_cast<int>(ni),
                                     in_data, out, FFTW_ESTIMATE);
  assert (p != 0);

  fftw_execute(p);
  fftw_destroy_plan(p);
}


template <>
void NRLibPrivate::ComputeFFT3D<float>(size_t ni, size_t nj, size_t nk,
                                       float* in, std::complex<float>* out)
{
  fftwf_complex* out_data = reinterpret_cast<fftwf_complex*>(out);
  fftwf_plan p = fftwf_plan_dft_r2c_3d(static_cast<int>(nk), static_cast<int>(nj), static_cast<int>(ni),
                                       in, out_data, FFTW_ESTIMATE);
  assert (p != 0);

  fftwf_execute(p);
  fftwf_destroy_plan(p);
}


template <>
void NRLibPrivate::ComputeFFT3DInverse<float>(size_t ni, size_t nj, size_t nk,
                                              std::complex<float>* in, float* out)
{
  fftwf_complex* in_data = reinterpret_cast<fftwf_complex*>(in);
  fftwf_plan p = fftwf_plan_dft_c2r_3d(static_cast<int>(nk), static_cast<int>(nj), static_cast<int>(ni),
                                       in_data, out, FFTW_ESTIMATE);
  assert (p != 0);

  fftwf_execute(p);
  fftwf_destroy_plan(p);
}

} // namespace NRLib
