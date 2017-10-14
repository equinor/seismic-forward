// $Id: fft.cpp 1168 2013-05-07 06:07:04Z anner $

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

#include "fft.hpp"

template <>
void NRLib::NRLibPrivate::ComputeFFT1D<double>(size_t n,
                                               double* in,
                                               std::complex<double>* out)
{
  fftw_complex* out_data = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_r2c_1d(static_cast<int>(n), in, out_data, FFTW_ESTIMATE);
  assert (p != 0);

  fftw_execute(p);
  fftw_destroy_plan(p);
}

template <>
void NRLib::NRLibPrivate::ComputeFFTInv1D<double>(size_t n,
                                                  std::complex<double>* in,
                                                  double* out)
{
  fftw_complex* in_data = reinterpret_cast<fftw_complex*>(in);
  fftw_plan p = fftw_plan_dft_c2r_1d(static_cast<int>(n), in_data, out, FFTW_ESTIMATE);
  assert (p != 0);

  fftw_execute(p);
  fftw_destroy_plan(p);
}

size_t
NRLib::FindNewSizeWithPadding(size_t minSize, bool must_be_even) {
  int i,j,k,l,m,n;
  size_t factor   =       1;

  int minant2 = 0;
  if (must_be_even)
    minant2 = 1;

  int maxant2    = static_cast<int>(ceil((double) log((float)(minSize)) / log(2.0f) ));
  int maxant3    = static_cast<int>(ceil((double) log((float)(minSize)) / log(3.0f) ));
  int maxant5    = static_cast<int>(ceil((double) log((float)(minSize)) / log(5.0f) ));
  int maxant7    = static_cast<int>(ceil((double) log((float)(minSize)) / log(7.0f) ));
  int maxant11   = 0;
  int maxant13   = 0;
  size_t closestprod= static_cast<size_t>(pow(2.0f,maxant2));

  // kan forbedres ved aa trekke fra i endepunktene.i for lokkene
  for (i = minant2; i < maxant2+1; i++)
    for (j = 0; j < maxant3+1; j++)
      for (k = 0; k < maxant5+1; k++)
        for (l = 0; l < maxant7+1; l++)
          for (m = 0; m < maxant11+1; m++)
            for (n = maxant11; n < maxant13+1; n++) {
              factor = static_cast<size_t> (pow(2.0f,i)*pow(3.0f,j)*pow(5.0f,k)*
                pow(7.0f,l)*pow(11.0f,m)*pow(13.0f,n));
              if ((factor >=  minSize) &&  (factor <  closestprod)) {
                closestprod = factor;
              }
            }

  return closestprod;
}
