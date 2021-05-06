// $Id: gaussianfield.cpp 1694 2017-09-21 13:03:38Z hauge $

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

#include <cmath>
#include <algorithm>
#include "gaussianfield.hpp"

#include "variogram.hpp"
#include "../grid/grid2d.hpp"
#include "../grid/grid.hpp"
#include "../random/random.hpp"

#include "../fft/fftgrid2d.hpp"
#include "../fft/fftgrid3d.hpp"
#include "../fft/fft.hpp"



using namespace NRLib;



// NEW
// Local function, not accessible elsewhere.
// Using int instead of size_t, since we want negative differences.
Grid<double> CovGridForFFT3D(const Variogram & variogram,
               int               nx,
               double            dx,
               int               ny,
               double            dy,
               int               nz,
               double            dz)
{
  Grid<double> cov(nx, ny, nz);
  int nxm = (nx + 1) / 2;
  int nym = (ny + 1) / 2;
  int nzm = (nz + 1) / 2;
  double ddx, ddy, ddz;
  int i, j, k;
  for (k = 0; k < nzm; k++) {
    ddz = k * dz;
    for (j = 0; j < nym; j++) {
      ddy = j * dy;
      for (i = 0; i < nxm; i++) {
        ddx = i * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }

      for (i = nxm; i < nx; i++) {
        ddx = (i - nx) * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }
    }
  }

  for (k = 0; k < nzm; k++) {
    ddz = k * dz;
    for (j = nym; j < ny; j++) {
      ddy = (j - ny) * dy;
      for (i = nxm; i < nx; i++) {
        ddx = (i - nx)* dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }

      for (i = 0; i < nxm; i++) {
        ddx = i*dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }
    }
  }

  for (k = nzm; k < nz; k++) {
    ddz = (k - nz) * dz;
    for (j = 0; j < nym; j++) {
      ddy = j * dy;
      for (i = 0; i < nxm; i++) {
        ddx = i * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }

      for (i = nxm; i < nx; i++) {
        ddx = (i - nx) * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }
    }
  }


  for (k = nzm; k < nz; k++) {
    ddz = (k - nz) * dz;
    for (j = nym; j < ny; j++) {
      ddy = (j - ny) * dy;
      for (i = 0; i < nxm; i++) {
        ddx = i * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }

      for (i = nxm; i < nx; i++) {
        ddx = (i - nx) * dx;
        cov(i, j, k) = variogram.GetCov(ddx, ddy, ddz);
      }
    }
  }

  return cov;
}




// NEW: 3D version
// Simulation of multiple fields with the same covariance function.
void NRLib::Simulate3DGaussianField(const Variogram              & variogram,
                  size_t                         nx,
                  double                         dx,
                  size_t                         ny,
                  double                         dy,
                  size_t                         nz,
                  double                         dz,
                  int                            n_fields,
                  std::vector<Grid<double> > &   grid_out)
{
  // Find grid size.
  double range_x, range_y, range_z;
  range_x = std::max(variogram.GetRangeX(), variogram.GetRangeY());
  range_y = range_x;
  range_z = variogram.GetRangeZ();   // added

  size_t n_pad_x = std::max(static_cast<size_t>(4.0 * range_x / dx), 2 * nx) - nx;
  size_t n_pad_y = std::max(static_cast<size_t>(4.0 * range_y / dy), 2 * ny) - ny;
  size_t n_pad_z = std::max(static_cast<size_t>(4.0 * range_z / dz), 2 * nz) - nz; // added

  //FFTGrid2D<double> fftgrid(nx, ny, n_pad_x, n_pad_y, true);
  FFTGrid3D<double> fftgrid(nx, ny, nz, n_pad_x, n_pad_y, n_pad_z, true); // added

  size_t nx_tot = fftgrid.GetNItot();
  size_t ny_tot = fftgrid.GetNJtot();
  size_t nz_tot = fftgrid.GetNKtot(); // added

  // Covariance grid.
  Grid<double> cov = CovGridForFFT3D(variogram,
                     static_cast<int>(nx_tot), dx,
                     static_cast<int>(ny_tot), dy,
                     static_cast<int>(nz_tot), dz);

  FFTGrid3D<double> filter(nx, ny, nz, n_pad_x, n_pad_y, n_pad_z, false); // added


  filter.Initialize(cov);
  filter.DoFFT();
  // Convolve with covariance grid.
  for (size_t i = 0; i < filter.GetComplexNI()*filter.GetComplexNJ()*filter.GetComplexNK(); i++)
    filter.ComplexData()[i] = std::sqrt(filter.ComplexData()[i]);


  // Grid with white noice.
  Grid<double> noise(nx_tot, ny_tot, nz_tot);

  ////if (seed == -999)
  ////  seed = NRLib::Random::DrawUint32();
  //// This is done to be able to do multithreading
  ////NRLib::RandomGenerator random;
  //// random.Initialize(seed);
  //if (rg == NULL) {
  //  rg = new RandomGenerator();
  //  rg->Initialize(NRLib::Random::DrawUint32());
  //}
  for (int m = 0; m < n_fields; m++) {
    for (size_t i = 0; i < nx_tot; i++)
      for (size_t j = 0; j < ny_tot; j++)
        for (size_t k = 0; k < nz_tot; k++)
          noise(i, j, k) = NRLib::Random::Norm01();

    fftgrid.Initialize(noise);

    fftgrid.DoFFT();
    for (size_t i = 0; i < filter.GetComplexNI()*filter.GetComplexNJ()*filter.GetComplexNK(); i++)
      fftgrid.ComplexData()[i] *= filter.ComplexData()[i];

    fftgrid.DoInverseFFT();
    grid_out.push_back(fftgrid.GetRealGrid());
  }
}




// Local function, not accessible elsewhere.
// Using int instead of size_t, since we want negative differences.
Grid2D<double> CovGridForFFT2D(const Variogram & variogram,
                               int               nx,
                               double            dx,
                               int               ny,
                               double            dy)
{
  Grid2D<double> cov(nx, ny);
  int nxm = (nx + 1)/2;
  int nym = (ny + 1)/2;
  double ddx, ddy;
  int i, j;
  for (j = 0; j < nym; j++) {
    ddy = j * dy;
    for (i = 0; i < nxm; i++) {
      ddx = i * dx;
      cov(i, j) = variogram.GetCov(ddx,ddy);
    }
     // -1 - i quadrant
    for (i = nxm; i < nx; i++) {
      ddx = (i - nx) * dx;
      cov(i,j) = variogram.GetCov(ddx, ddy);
    }
  }

  // 1 + i quadrant
  for (j = nym; j < ny; j++) {
    ddy = (j - ny) * dy;
    for (i = nxm; i < nx; i++) {
      ddx = (i - nx)* dx;
      cov(i,j) = variogram.GetCov(ddx, ddy);
    }
      // -1 + i quadrant
     for (i = 0; i < nxm; i++) {
      ddx = i*dx;
      cov(i,j) = variogram.GetCov(ddx, ddy);
    }
  }


  return cov;
}

void NRLib::Simulate2DGaussianField(const Variogram& variogram,
                                    size_t nx, double dx,
                                    size_t ny, double dy,
                                    // double padding_fraction,
                                    // bool   user_defined_padding,
                                    Grid2D<double> & grid_out)
{
  // Find grid size.
  double range_x, range_y;
  if (variogram.GetAzimuthAngle() == 0.0) {
    range_x = variogram.GetRangeX();
    range_y = variogram.GetRangeY();
  }
  else {
    range_x = std::max(variogram.GetRangeX(), variogram.GetRangeY());
    range_y = range_x;
  }
  size_t n_pad_x = std::max(static_cast<size_t>(4.0 * range_x / dx), 2*nx) - nx;
  size_t n_pad_y = std::max(static_cast<size_t>(4.0 * range_y / dy), 2*ny) - ny;

  FFTGrid2D<double> fftgrid(nx, ny, n_pad_x, n_pad_y, true);

  size_t nx_tot = fftgrid.GetNItot();
  size_t ny_tot = fftgrid.GetNJtot();

  // Covariance grid.
  Grid2D<double> cov = CovGridForFFT2D(variogram,
                                       static_cast<int>(nx_tot), dx,
                                       static_cast<int>(ny_tot), dy);

  // Grid with white noice.
  Grid2D<double> noise(nx_tot,ny_tot);
  size_t i,j;
    for (i = 0; i < nx_tot; i++)
      for (j = 0; j < ny_tot; j++)
        noise(i,j) = NRLib::Random::Norm01();

  fftgrid.Initialize(noise);
  FFTGrid2D<double> filter(nx, ny, n_pad_x, n_pad_y, false);
  filter.Initialize(cov);
  fftgrid.ConvolveCovariance(filter);
  grid_out = fftgrid.GetRealGrid();
}

// Simulation of multiple fields with the same covariance function.
void NRLib::Simulate2DGaussianField(const Variogram              & variogram,
                                    size_t                         nx,
                                    double                         dx,
                                    size_t                         ny,
                                    double                         dy,
                                    int                            n_fields,
                                    std::vector<Grid2D<double> > & grid_out,
                                    NRLib::RandomGenerator        *rg)
{
  // Find grid size.
  double range_x, range_y;
  if (variogram.GetAzimuthAngle() == 0.0) {
    range_x = variogram.GetRangeX();
    range_y = variogram.GetRangeY();
  }
  else {
    range_x = std::max(variogram.GetRangeX(), variogram.GetRangeY());
    range_y = range_x;
  }
  size_t n_pad_x = std::max(static_cast<size_t>(4.0 * range_x / dx), 2*nx) - nx;
  size_t n_pad_y = std::max(static_cast<size_t>(4.0 * range_y / dy), 2*ny) - ny;

  FFTGrid2D<double> fftgrid(nx, ny, n_pad_x, n_pad_y, true);

  size_t nx_tot = fftgrid.GetNItot();
  size_t ny_tot = fftgrid.GetNJtot();

  // Covariance grid.
  Grid2D<double> cov = CovGridForFFT2D(variogram,
                                       static_cast<int>(nx_tot), dx,
                                       static_cast<int>(ny_tot), dy);


  FFTGrid2D<double> filter(nx, ny, n_pad_x, n_pad_y, false);
  filter.Initialize(cov);
  filter.DoFFT();
  // Convolve with covariance grid.
  for (size_t i = 0; i < filter.GetComplexNI()*filter.GetComplexNJ(); i++)
    filter.ComplexData()[i] = std::sqrt(filter.ComplexData()[i]);


  // Grid with white noice.
  Grid2D<double> noise(nx_tot,ny_tot);

  //if (seed == -999)
  //  seed = NRLib::Random::DrawUint32();
  // This is done to be able to do multithreading
  //NRLib::RandomGenerator random;
 // random.Initialize(seed);
  if(rg == NULL) {
    rg = new RandomGenerator();
    rg->Initialize(NRLib::Random::DrawUint32());
  }
  for (int k = 0; k < n_fields; k++) {
    for (size_t i = 0; i < nx_tot; i++)
      for (size_t j = 0; j < ny_tot; j++)
        noise(i,j) = rg->Norm01();

    fftgrid.Initialize(noise);

    fftgrid.DoFFT();
    for(size_t i = 0; i < filter.GetComplexNI()*filter.GetComplexNJ(); i++)
      fftgrid.ComplexData()[i] *=filter.ComplexData()[i];

    fftgrid.DoInverseFFT();
    grid_out.push_back(fftgrid.GetRealGrid());
  }
}

void
NRLib::Simulate1DGaussianField(const Variogram       & variogram,
                               size_t                  nx,
                               double                  dx,
                               std::vector<double>   & grid_out,
                               NRLib::RandomGenerator *rg)
{
  double range = variogram.GetRangeX();
  size_t nx_pad = std::min(static_cast<size_t>(2.0 * range / dx), nx);
  size_t nxp = NRLib::FindNewSizeWithPadding(nx + nx_pad);
  nx_pad =nxp-nx;

  std::vector<double> cov(nxp);
  size_t nxm = (nxp+1)/2;
  for(size_t i = 0; i < nxm; i++) {
    cov[i] = variogram.GetCov(i*dx);
  }
  for(size_t i=nxm; i < nxp; i++) {
    double ddx = dx*(nxp - i);
    cov[i] = variogram.GetCov(ddx);
  }
  if(rg == NULL) {
    rg = new RandomGenerator();
    rg->Initialize(NRLib::Random::DrawUint32());
  }
  std::vector<double> noise(nxp);
  for(size_t i = 0; i < nxp; i++)
    noise[i] = rg->Norm01();
  std::vector<std::complex<double> > cov_fft;
  std::vector<std::complex<double> > noise_fft;
  NRLib::ComputeFFT1D(cov, cov_fft, false,0);
  NRLib::ComputeFFT1D(noise, noise_fft, true,0);
  std::vector<std::complex<double> > convolve(cov_fft.size());

  for(size_t i = 0; i < cov_fft.size(); i++)
    convolve[i] = std::sqrt(cov_fft[i])*noise_fft[i];

  std::vector<double> tmp(nxp);

  NRLib::ComputeFFTInv1D(convolve, tmp, true);

  for(size_t i = 0; i < grid_out.size(); i++)
    grid_out[i] = tmp[i];
}


/*
void
NRLib::Simulate1DGaussianField(NRLib::Matrix cov_in,
                               size_t nx,
                               std::vector<double> & grid_out)
{
 // double range = variogram.GetRangeX();
 // size_t nx_pad = std::max(static_cast<size_t>(4.0 * range / dx), 2*nx) - nx;
  size_t nx_pad = 2*nx;
  std::vector<double> cov(nx_pad);
 // size_t nxm = (nx_pad+1)/2;
  for(size_t i = 0; i < nx; i++){
    cov[i] = cov_in(static_cast<int> (i),0);
  }
  for(size_t i=nx; i < nx_pad; i++){
    cov[i] = cov_in(static_cast<int> (nx_pad-i),0);
  }

  std::vector<double> noise(2*nx_pad);
  for(size_t i = 0; i < 2*nx_pad; i++)
    noise[i] = NRLib::Random::Norm01();
  std::vector<std::complex<double> > cov_fft;
  std::vector<std::complex<double> > noise_fft;
  NRLib::ComputeFFT1D(cov, cov_fft, true);
  NRLib::ComputeFFT1D(noise, noise_fft, true);
  std::vector<std::complex<double> > convolve(cov_fft.size());

  for(size_t i = 0; i < cov_fft.size(); i++)
    convolve[i] = std::sqrt(cov_fft[i])*noise_fft[i];


  NRLib::ComputeFFTInv1D(convolve, grid_out, true);
}
*/
