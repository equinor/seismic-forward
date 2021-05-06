// $Id: kriging.cpp 1328 2016-04-26 13:40:28Z perroe $

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

#include <vector>
#include <cstring>
#include <typeinfo>


#include "kriging.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/variogram/covgrid2d.hpp"
#include "nrlib/variogram/variogramtypes.hpp"
#include "nrlib/statistics/nodalkrigingdata2d.hpp"
#include "nrlib/statistics/continuouskrigingdata2d.hpp"
#include "nrlib/statistics/krigingdata2d.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#ifdef WITH_OMP
#include <omp.h>
#endif

////for debug print
//#include <nrlib/surface/regularsurface.hpp>
//#include <nrlib/surface/surfaceio.hpp>

#ifdef PARALLEL
#include <omp.h>
#endif

using namespace NRLib;

void Kriging::Krig1D(std::vector<double>       &field,
                     const std::vector<bool>   &is_known,
                     const std::vector<double> &obs,
                     double                     dx,
                     const Variogram           &vario)
{
  assert(field.size() == is_known.size());
  assert(obs.size() <= field.size());
  int n_obs = static_cast<int>(obs.size());

  NRLib::Vector residual(static_cast<int> (n_obs));
  int j = 0;
  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true){
      residual(j) = (obs[j]-field[i]);
      j++;
    }

  }

  NRLib::SymmetricMatrix K(n_obs);
  NRLib::Vector k_vec(n_obs);
  std::vector<double> x_known;
  std::vector<double> x_unknown;
  NRLib::Vector K_res(n_obs);

  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true)
      x_known.push_back(i*dx);
    else
      x_unknown.push_back(i*dx);
  }


  for(int i = 0; i < n_obs; i++)
    for(j = 0; j <= i; j++)
      K(j,i) = vario.GetCorr(x_known[i]-x_known[j]);

  CholeskyInvert(K);
  K_res = K * residual;
  j = 0;
  size_t kk = 0;
  for(size_t k = 0; k < field.size(); k++){
    if(is_known[k] == false){
      for(int i = 0; i < n_obs; i++)
        k_vec(i) = vario.GetCorr(x_known[i]-x_unknown[j]);

      double prod = k_vec * K_res;
      field[k] += prod;
      j++;
    }
    else{
      field[k] = obs[kk];
      kk++;
    }
  }
}


//----------------------------------------------------------
void Kriging::Krig2D(Grid2D<double>           & trend,
                     const NodalKrigingData2D & kriging_data,
                     const CovGrid2D          & cov,
                     bool                       get_residuals,
                     bool                       krig_all,
                     double                     P)
//----------------------------------------------------------
{
  // Unused parameters needed in common function
  ConstVario             dummy_variogram(0, 0, 0, 0, 0, 0);
  RegularSurface<double> dummy_trend(0, 0, 0, 0, 0, 0, 0);

  Krig2D(trend,
         dummy_trend,
         (*dynamic_cast<const KrigingData2D *>(&kriging_data)),
         cov,
         dummy_variogram,
         get_residuals,
         cov.GetDX(),
         cov.GetDY(),
         cov.GetRangeX(),
         cov.GetRangeY(),
         P,
         krig_all);

}

//------------------------------------------
bool compare(const std::pair<int, int> & i1,
             const std::pair<int, int> & i2)
//------------------------------------------
{
  return (i1.second > i2.second);
}


//----------------------------------------------------------
void Kriging::Krig2D(RegularSurface<double>      & trend,
                     const std::vector<Point* >  & data_points,
                     const Variogram             & var,
                     bool                          get_residuals,
                     bool                          use_index_grid,
                     bool                          krig_all,
                     double                        P)
//----------------------------------------------------------
{
  // When using index based grid, then we convert the vector of data points
  // to KrigingData2D and the Variogram to CovGrid2D.
  if(use_index_grid) {
    const NodalKrigingData2D kriging_data(data_points, trend);
    const CovGrid2D cov(var,
                        static_cast<int>(trend.GetNI()),
                        static_cast<int>(trend.GetNJ()),
                        trend.GetDX(),
                        trend.GetDY());

    Krig2D(trend, kriging_data, cov, get_residuals, krig_all, P);

  }
  else {
    ContinuousKrigingData2D kriging_data_continuous(data_points, trend);
    KrigingData2D * kriging_data = &kriging_data_continuous;

    Grid2D<double>  * trend_grid2d = dynamic_cast<Grid2D<double> *>(&trend);
    const CovGrid2D dummy_cov(var, 0, 0, 0, 0);

    Krig2D(*trend_grid2d,
           trend,
           *kriging_data,
           dummy_cov,
           var,
           get_residuals,
           trend.GetDX(),
           trend.GetDY(),
           var.GetRangeX(),
           var.GetRangeY(),
           P,
           krig_all);

  }
}

// Common function for both index based kriging and index kriging and continuous data and continuous kriging.
void Kriging::Krig2D(Grid2D<double>          & trend,
                     RegularSurface<double>  & trend_surface,
                     const KrigingData2D     & kriging_data,
                     const CovGrid2D         & cov,
                     const Variogram         & var,
                     bool                      get_residuals,
                     double                    dx,
                     double                    dy,
                     double                    Rx,
                     double                    Ry,
                     double                    P,
                     bool                      krig_all)
{
  // blocking
  size_t nxb, nyb, n_blocks_x, n_blocks_y;
  std::vector<size_t> block_x, block_y;
  std::vector<KrigingData2D *> kriging_data_blocks;
  size_t n_data        = kriging_data.GetNumberOfData();
  size_t n_data_blocks = n_data;

  int n_threads            = 1;

#ifdef PARALLEL
  n_threads = omp_get_num_procs();
#endif

  if (krig_all == false){
    FindOptimalBlockP(trend,
                     Rx,
                     Ry,
                     dx,
                     dy,
                     n_data,
                     n_threads,
                     n_blocks_x,
                     n_blocks_y,
                     P,
                     n_data_blocks); // returns n_block_x and n_block_y and n_data_blocks
  }
  else {
    n_blocks_x = 1;
    n_blocks_y = 1;
  }
  size_t nx = trend.GetNI();
  size_t ny = trend.GetNJ();

  // number of grid cells in block in x and y direction
  nxb = static_cast<size_t>(floor(nx/static_cast<double>(n_blocks_x)+0.5));
  nyb = static_cast<size_t>(floor(ny/static_cast<double>(n_blocks_y)+0.5));

  size_t n_blocks = n_blocks_x * n_blocks_y;
  kriging_data_blocks.resize(n_blocks);
  block_x.resize(n_blocks);
  block_y.resize(n_blocks);

  if (krig_all == false){
    FindOptimalDataInBlock(trend,
                           kriging_data,
                           Rx,
                           Ry,
                           dx,
                           dy,
                           n_blocks_x,
                           n_blocks_y,
                           P,
                           n_data_blocks,
                           kriging_data_blocks,
                           block_x,
                           block_y);             // return kriging_data_blocks, block_x and block_y
  }
  else {
    kriging_data_blocks[0] = kriging_data.Clone();  // Pointers deleted in Krig2D()
    block_x[0] = 0;
    block_y[0] = 0;
  }

  // kriging in parallel
  Grid2D<double> filled(nx,ny,0);
  const Grid2D<double> trend_orig(trend);
  const RegularSurface<double> trend_orig_surface(trend_surface);

  bool nodal_data      = false;
  bool continuous_data = false;

  // For portable comparisons of class types
  NodalKrigingData2D      * nodal_class      = new NodalKrigingData2D();
  ContinuousKrigingData2D * continuous_class = new ContinuousKrigingData2D();

   if(std::strcmp(typeid(kriging_data).name(), typeid(*nodal_class).name()) == 0) {
     nodal_data = true;
   }
   else if(std::strcmp(typeid(kriging_data).name(), typeid(*continuous_class).name()) == 0) {
     continuous_data = true;
   }

   delete nodal_class;
   delete continuous_class;

   int  chunk_size;
   chunk_size = 1;
#ifdef PARALLEL
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
   for (size_t i = 0; i < kriging_data_blocks.size(); ++i)
   {
     if(nodal_data) {
       Krig2DBlock(trend_orig,
                   *dynamic_cast<const NodalKrigingData2D *>(kriging_data_blocks[i]),
                   cov,
                   block_x[i],
                   block_y[i],
                   nxb,
                   nyb,
                   n_blocks_x,
                   n_blocks_y,
                   trend,
                   filled,
                   get_residuals);
     }
     else if (continuous_data) {
       Krig2DBlock(trend_orig_surface,
                   *dynamic_cast<const ContinuousKrigingData2D *>(kriging_data_blocks[i]),
                   var,
                   block_x[i],
                   block_y[i],
                   nxb,
                   nyb,
                   n_blocks_x,
                   n_blocks_y,
                   trend_surface,
                   get_residuals);
     }
   }

   // Release memory
   for(int i = 0; i < static_cast<int>(kriging_data_blocks.size()); i++){
     delete kriging_data_blocks[i];
   }
}

//----------------------------------------------------------------------
void Kriging::FindOptimalBlockP(const Grid2D<double> & trend,
                               double                 Rx,
                               double                 Ry,
                               double                 dx,
                               double                 dy,
                               size_t                 n_data,
                               size_t                 n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y,
                               double               & P,
                               size_t               & n_data_blocks)
//----------------------------------------------------------------------
{
  // optimize time with respect to PETROSIM document

  const size_t nx      = trend.GetNI();
  const size_t ny      = trend.GetNJ();
  const double xlength = nx*dx;
  const double ylength = ny*dy;

  const double Vcell   = dx*dy;
  const double V       = xlength * ylength;
  const double R       = Rx * Ry;
  double       mintime = std::numeric_limits<double>::infinity();


  size_t max_points = 200;
  n_blocks_x = 10;
  n_blocks_y = 10;


  if (n_data > max_points){

    size_t min_grid_cells = 5;                                 // minimum size of blocks
    size_t min_blocks     = 5;                                 // minimum number of blocks in each direction
    double Smin           = std::max((min_grid_cells * dx / Rx),(min_grid_cells * dy / Ry));
    double Smax           = std::min((xlength / (min_blocks * Rx)),(ylength / (min_blocks * Ry)));
    double minatS         = Smax;
    size_t NS             = 200;
    double factor_x       = 1.0;
    double factor_y       = 1.0;

    for (size_t k = 0; k <= NS; ++k){                          // interval for S = [Smin, Smax] is divided into NS equal parts
      double       S       = Smin + k * (Smax - Smin) / NS;
      double       v       = std::max(R*S*S, Vcell);           // volume of segment
      double       vd      = R*(2*P + S)*(2*P + S);                // volume of neighbourhood
      size_t       Ns      = static_cast<size_t>(V/v);
      size_t       nd      = static_cast<size_t>(n_data * (vd / V));//      size_t       nd      = static_cast<size_t>(n_data * std::min((vd / V), 1.0));
      double       totTime = OverallTime(nd, Ns, nx, ny);
      n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(S*Rx)));
      n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(S*Ry)));
      if (totTime <= mintime){
        if ((n_blocks_x * n_blocks_y) >= n_threads){
          minatS  = S;
          mintime = totTime;
          n_data_blocks = nd;
        }
      }
    }
    n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(minatS*Rx)));
    n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(minatS*Ry)));
  }
}

//----------------------------------------------------------------------
void Kriging::FindOptimalBlock(const Grid2D<double> & trend,
                               double                 Rx,
                               double                 Ry,
                               double                 dx,
                               double                 dy,
                               size_t                 n_data,
                               size_t                 n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y)
//----------------------------------------------------------------------
{
  // optimize time with respect to PETROSIM document

  const size_t nx      = trend.GetNI();
  const size_t ny      = trend.GetNJ();
  const double xlength = nx*dx;
  const double ylength = ny*dy;

  const double Vcell   = dx*dy;
  const double V       = xlength * ylength;
  const double R       = Rx * Ry;
  double       mintime = std::numeric_limits<double>::infinity();

  size_t max_points = 200;
  n_blocks_x = 10;
  n_blocks_y = 10;

  if (n_data > max_points){

    size_t min_grid_cells = 5;                                 // minimum size of blocks
    size_t min_blocks     = 5;                                 // minimum number of blocks in each direction
    double Smin           = std::max((min_grid_cells * dx / Rx),(min_grid_cells * dy / Ry));
    double Smax           = std::min((xlength / (min_blocks * Rx)),(ylength / (min_blocks * Ry)));
    double minatS         = Smax;
    size_t NS             = 200;
    double factor_x       = 1.0;
    double factor_y       = 1.0;

    for (size_t k = 0; k <= NS; ++k){                          // interval for S = [Smin, Smax] is divided into NS equal parts
      double       S       = Smin + k * (Smax - Smin) / NS;
      double       v       = std::max(R*S*S, Vcell);           // volume of segment
      double       vd      = R*(2 + S)*(2 + S);                // volume of neighbourhood
      size_t       Ns      = static_cast<size_t>(V/v);
      size_t       nd      = static_cast<size_t>(n_data * vd / V);
      double       totTime = OverallTime(nd, Ns, nx, ny);
      if (totTime <= mintime){
        n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(S*Rx)));
        n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(S*Ry)));
        if ((n_blocks_x * n_blocks_y) >= n_threads){
          minatS  = S;
          mintime = totTime;
        }
      }
    }
    n_blocks_x = static_cast<unsigned int>(std::max(1.0, factor_x*xlength/(minatS*Rx)));
    n_blocks_y = static_cast<unsigned int>(std::max(1.0, factor_y*ylength/(minatS*Ry)));
  }
}


//------------------------------------
double Kriging::OverallTime(size_t n,
                            size_t Ns,
                            size_t nx,
                            size_t ny)
//------------------------------------
{
  // Relative times (i.e. divided by the time it takes to assemble one k vector)
  // Relative times esimtated based on 6 tests
  double       V     = static_cast<double>(nx * ny);
  double       tChol = 0.00032;
  double       tK    = 0.75088;
  double       tSol  = 0.12702;

  double       Tchol = tChol * Ns * n*n*n;
  double       TK    = tK    * Ns * n*n;
  double       Tsol  = tSol  * Ns * n*n;
  double       Tk    = 1.0   * V  * n;

  double       T     = Tchol + TK + Tsol + Tk;

  return T;
}

//----------------------------------------------------
void Kriging::SetConstants(size_t & max_data_in_range,
                           size_t & min_data_in_block,
                           size_t & min_data_in_range,
                           size_t   nxb,
                           size_t   nyb,
                           size_t   rx,
                           size_t   ry,
                           double   P)
//----------------------------------------------------
{
  max_data_in_range = 100000;
  min_data_in_range = std::min(static_cast<int>(min_data_in_range), static_cast<int>(max_data_in_range));

  double l_x = static_cast<double>(nxb)+2*P*static_cast<double>(rx);
  double l_y = static_cast<double>(nyb)+2*P*static_cast<double>(ry);
  double factor = (l_x*l_y - static_cast<double>(nxb)*static_cast<double>(nyb))/(l_x*l_y);

  min_data_in_block = static_cast<size_t>(std::floor(static_cast<double>(min_data_in_range)*factor/8));
}

//-------------------------------------------------------------------------------------------


void Kriging::DataInBlockBinarySearch(double              & P,
                                       double                P_min,
                                       double                P_max,
                                       double                rx,
                                       double                ry,
                                       size_t                nx,
                                       size_t                ny,
                                       size_t                imin,
                                       size_t                imax,
                                       size_t                jmin,
                                       size_t                jmax,
                                       bool                  imin_b,
                                       bool                  imax_b,
                                       bool                  jmin_b,
                                       bool                  jmax_b,
                                       const size_t          n_target,
                                       const KrigingData2D & kriging_data,
                                       const double          tol)
{
  size_t imin_in = imin;
  size_t imax_in = imax;
  size_t jmin_in = jmin;
  size_t jmax_in = jmax;
  if (imin_b == false && imax_b == false && jmin_b == false && jmax_b == false)
    return;
  P = (P_max + P_min) / 2;
  if ((P_max - P_min) < tol)
    return;

  if (imin_b == true){
    imin = std::max(static_cast<int>(std::floor(static_cast<double>(imin_in) - P * rx)), 0);
    imax = imin_in;
  }
  if (imax_b == true) {
    imax = std::min(static_cast<int>(std::floor(static_cast<double>(imax_in) + P * rx)), static_cast<int>(nx));
    if (imin_b == false)
      imin = imax_in;
  }

  if (jmin_b == true) {
    jmin = std::max(static_cast<int>(std::floor(static_cast<double>(jmin_in) - P * ry)), 0);
    jmax = jmin_in;
  }
  if (jmax_b == true) {
    jmax = std::min(static_cast<int>(std::floor(static_cast<double>(jmax_in) + P * ry)), static_cast<int>(ny));
    if (jmin_b == false)
      jmin = jmax_in;
  }
  if (!(imin_b == true && imax_b == true && jmin_b == true && jmax_b == true)){
    if (imin == imin_in) imin_b = false;
    if (imax == imax_in) imax_b = false;
    if (jmin == jmin_in) jmin_b = false;
    if (jmax == jmax_in) jmax_b = false;
  }

  size_t n_data = kriging_data.CountDataInBlock(imin, imax, jmin, jmax);
  if (n_data < n_target){
    P_min = P;
    DataInBlockBinarySearch(P, P_min, P_max, rx, ry, nx, ny, imin_in, imax_in, jmin_in, jmax_in, imin_b, imax_b, jmin_b, jmax_b, n_target, kriging_data, tol);
  }
  else if (n_data > n_target) {
    P_max = P;
    DataInBlockBinarySearch(P, P_min, P_max, rx, ry, nx, ny, imin_in, imax_in, jmin_in, jmax_in, imin_b, imax_b, jmin_b, jmax_b, n_target, kriging_data, tol);
  }
  else
    return;
}




//-------------------------------------------------------------------------------------------
void Kriging::FindOptimalDataInBlock(const Grid2D<double>          & trend,
                                     const KrigingData2D           & kriging_data,
                                     double                          Rx,
                                     double                          Ry,
                                     double                          dx,
                                     double                          dy,
                                     size_t                          n_blocks_x,
                                     size_t                          n_blocks_y,
                                     double                          P,
                                     size_t                          n_avg_data_blocks,
                                     std::vector<KrigingData2D *>  & kriging_data_blocks,
                                     std::vector<size_t>           & block_index_x,
                                     std::vector<size_t>           & block_index_y)
//-------------------------------------------------------------------------------------------
{
  // find number of data in each block and block+range
  size_t nx  = trend.GetNI();
  size_t ny  = trend.GetNJ();
  size_t nxb = static_cast<size_t>(floor(nx/static_cast<double>(n_blocks_x)+0.5));
  size_t nyb = static_cast<size_t>(floor(ny/static_cast<double>(n_blocks_y)+0.5));
  int    rx  = static_cast<int>(Rx/dx);
  int    ry  = static_cast<int>(Ry/dy);

  //A:
  //double P_max = P;
  //double P_min = P;
  //bool binary_search_dir = false;
  //
  //B:
  double P_max = 2*P;
  double P_min = P/2;
  bool binary_search_dir = false;
  //
  //C:
  //double P_max = 2*P;
  //double P_min = P/2;
  //bool binary_search_dir = true;

  double P_est_all = 0;
  double P_min_total = 1000;
  double P_max_total = 0;

  size_t n_blocks = n_blocks_x * n_blocks_y;

  size_t imin,   imax,   jmin,   jmax;

  size_t max_data_in_range = 0;
  size_t min_data_in_block = 0;
  size_t min_data_in_range = n_avg_data_blocks;

  SetConstants(max_data_in_range, min_data_in_block, min_data_in_range, nxb, nyb, rx, ry, P);

  size_t index;
  for (size_t i = 0; i < n_blocks_x; ++i){
    for (size_t j = 0; j < n_blocks_y; ++j){
      index  = i*n_blocks_y + j;
      block_index_x[index] = i;
      block_index_y[index] = j;
    }
  }

  if (n_blocks_x > 1 || n_blocks_y > 1) {
    for (size_t i = 0; i < n_blocks_x; ++i) {
      for (size_t j = 0; j < n_blocks_y; ++j) {

        index = i*n_blocks_y + j;
        kriging_data_blocks[index] = kriging_data.NewInstance();

        size_t imin_block = i*nxb;
        size_t imax_block = std::min(static_cast<int>((i+1)*nxb), static_cast<int>(nx));
        size_t jmin_block = j*nyb;
        size_t jmax_block = std::min(static_cast<int>((j+1)*nyb), static_cast<int>(ny));
        double P_est;

        DataInBlockBinarySearch(P_est, P_min, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, true, true, true, true, min_data_in_range, kriging_data, 0.05);

        size_t imin_range = std::max(static_cast<int>(std::floor(static_cast<double>(i*nxb)-P_est*static_cast<double>(rx))), 0);
        size_t imax_range = std::min(static_cast<int>(std::ceil(static_cast<double>((i+1)*nxb)+P_est*static_cast<double>(rx))), static_cast<int>(nx));
        size_t jmin_range = std::max(static_cast<int>(std::floor(static_cast<double>(j*nyb)-P_est*static_cast<double>(ry))), 0);
        size_t jmax_range = std::min(static_cast<int>(std::ceil(static_cast<double>((j+1)*nyb)+P_est*static_cast<double>(ry))), static_cast<int>(ny));

        kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_range, imax_range, jmin_range, jmax_range);

        size_t imin_new, imax_new, jmin_new, jmax_new;



        if (binary_search_dir){
          double P_dir;
          size_t index_block;
          //-----------------------------------------
          // y = j
          //-----------------------------------------
          jmin = j*nyb;
          if (j == n_blocks_y - 1)
            jmax = ny;
          else
            jmax = (j+1)*nyb;
          // x = i-1
          if (i > 0){
            index_block = (i-1)*n_blocks_y + j;
            imax = i*nxb;
            imin = std::max(static_cast<int>(i*nxb - P_est*rx), 0);
            if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
              DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, true, false, false, false, min_data_in_block, kriging_data, 0.05);
              imin_new = std::max(static_cast<int>(i*nxb - P_dir*rx), 0);
              if (imin_new < imin_range)
                kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_new, imin_range, jmin_block, jmax_block);
            }
          }
          // x = i+1
          jmin = j*nyb;
          if (j == n_blocks_y - 1)
            jmax = ny;
          else
            jmax = (j+1)*nyb;
          if ((i+1) < n_blocks_x){
            index_block = (i+1)*n_blocks_y + j;
            imin = (i+1)*nxb;
            if ((i+1) == (n_blocks_x-1))
              imax = nx;
            else
              imax = std::min((static_cast<int>((i+1)*nxb+(P_est*rx))), static_cast<int>(nx));
            if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
              DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, false, true, false, false, min_data_in_block, kriging_data, 0.05);
              imax_new = std::min((static_cast<int>((i+1)*nxb+P_dir*rx)), static_cast<int>(nx));
              if (imax_new > imax_range)
                kriging_data.AddDataToBlock(*kriging_data_blocks[index], imax_range, imax_new, jmin_block, jmax_block);
            }
          }
          //-----------------------------------------
          //y = j+1
          //-----------------------------------------
          if (j < (n_blocks_y-1)){
            jmin = (j+1)*nyb;
            // x = i-1
            if (i > 0){
              index_block = (i-1)*n_blocks_y + j+1;
              imax = i*nxb;
              imin = std::max((static_cast<int>(i*nxb-(P_est*rx))), 0);
              if ((j+1) == (n_blocks_y-1))
                jmax = ny;
              else
                jmax = std::min((static_cast<int>((j+1)*nyb+(P_est*ry))), static_cast<int>(ny));
              if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
                DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, true, false, false, true, min_data_in_block, kriging_data, 0.05);
                imin_new = std::max((static_cast<int>(i*nxb-P_dir*rx)), 0);
                jmax_new = std::min((static_cast<int>((j+1)*nyb+P_dir*ry)), static_cast<int>(ny));
                if (imin_range > imin_new)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_new,   imin_range, jmax_block, jmax_new);
                if (jmax_new > jmax_range)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_range, imin_block, jmax_range, jmax_new);
              }
            }
            // x = i
            jmin = (j+1)*nyb;
            index_block = i*n_blocks_y + j+1;
            imin = i*nxb;
            if (i == n_blocks_x - 1)
              imax = nx;
            else
              imax = (i+1)*nxb;
            if ((j+1) == (n_blocks_y-1))
              jmax = ny;
            else
              jmax = std::min((static_cast<int>((j+1)*nyb+P_est*ry)), static_cast<int>(ny));
            if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
              DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, false, false, false, true, min_data_in_block, kriging_data, 0.05);
              jmax_new = std::min((static_cast<int>((j+1)*nyb+P_dir*rx)), static_cast<int>(ny));
              if (jmax_new > jmax_range)
                kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_block, imax_block, jmax_range, jmax_new);
            }
            // x = i+1
            jmin = (j+1)*nyb;
            if ((i+1) < n_blocks_x){
              index_block = (i+1)*n_blocks_y + j+1;
              imin = (i+1)*nxb;
              if ((j+1) == (n_blocks_y-1))
                jmax = ny;
              else
                jmax = std::min((static_cast<int>((j+1)*nyb+P_est*ry)), static_cast<int>(ny));
              if ((i+1) == (n_blocks_x-1))
                imax = nx;
              else
                imax = std::min((static_cast<int>((i+1)*nxb+P_est*rx)), static_cast<int>(nx));
              if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
                DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, false, true, false, true, min_data_in_block, kriging_data, 0.05);
                imax_new = std::min((static_cast<int>((i+1)*nxb+P_dir*rx)), static_cast<int>(nx));
                jmax_new = std::min((static_cast<int>((j+1)*nyb+P_dir*ry)), static_cast<int>(ny));
                if (imax_new > imax_range)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imax_range, imax_new,   jmax_block, jmax_new);
                if (jmax_new > jmax_range)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imax_block, imax_range, jmax_range, jmax_new);
              }
            }
          }
          //-----------------------------------------
          // y = j-1
          //-----------------------------------------
          if (j > 0){
            jmax = j*nyb;
            // x = i-1
            if (i > 0){
              index_block = (i-1)*n_blocks_y + j-1;
              imax = i*nxb;
              jmin = std::max((static_cast<int>(j*nyb-P_est*ry)), 0);
              imin = std::max((static_cast<int>(i*nxb-P_est*rx)), 0);
              if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
                DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, true, false, true, false, min_data_in_block, kriging_data, 0.05);
                jmin_new = std::max((static_cast<int>(j*nyb-P_dir*rx)), 0);
                imin_new = std::max((static_cast<int>(i*nxb-P_dir*ry)), 0);
                if (imin_range > imin_new)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_new,   imin_range, jmin_new, jmin_block);
                if (jmin_range > jmin_new)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_range, imin_block, jmin_new, jmin_range);
              }
            }
            // x = i
            jmax = j*nyb;
            index_block = i*n_blocks_y + j-1;
            imin = i*nxb;
            if (i == n_blocks_x - 1)
              imax = nx;
            else
              imax = (i+1)*nxb;
            jmin = std::max((static_cast<int>(j*nyb-P_est*ry)), 0);
            if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
              DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, false, false, true, false, min_data_in_block, kriging_data, 0.05);
              jmin_new = std::max((static_cast<int>(j*nyb-P_dir*ry)), 0);
              if (jmin_range > jmin_new)
                kriging_data.AddDataToBlock(*kriging_data_blocks[index], imin_block, imax_block, jmin_new, jmin_range);
            }
            // x = i+1
            jmax = j*nyb;
            if ((i+1) < n_blocks_x){
              index_block = (i+1)*n_blocks_y + j-1;
              imin = (i+1)*nxb;
              jmin = std::max((static_cast<int>(j*nyb-P_est*ry)), 0);
              if ((i+1) == (n_blocks_x-1))
                imax = nx;
              else
                imax = std::min((static_cast<int>((i+1)*nxb+P_est*rx)), static_cast<int>(nx));
              if (kriging_data.CountDataInBlock(imin, imax, jmin, jmax) < min_data_in_block){
                DataInBlockBinarySearch(P_dir, P_est, P_max, rx, ry, nx, ny, imin_block, imax_block, jmin_block, jmax_block, false, true, true, false, min_data_in_block, kriging_data, 0.05);
                jmin_new = std::max((static_cast<int>(j*nyb -    P_dir*ry)), 0);
                imax_new = std::min((static_cast<int>((i+1)*nxb + P_dir*rx)), static_cast<int>(nx));
                if (imax_new > imax_range)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imax_range, imax_new,   jmin_new, jmin_block);
                if (jmin_range > jmin_new)
                  kriging_data.AddDataToBlock(*kriging_data_blocks[index], imax_block, imax_range, jmin_new, jmin_range);
              }
            }
          }
        }
        P_est_all += P_est;
        if (P_est < P_min_total) P_min_total = P_est;
        if (P_est > P_max_total) P_max_total = P_est;
      }
    }
  }
  else // only one block
    kriging_data_blocks[0] = kriging_data.Clone();  // Pointers deleted in Krig2D()

  //Return average of P_est from all blocks:
  P = P_est_all * (1/static_cast<double>(n_blocks));


  // sort blocks after largest number of data, for parallelizing
  typedef std::pair<int, int> Pair;
  std::vector<Pair> nums;

  size_t avg_data_in_blocks = 0;

  NodalKrigingData2D * nodal_class = new NodalKrigingData2D();

  for (size_t i = 0 ; i < n_blocks ; ++i) {
    if(strcmp(typeid(*kriging_data_blocks[i]).name(), typeid(*nodal_class).name()) == 0) {
      NodalKrigingData2D & data = dynamic_cast<NodalKrigingData2D &>(*kriging_data_blocks[i]);
      if(&data != 0) {
        data.FindMeanValues();
      }
    }
    nums.push_back(Pair(static_cast<int>(i), kriging_data_blocks[i]->GetNumberOfData()));
    avg_data_in_blocks += kriging_data_blocks[i]->GetNumberOfData();
  }

  n_avg_data_blocks = static_cast<size_t>(std::floor(static_cast<double>(avg_data_in_blocks)/static_cast<double>(n_blocks)));
  delete nodal_class;

  std::sort(nums.begin(), nums.end(), compare);

  std::vector<KrigingData2D *>     kriging_data_blocks_sort;
  std::vector<size_t>              block_x_sort;
  std::vector<size_t>              block_y_sort;

  for (size_t i = 0 ; i < n_blocks ; ++i) {
    kriging_data_blocks_sort.push_back( kriging_data_blocks[nums[i].first] );
    block_x_sort.push_back( block_index_x[nums[i].first] );
    block_y_sort.push_back( block_index_y[nums[i].first] );
  }
  kriging_data_blocks = kriging_data_blocks_sort;
  block_index_x = block_x_sort;
  block_index_y = block_y_sort;
}

//----------------------------------------------------------------
void Kriging::Krig2DBlock(const Grid2D<double>      & trend_orig,
                          const NodalKrigingData2D  & kriging_data,
                          const CovGrid2D           & cov,
                          size_t                      block_x,
                          size_t                      block_y,
                          size_t                      nxb,
                          size_t                      nyb,
                          size_t                      n_blocks_x,
                          size_t                      n_blocks_y,
                          Grid2D<double>            & trend,
                          Grid2D<double>            & filled,
                          bool                        get_residuals)
//----------------------------------------------------------------
{
  //
  // This routine returns z(x) = m(x) + k(x)K^{-1}(d - m).
  //
  int                      md     = kriging_data.GetNumberOfData();
  const std::vector<int> & indexi = kriging_data.GetIndexI();
  const std::vector<int> & indexj = kriging_data.GetIndexJ();
  std::vector<float>       data   = kriging_data.GetData();   // Take an editable copy

  int                      nx     = static_cast<int>(trend.GetNI());
  int                      ny     = static_cast<int>(trend.GetNJ());


  // find min and max of block
  int imin = std::max(static_cast<int>(block_x*nxb), 0);
  int imax;
  if (block_x == (n_blocks_x-1))
    imax = nx;
  else
    imax = std::min(static_cast<int>((block_x+1)*nxb),nx);

  int jmin = std::max(static_cast<int>(block_y*nyb), 0);
  int jmax;
  if (block_y == (n_blocks_y-1))
    jmax = ny;
  else
    jmax = std::min(static_cast<int>((block_y+1)*nyb),ny);


  if (md > 0 && md <= nx*ny) {

    SymmetricMatrix K(md);
    Vector          residual(md);
    Vector          k(md);
    Vector          x(md);

    // subtract trend from original trend grid
    for (int i = 0 ; i < md ; i++)
      residual(i) = data[i] - static_cast<float>(trend_orig(indexi[i],indexj[i]));

    FillKrigingMatrix(K, cov, indexi, indexj);

    CholeskySolve(K, residual, x);

    for(size_t i=0; i < static_cast<size_t>(md); i++) {
      if (indexi[i] < imax)
        if (indexi[i] >= imin)
          if (indexj[i] < jmax)
            if (indexj[i] >= jmin) {
              if(get_residuals)
                trend(indexi[i],indexj[i]) = residual(static_cast<int>(i));
              else
                trend(indexi[i],indexj[i]) += residual(static_cast<int>(i));
              filled(indexi[i],indexj[i]) = 1.0;
            }
    }

    for (int i = imin ; i < imax ; i++) {
      for (int j = jmin ; j < jmax ; j++) {
        if(!(filled(i,j) > 0.0)){ // if this is not a datapoint
          FillKrigingVector(k, cov, indexi, indexj, i, j);
          if (get_residuals)
            trend(i,j) = k * x;
          else
            trend(i,j) += k * x;
        }
      }
    }
  }
}

//---------------------------------------------------------------
void Kriging::FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj)
//---------------------------------------------------------------
{
  int n = K.dim();
  //Grid2D<double> K_test(n, n, 0); //for debug print

  for(int i=0 ; i < n  ; i++) {
    for(int j=0 ; j <= i ; j++) {
      int deltai = indexi[i] - indexi[j];
      int deltaj = indexj[i] - indexj[j];
      K(j,i) = static_cast<double>(cov.GetCov(deltai,deltaj));
      //K_test(j,i) = static_cast<double>(cov.GetCov(deltai,deltaj)); //for debug print
    }
  }
  //RegularSurface<double> K_surface(0,0,static_cast<double>(n-1),static_cast<double>(n-1),K_test); //for debug print
  //std::string filename = "K_surface.irap";
  //WriteIrapClassicAsciiSurf(K_surface, 0.0, filename);
}

//--------------------------------------------------------------
void Kriging::FillKrigingVector(NRLib::Vector          & k,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj,
                                int i, int j)
//--------------------------------------------------------------
{
  for(int ii=0 ; ii < k.length() ; ii++) {
    int deltai = indexi[ii] - i;
    int deltaj = indexj[ii] - j;
    k(ii) = static_cast<double>(cov.GetCov(deltai,deltaj));
  }
}

//----------------------------------------------------------------
void Kriging::Krig2DBlock(const RegularSurface<double>   & trend_orig,
                          const ContinuousKrigingData2D  & kriging_data,
                          const Variogram                & var,
                          size_t                           block_x,
                          size_t                           block_y,
                          size_t                           nxb,
                          size_t                           nyb,
                          size_t                           n_blocks_x,
                          size_t                           n_blocks_y,
                          RegularSurface<double>         & trend,
                          bool                             get_residuals)
//----------------------------------------------------------------
{
  //
  // This routine returns z(x) = m(x) + k(x)K^{-1}(d - m).
  //

  int                      md     = kriging_data.GetNumberOfData();
  const std::vector<double> & xvector = kriging_data.GetXVector();
  const std::vector<double> & yvector = kriging_data.GetYVector();
  std::vector<double>       data   = kriging_data.GetData();   // Take an editable copy

  int                      nx     = static_cast<int>(trend.GetNI());
  int                      ny     = static_cast<int>(trend.GetNJ());

  double dx = trend.GetDX();
  double dy = trend.GetDY();
  double xmin = trend.GetXMin();
  double ymin = trend.GetYMin();

  // find min and max of block
  int imin = std::max(static_cast<int>(block_x*nxb), 0);
  int imax;
  if (block_x == (n_blocks_x-1))
    imax = nx;
  else
    imax = std::min(static_cast<int>((block_x+1)*nxb),nx);

  int jmin = std::max(static_cast<int>(block_y*nyb), 0);
  int jmax;
  if (block_y == (n_blocks_y-1))
    jmax = ny;
  else
    jmax = std::min(static_cast<int>((block_y+1)*nyb),ny);


  if (md > 0 && md <= nx*ny) {

    SymmetricMatrix K(md);
    Vector          residual(md);
    Vector          k(md);
    Vector          x(md);


    // NB dette kunne vært gjort bare en gang for alle blockene. Dette er blockuavhening
    // subtract trend from original trend grid
      for (int i = 0 ; i < md ; i++)
        residual(i) = data[i] - static_cast<float>(trend_orig.GetZ(xvector[i], yvector[i]));

    FillKrigingMatrix(K, var, xvector, yvector);
    CholeskySolve(K, residual, x);
    for (int i = imin ; i < imax ; i++) {
      for (int j = jmin ; j < jmax ; j++) {
        double x_coordinate = xmin + i*dx;
        double y_coordinate = ymin + j*dy;
        FillKrigingVector(k, var, xvector, yvector, x_coordinate, y_coordinate);
        if (get_residuals)
          trend(i,j) = k * x;
        else
          trend(i,j) += k * x;
      }
    }
  }
}


//---------------------------------------------------------------
void Kriging::FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const Variogram        & var,
                                const std::vector<double> & xvector,
                                const std::vector<double> & yvector)
  //---------------------------------------------------------------
{
  int n = K.dim();
  //Grid2D<double> K_test(n, n, 0); //for debug print

  for(int i=0 ; i < n  ; i++) {
    for(int j=0 ; j <= i ; j++) {
      double deltax = xvector[i] - xvector[j];
      double deltay = yvector[i] - yvector[j];
      K(j,i) = var.GetCov(deltax, deltay);
      //K_test(j,i) = K(j,i);

    }
  }
  //RegularSurface<double> K_surface(0,0,n-1,n-1,K_test); //for debug print
  //std::string filename = "K_surface3.irap";
  //WriteIrapClassicAsciiSurf(K_surface, 0.0, filename);
}

//--------------------------------------------------------------
void Kriging::FillKrigingVector(NRLib::Vector             & k,
                                const Variogram           & var,
                                const std::vector<double> & xvector,
                                const std::vector<double> & yvector,
                                double x, double y)
  //--------------------------------------------------------------
{
  for(int i=0 ; i < k.length() ; i++) {
    double deltax = xvector[i] - x;
    double deltay = yvector[i] - y;
    k(i) = var.GetCov(deltax, deltay);
  }
}

