// $Id: kriging.hpp 1301 2015-06-02 13:51:17Z vigsnes $

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

#ifndef NRLIB_STATISTICS_KRIGING_HPP
#define NRLIB_STATISTICS_KRIGING_HPP

#include <cstdlib>
#include <vector>

#include "../flens/nrlib_flens.hpp"
#include "../variogram/variogram.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/variogram/covgrid2d.hpp"
#include "nrlib/variogram/variogram.hpp"
#include "nrlib/statistics/krigingdata2d.hpp"
#include "nrlib/statistics/continuouskrigingdata2d.hpp"
#include "nrlib/statistics/nodalkrigingdata2d.hpp"
#include "nrlib/geometry/point.hpp"
#include "nrlib/surface/regularsurface.hpp"


namespace NRLib {
class Kriging
{
public:
  static void Krig1D(std::vector<double> &field,
              const std::vector<bool>    &is_known,
              const std::vector<double>  &obs,
              double                      dx,
              const Variogram            &vario);

  static void Krig2D(Grid2D<double>           & trend,
                     const NodalKrigingData2D & krigingData,
                     const CovGrid2D          & cov,
                     bool                       get_residuals = false,
                     bool                       krig_all = false,
                     double                     P = 2.0);

  static void Krig2D(RegularSurface<double>     & trend,
                     const std::vector<Point* > & data_points,
                     const Variogram            & var,
                     bool                         get_residuals  = false,
                     bool                         use_index_grid = false,
                     bool                         krig_all = false,
                     double                       P = 2.0);

private:
  static void Krig2D(Grid2D<double>             & trend,
                     RegularSurface<double>     & trend_surface,
                     const KrigingData2D        & kriging_data,
                     const CovGrid2D            & cov,
                     const Variogram            & var,
                     bool                         get_residuals,
                     double                       dx,
                     double                       dy,
                     double                       Rx,
                     double                       Ry,
                     double                       P,
                     bool                         krig_all = false);


  static void DataInBlockBinarySearch(double              & P,
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
                                      const double          tol);

    static void FindOptimalDataInBlock(const Grid2D<double>         & trend,
                                       const KrigingData2D          & kriging_data,
                                       double                         Rx,
                                       double                         Ry,
                                       double                         dx,
                                       double                         dy,
                                       size_t                         n_blocks_x,
                                       size_t                         n_blocks_y,
                                       double                         P,
                                       size_t                         n_avg_data_blocks,
                                       std::vector<KrigingData2D *> & kriging_data_blocks,
                                       std::vector<size_t>           & block_index_x,
                                       std::vector<size_t>           & block_index_y);


  static void FindOptimalBlock(const Grid2D<double> & trend,
                               double                 Rx,
                               double                 Ry,
                               double                 dx,
                               double                 dy,
                               size_t                 n_data,
                               size_t                 n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y);

  static void FindOptimalBlockP(const Grid2D<double> & trend,
                               double                 Rx,
                               double                 Ry,
                               double                 dx,
                               double                 dy,
                               size_t                 n_data,
                               size_t                 n_threads,
                               size_t               & n_blocks_x,
                               size_t               & n_blocks_y,
                               double               & P,
                               size_t               & n_data_blocks);

  static double OverallTime(size_t n,
                            size_t Ns,
                            size_t nx,
                            size_t ny);

  static void SetConstants(size_t & max_data_in_range,
                           size_t & min_data_in_block,
                           size_t & min_data_in_range,
                           size_t   nxb,
                           size_t   nyb,
                           size_t   rx,
                           size_t   ry,
                           double   P);

  static void Krig2DBlock(const Grid2D<double>      & trend_orig,
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
                          bool                        get_residuals);

  static void FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj);

  static void FillKrigingVector(NRLib::Vector          & k,
                                const CovGrid2D        & cov,
                                const std::vector<int> & indexi,
                                const std::vector<int> & indexj,
                                int i,
                                int j);

  static void Krig2DBlock(const RegularSurface<double>   & trend_orig,
                          const ContinuousKrigingData2D  & kriging_data,
                          const Variogram                & var,
                          size_t                           block_x,
                          size_t                           block_y,
                          size_t                           nxb,
                          size_t                           nyb,
                          size_t                           n_blocks_x,
                          size_t                           n_blocks_y,
                          RegularSurface<double>         & trend,
                          bool                             get_residuals);

  static void FillKrigingMatrix(NRLib::SymmetricMatrix & K,
                                const Variogram        & var,
                                const std::vector<double> & xvector,
                                const std::vector<double> & yvector);

  static void FillKrigingVector(NRLib::Vector          & k,
                                const Variogram        & var,
                                const std::vector<double> & xvector,
                                const std::vector<double> & yvector,
                                double x,
                                double y);
};
}
#endif

