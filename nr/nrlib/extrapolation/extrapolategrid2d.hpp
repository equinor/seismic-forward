#ifndef NRLIB_EXTRAPOLATEGRID2D_H
#define NRLIB_EXTRAPOLATEGRID2D_H

#include "../grid/grid2d.hpp"

namespace NRLib {

  namespace ExtrapolateGrid2D
  {
    void   InverseDistanceWeightingExtrapolation(Grid2D<double>     & grid,
                                                 const Grid2D<bool> & extrapolate,
                                                 const double         x0,
                                                 const double         y0,
                                                 const double         xinc,
                                                 const double         yinc,
                                                 const double         missing);

    void   ClassifyPoints(std::vector<std::pair<size_t, size_t> > & edge_indices,
                          std::vector<std::pair<size_t, size_t> > & inside_indices,
                          std::vector<std::pair<size_t, size_t> > & outside_indices,
                          std::vector<std::pair<size_t, size_t> > & stationary_indices,
                          std::vector<std::pair<size_t, size_t> > & regular_inside_indices,
                          std::vector<std::pair<size_t, size_t> > & control_indices,
                          const Grid2D<double>                    & grid,
                          const Grid2D<bool>                      & extrapolate,
                          const double                              missing);

    double FindSquareDistanceInEclipseGrid(const int    dip,
                                           const int    djp,
                                           const double xinc,
                                           const double yinc,
                                           const double cosA,
                                           const double sinA);

    void   DumpResults(const std::vector<std::pair<size_t, size_t> > & edge_indices,
                       const std::vector<std::pair<size_t, size_t> > & stationary_indices,
                       const std::vector<std::pair<size_t, size_t> > & regular_inside_indices,
                       const std::vector<std::pair<size_t, size_t> > & control_indices,
                       double                                          x0,
                       double                                          y0,
                       double                                          xinc,
                       double                                          yinc);

    void   WritePointsToFile(const std::string                             & outfile,
                             const std::vector<std::pair<size_t, size_t> > & indices,
                             double                                          x0,
                             double                                          y0,
                             double                                          xinc,
                             double                                          yinc);

  }
}

#endif // NRLIB_EXTRAPOLATEGRID2D_END
