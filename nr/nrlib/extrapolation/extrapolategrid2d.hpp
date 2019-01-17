#ifndef NRLIB_EXTRAPOLATEGRID2D_H
#define NRLIB_EXTRAPOLATEGRID2D_H

#include "../grid/grid2d.hpp"

namespace NRLib {

  namespace ExtrapolateGrid2D
  {
    void ExtrapolateLayer(NRLib::Grid2D<double>                         & grid,
                          const std::vector<std::pair<size_t, size_t> > & data_indices,
                          const std::vector<std::pair<size_t, size_t> > & miss_indices,
                          const double                                    xinc,
                          const double                                    yinc);

    void   InverseDistanceWeightingExtrapolation(Grid2D<double>                                & grid,
                                                 const std::vector<std::pair<size_t, size_t> > & miss_indices,
                                                 const std::vector<std::pair<size_t, size_t> > & data_indices,
                                                 const double                                    xinc,
                                                 const double                                    yinc);

    void   ClassifyPoints(std::vector<std::pair<size_t, size_t> > & outside_indices,
                          std::vector<std::pair<size_t, size_t> > & control_indices,
                          const Grid2D<double>                    & grid,
                          const Grid2D<bool>                      & mask,
                          const double                              missing);

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
