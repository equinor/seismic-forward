#include "extrapolategrid2d.hpp"

#include "../iotools/fileio.hpp"

#include <iomanip>
#include <fstream>

#include <math.h>

using namespace NRLib;

//-------------------------------------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::InverseDistanceWeightingExtrapolation(NRLib::Grid2D<double>                         & grid,
                                                              const std::vector<std::pair<size_t, size_t> > & miss_indices,
                                                              const std::vector<std::pair<size_t, size_t> > & data_indices,
                                                              const double                                    xinc,
                                                              const double                                    yinc)
//-------------------------------------------------------------------------------------------------------------------------
{
  int n_data = std::min(3, static_cast<int>(data_indices.size()));  // Require at least 3 grid cells to do extrapolation

  for (size_t k = 0 ; k < miss_indices.size() ; k++) {
    int i     = miss_indices[k].first;
    int j     = miss_indices[k].second;

    int imax  = 5;  // Start looking for data +- 9 grid cells in i-direction (5*2 - 1)
    int jmax  = 5;  // Start looking for data +- 9 grid cells in j-direction (5*2 - 1)

    //
    // Find search radius with enough data and collect indices
    //
    std::vector<size_t> indices;

    while (indices.size() < n_data) {
      indices.clear();

      imax *= 2; // Double search radius
      jmax *= 2; // Double search radius

      for (size_t p = 0 ; p < data_indices.size() ; p++) {
        int ip = data_indices[p].first;
        int jp = data_indices[p].second;

        if ((std::abs(i - ip) < imax) && (std::abs(j - jp) < jmax)) {
          indices.push_back(p);
        }
      }
    }

    //
    // Extrapolate
    //
    if (indices.size() > 0) {
      double z    = 0.0;
      double norm = 0.0;

      for (size_t p = 0 ; p < indices.size() ; p++) {
        int    ip     = data_indices[indices[p]].first;
        int    jp     = data_indices[indices[p]].second;
        double dipx   = static_cast<double>(i - ip)*xinc;
        double djpy   = static_cast<double>(j - jp)*yinc;
        double d2     = dipx*dipx + djpy*djpy;
        double weight = 1.0 / d2;       // Use inverse square distance for simplicity
        double zp     = grid(ip, jp);
        z            += zp * weight;
        norm         += weight;
      }
      grid(i, j) = z/norm;
    }
  }
}

//-------------------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::ClassifyPoints(std::vector<std::pair<size_t, size_t> > & missing_indices,
                                       std::vector<std::pair<size_t, size_t> > & control_indices,
                                       const Grid2D<double>                    & grid,
                                       const Grid2D<bool>                      & mask,
                                       const double                              missing)
//-------------------------------------------------------------------------------------------------------
{
  std::vector<std::pair<size_t, size_t> > edge_indices;
  std::vector<std::pair<size_t, size_t> > inside_indices;
  std::vector<std::pair<size_t, size_t> > regular_inside_indices;
  std::vector<std::pair<size_t, size_t> > stationary_indices;

  missing_indices.clear();
  control_indices.clear();

  // Find edges
  size_t nx = grid.GetNI();
  size_t ny = grid.GetNJ();

  size_t inside_imin = 0;
  size_t inside_imax = nx;
  size_t inside_jmin = 0;
  size_t inside_jmax = ny;

  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      if (mask(i, j)) {
        if (grid(i, j) == missing)
          missing_indices.push_back(std::pair<size_t, size_t>(i, j));
        else if (grid.IsEdge(i, j, missing))
          edge_indices.push_back(std::pair<size_t, size_t>(i, j));
        else {
          inside_indices.push_back(std::pair<size_t, size_t>(i, j));
          inside_imin = std::min<size_t>(i, inside_imin);
          inside_imax = std::max<size_t>(i, inside_imax);
          inside_jmin = std::min<size_t>(j, inside_jmin);
          inside_jmax = std::max<size_t>(j, inside_jmax);
        }
      }
    }
  }

  // Find gradient on indices inside

  std::vector<std::vector<double> > dx(nx - 1, std::vector<double>(ny, missing));
  std::vector<std::vector<double> > dy(ny - 1, std::vector<double>(nx, missing));

  for (size_t p = 0; p < inside_indices.size(); ++p) {
    size_t i = inside_indices[p].first;
    size_t j = inside_indices[p].second;
    dx[i][j] = grid(i + 1, j) - grid(i, j);
    dy[j][i] = grid(i, j + 1) - grid(i, j);
  }

  /*
  // Find stationary points

  for (size_t i = 1; i < nx - 1; ++i) {
    for (size_t j = 1; j < ny - 1; ++j) {
      if ((dx[i - 1][j] < 0.0) != (dx[i][j] < 0.0)  && (dy[j - 1][i] < 0.0) != (dy[j][i] < 0.0))
        stationary_indices.push_back(std::pair<size_t, size_t>(i, j));
    }
  }
   */

  // Final control indices

  control_indices = edge_indices;
  if (stationary_indices.empty()) { // This may occur, need some points inside anyway. Just pick the middle point
    size_t imid = (inside_imax + inside_imin) / 2;
    size_t jmid = (inside_jmax + inside_jmin) / 2;
    for (size_t p = 0; p < inside_indices.size(); ++p) {
      size_t i = inside_indices[p].first;
      size_t j = inside_indices[p].second;
      if (i == imid && j == jmid)
       regular_inside_indices.push_back(std::pair<size_t, size_t>(i, j));
    }
    // Join edge indices and stationary indices (if any) to form control points
    control_indices.insert(control_indices.begin(), regular_inside_indices.begin(), regular_inside_indices.end());
  }
  else {
    control_indices.insert(control_indices.begin(), stationary_indices.begin(), stationary_indices.end());
  }

  /*
  DumpResults(edge_indices,
              stationary_indices,
              regular_inside_indices,
              control_indices,
              x0,
              y0,
              xinc,
              yinc);
              */
}



//---------------------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::DumpResults(const std::vector<std::pair<size_t, size_t> > & edge_indices,
                                    const std::vector<std::pair<size_t, size_t> > & stationary_indices,
                                    const std::vector<std::pair<size_t, size_t> > & regular_inside_indices,
                                    const std::vector<std::pair<size_t, size_t> > & control_indices,
                                    double                                          x0,
                                    double                                          y0,
                                    double                                          xinc,
                                    double                                          yinc)
//---------------------------------------------------------------------------------------------------------
{
  // Edge
  std::string outfile = "extrapolation_points_edges.xyz";
  WritePointsToFile(outfile, edge_indices, x0, y0, xinc, yinc);

  // Stationary
  outfile =  "extrapolation_points_stationary.xyz";
  WritePointsToFile(outfile, stationary_indices, x0, y0, xinc, yinc);

  // Regular inside
  outfile =  "extrapolation_points_regular_inside.xyz";
  WritePointsToFile(outfile, regular_inside_indices, x0, y0, xinc, yinc);

  // Control points
  outfile =  "extrapolation_points_control.xyz";
  WritePointsToFile(outfile, control_indices, x0, y0, xinc, yinc);
}

//-------------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::WritePointsToFile(const std::string                             & filename,
                                          const std::vector<std::pair<size_t, size_t> > & indices,
                                          double                                          x0,
                                          double                                          y0,
                                          double                                          xinc,
                                          double                                          yinc)
//-------------------------------------------------------------------------------------------------
{
  if (indices.empty())
    return;

  std::ofstream file;
  OpenWrite(file, filename);

  for (size_t k = 0; k < indices.size(); ++k) {
    size_t i = indices[k].first;
    size_t j = indices[k].second;
    double x = x0 + i*xinc;
    double y = y0 + j*yinc;

    file << std::fixed
         << std::setprecision(2)
         << x << " " << y << " " << 0.00 << "\n";

  }
  file.close();
}
