#include "extrapolategrid2d.hpp"

#include "../iotools/fileio.hpp"

#include <iomanip>
#include <fstream>

#include <math.h>

using namespace NRLib;

//---------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::InverseDistanceWeightingExtrapolation(Grid2D<double>     & grid,
                                                              const Grid2D<bool> & mask, // Extrapolate where grid cells are true
                                                              const double         x0,
                                                              const double         y0,
                                                              const double         xinc,
                                                              const double         yinc,
                                                              const double         missing)
//---------------------------------------------------------------------------------------------
{
  std::vector<std::pair<size_t, size_t> > edge_indices;
  std::vector<std::pair<size_t, size_t> > inside_indices;
  std::vector<std::pair<size_t, size_t> > regular_inside_indices;
  std::vector<std::pair<size_t, size_t> > missing_indices;
  std::vector<std::pair<size_t, size_t> > stationary_indices;
  std::vector<std::pair<size_t, size_t> > control_indices;

  ClassifyPoints(edge_indices,
                 inside_indices,
                 missing_indices,
                 stationary_indices,
                 regular_inside_indices,
                 control_indices,
                 grid,
                 mask,
                 missing);

  // Do extrapolation with normed inverse distance
  double power = 1.0;

  for (size_t k = 0; k < missing_indices.size(); ++k) {
    size_t i        = missing_indices[k].first;
    size_t j        = missing_indices[k].second;

    double z_interp = 0.0;
    double norm     = 0.0;
    double d2;

    for (size_t p = 0; p < control_indices.size(); ++p) {
      int    ip     = control_indices[p].first;
      int    jp     = control_indices[p].second;
      int    dipx   = (static_cast<int>(i) - ip)*xinc;
      int    djpy   = (static_cast<int>(j) - jp)*yinc;
      double d2     = dipx*dipx + djpy*djpy;
      //double dist   = std::sqrt(d2);
      //double weight = std::pow(dist, -power);
      double weight = 1.0/d2;       // Use inverse square distance for simplicity
      double zp     = grid(ip, jp);
      z_interp     += zp * weight;
      norm         += weight;
    }
    grid(i, j) = z_interp/norm;
  }

  DumpResults(edge_indices,
              stationary_indices,
              regular_inside_indices,
              control_indices,
              x0,
              y0,
              xinc,
              yinc);
}

//-------------------------------------------------------------------------------------------------------
void ExtrapolateGrid2D::ClassifyPoints(std::vector<std::pair<size_t, size_t> > & edge_indices,
                                       std::vector<std::pair<size_t, size_t> > & inside_indices,
                                       std::vector<std::pair<size_t, size_t> > & missing_indices,
                                       std::vector<std::pair<size_t, size_t> > & stationary_indices,
                                       std::vector<std::pair<size_t, size_t> > & regular_inside_indices,
                                       std::vector<std::pair<size_t, size_t> > & control_indices,
                                       const Grid2D<double>                    & grid,
                                       const Grid2D<bool>                      & mask,
                                       const double                              missing)
//-------------------------------------------------------------------------------------------------------
{
  edge_indices.clear();
  inside_indices.clear();
  missing_indices.clear();
  stationary_indices.clear();

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

  // Find stationary points

  for (size_t i = 1; i < nx - 1; ++i) {
    for (size_t j = 1; j < ny - 1; ++j) {
      if ((dx[i - 1][j] < 0.0) != (dx[i][j] < 0.0)  && (dy[j - 1][i] < 0.0) != (dy[j][i] < 0.0))
        stationary_indices.push_back(std::pair<size_t, size_t>(i, j));
    }
  }

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
