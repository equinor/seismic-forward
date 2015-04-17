// $Id: eclipsecarfingroup.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "eclipsecarfingroup.hpp"

#include "eclipsecarfin.hpp"

#include <fstream>
#include <vector>


using namespace NRLib;

EclipseCarFinGroup::EclipseCarFinGroup(const std::string& name)
  : name_(name), refinement_x_(0), refinement_y_(0), refinement_z_(0)
{}


void EclipseCarFinGroup::WriteToFile(std::ofstream& of) const
{
  std::vector<const EclipseCarFin*>::const_iterator it = cartesian_refinements_.begin();
  size_t line_len = 0;
  /// @todo Move to a suitable location.
  const size_t max_line_len = 77;
  of << "  ";
  for ( ; it < cartesian_refinements_.end(); ++it) {
    std::string name = (*it)->GetName();
    if (line_len + name.size() + 3 > max_line_len) {
      of << "\n  ";
      line_len = 0;
    }
    of << "\'" << name << "\' ";
    line_len += name.size() + 3;
  }
  if (line_len + 3 > max_line_len)
    of << "\n  ";
  of << "/\n\n";
}


void EclipseCarFinGroup::AddCarFin(const EclipseCarFin* cartesian_refinement)
{
  size_t refinement_x, refinement_y, refinement_z;
  cartesian_refinement->GetRefinements(refinement_x, refinement_y, refinement_z);
  if (cartesian_refinements_.size() == 0) {
    refinement_x_ = refinement_x;
    refinement_y_ = refinement_y;
    refinement_z_ = refinement_z;
  } else {
    assert(refinement_x == refinement_x_ &&
           refinement_y == refinement_y_ &&
           refinement_z == refinement_z_    );
  }
  cartesian_refinements_.push_back(cartesian_refinement);
}


void EclipseCarFinGroup::GetLocation(size_t& i0, size_t& i1,
                                     size_t& j0, size_t& j1,
                                     size_t& k0, size_t& k1) const
{
  std::vector<const EclipseCarFin*>::const_iterator it = cartesian_refinements_.begin();
  (*it)->GetLocation(i0, i1, j0, j1, k0, k1);
  ++it;
  size_t tmpi0, tmpi1, tmpj0, tmpj1, tmpk0, tmpk1;
  for ( ; it < cartesian_refinements_.end(); ++it) {
    (*it)->GetLocation(tmpi0, tmpi1, tmpj0, tmpj1, tmpk0, tmpk1);
    if (tmpi0 < i0)  i0 = tmpi0;
    if (tmpi1 > i1)  i1 = tmpi1;
    if (tmpj0 < j0)  j0 = tmpj0;
    if (tmpj1 > j1)  j1 = tmpj1;
    if (tmpk0 < k0)  k0 = tmpk0;
    if (tmpk1 > k1)  k1 = tmpk1;
  }
}


void EclipseCarFinGroup::ConvertToEclipseGrid(const EclipseGrid & mother_grid,
                                              EclipseGrid       & grid_out) const
{
  assert(cartesian_refinements_.size() != 0);

  size_t nx, ny, nz;
  size_t i0, i1, j0, j1, k0, k1;
  GetLocation(i0, i1, j0, j1, k0, k1);

  nx = (i1 - i0 + 1) * refinement_x_;
  ny = (j1 - j0 + 1) * refinement_y_;
  nz = (k1 - k0 + 1) * refinement_z_;

  grid_out.Resize(nx, ny, nz);
  for (unsigned int i = 0; i < nx*ny*nz; ++i) {
    grid_out.SetActive(i, false);
  }

  // Make new grid with updated pillars.
  RefineGrid(mother_grid, i0, i1, j0, j1, k0, k1, grid_out);

  // Copy properties from the local grids into the result grid.

  RefinementConstIterator it = RefinementsBegin();

  assert (it != RefinementsEnd()); // Group without refinements is not valid.

  // Add parameters to output grid
  EclipseGrid::ParameterConstIterator p_it = (*it)->GetGrid().ParametersBegin();
  for ( ; p_it != (*it)->GetGrid().ParametersEnd(); ++p_it) {
    grid_out.AddParameter(p_it->first);
  }

  for ( ; it != RefinementsEnd(); ++it) {
    size_t ref_i0, ref_i1, ref_j0, ref_j1, ref_k0, ref_k1;
    (*it)->GetLocation(ref_i0, ref_i1, ref_j0, ref_j1, ref_k0, ref_k1);

    const EclipseGrid& local_grid = (*it)->GetGrid();

    size_t lnx = local_grid.GetNI();
    size_t lny = local_grid.GetNJ();
    size_t lnz = local_grid.GetNK();

    size_t ioff = (ref_i0 - i0) * refinement_x_;
    size_t joff = (ref_j0 - j0) * refinement_y_;
    size_t koff = (ref_k0 - k0) * refinement_z_;

    // ZCORN
    //for (size_t k = 0; k < lnz; ++k) {
    //  for (size_t j = 0; j < lny; ++j) {
    //    for (size_t i = 0; i < lnx; ++i) {
    //      for (int c = 0; c < 2; ++c) {
    //        for (int b = 0; b < 2; ++b) {
    //          for (int a = 0; a < 2; ++a) {
    //            double z = local_grid.GetZCorner(i, j, k, a, b, c);
    //            grid_out.SetZCorner(i + ioff, j + joff, k + koff, a, b, c, z);
    //          }
    //        }
    //      }
    //    }
    //  }
    //}

    // ACTNUM
    for (size_t k = 0; k < lnz; ++k) {
      for (size_t j = 0; j < lny; ++j) {
        for (size_t i = 0; i < lnx; ++i) {
          bool active = local_grid.IsActive(i, j, k);
          grid_out.SetActive(i + ioff, j + joff, k + koff, active);
        }
      }
    }

    EclipseGrid::ParameterConstIterator para_it = local_grid.ParametersBegin();
    for ( ; para_it != local_grid.ParametersEnd(); ++para_it) {
      assert (grid_out.HasParameter(para_it->first));
      const NRLib::Grid<double> & local_parameter = para_it->second;
      NRLib::Grid<double>       & out_parameter   = grid_out.GetParameter(para_it->first);
      for (size_t k = 0; k < lnz; ++k) {
        for (size_t j = 0; j < lny; ++j) {
          for (size_t i = 0; i < lnx; ++i) {
            out_parameter(i + ioff, j + joff, k + koff) = local_parameter(i, j, k);
          }
        }
      }
    }
  }
}


bool EclipseCarFinGroup::IsConnected(const EclipseCarFin& car_fin) const
{
  std::vector<const EclipseCarFin*>::const_iterator it = cartesian_refinements_.begin();
  for (; it != cartesian_refinements_.end(); ++it) {
    if ((*it)->IsConnected(car_fin)) {
      return true;
    }
  }
  return false;
}


void EclipseCarFinGroup::RefineGrid(const EclipseGrid & mother_grid,
                                    size_t              i_from,
                                    size_t              i_to,
                                    size_t              j_from,
                                    size_t              j_to,
                                    size_t              k_from,
                                    size_t              k_to,
                                    EclipseGrid       & grid_out)
{
  // Sets the coord-lines. Works for general local mother_grid refinements.
  size_t mother_nx = i_to - i_from + 1;
  size_t mother_ny = j_to - j_from + 1;

  size_t nx = grid_out.GetNI();
  size_t ny = grid_out.GetNJ();
  size_t nz = grid_out.GetNK();

  for (size_t j = 0; j <= grid_out.GetNJ(); ++j) {
    for (size_t i = 0; i <= grid_out.GetNI(); ++i) {
      size_t i_corner = (mother_nx*i) / nx + i_from;
      size_t i_mod    = (mother_nx*i) % nx;
      size_t j_corner = (mother_ny*j) / ny + j_from;
      size_t j_mod    = (mother_ny*j) % ny;
      NRLib::Line corner1 =
        mother_grid.GetPillar(i_corner, j_corner);
      if ( (0 == i_mod) && (0 == j_mod) ) {
        grid_out.SetPillar(i, j, corner1);
      } else if (0 == j_mod) {
        NRLib::Line corner2 =
          mother_grid.GetPillar(i_corner + 1, j_corner);
        Point top1 = corner1.GetPt1();
        Point bot1 = corner1.GetPt2();
        Point top2 = corner2.GetPt1();
        Point bot2 = corner2.GetPt2();
        double i_fac = (1.0*mother_nx*i)/nx - (mother_nx*i)/nx;
        Point top = i_fac * top2 + (1-i_fac) * top1;
        Point bot = i_fac * bot2 + (1-i_fac) * bot1;
        grid_out.SetPillar(i, j, NRLib::Line(top, bot));
      } else if (0 == i_mod) {
        NRLib::Line corner3 =
          mother_grid.GetPillar(i_corner, j_corner + 1);
        Point top1 = corner1.GetPt1();
        Point bot1 = corner1.GetPt2();
        Point top3 = corner3.GetPt1();
        Point bot3 = corner3.GetPt2();
        double j_fac = (1.0*mother_ny*j)/ny - (mother_ny*j)/ny;
        Point top = j_fac * top3 + (1-j_fac) * top1;
        Point bot = j_fac * bot3 + (1-j_fac) * bot1;
        grid_out.SetPillar(i, j, NRLib::Line(top, bot));
      } else {
        NRLib::Line corner2 =
          mother_grid.GetPillar(i_corner + 1, j_corner);
        NRLib::Line corner3 =
          mother_grid.GetPillar(i_corner, j_corner +1);
        NRLib::Line corner4 =
          mother_grid.GetPillar(i_corner + 1, j_corner + 1);
        Point top1 = corner1.GetPt1();
        Point bot1 = corner1.GetPt2();
        Point top2 = corner2.GetPt1();
        Point bot2 = corner2.GetPt2();
        Point top3 = corner3.GetPt1();
        Point bot3 = corner3.GetPt2();
        Point top4 = corner4.GetPt1();
        Point bot4 = corner4.GetPt2();
        double i_fac = (1.0*mother_nx*i)/nx - (mother_nx*i)/nx;
        double j_fac = (1.0*mother_ny*j)/ny - (mother_ny*j)/ny;
        Point top = j_fac * (i_fac * top4 + (1-i_fac) * top3)
          + (1-j_fac) * (i_fac * top2 + (1-i_fac) * top1);
        Point bot = j_fac * (i_fac * bot4 + (1-i_fac) * bot3)
          + (1-j_fac) * (i_fac * bot2 + (1-i_fac) * bot1);
        grid_out.SetPillar(i, j, NRLib::Line(top, bot));
      }
    }
  }

  // Update the ZCorn values.
  // Only works for grids where each local cell is completely inside a
  // single cell in the mother grid.
  size_t rx = nx/(i_to - i_from + 1);
  size_t ry = ny/(j_to - j_from + 1);
  size_t rz = nz/(k_to - k_from + 1);

  for (size_t k = 0; k < nz; ++k) {
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        size_t mother_i = i_from + i/rx;
        size_t mother_j = j_from + j/ry;
        size_t mother_k = k_from + k/rz;
        size_t mother_nz = k_to - k_from + 1;
        double w_i[2]; double w_j[2]; double w_k[2];
        int ii;
        for (ii = 0; ii < 2; ++ii) {
          w_i[ii] = (1.0*mother_nx*(i+ii))/nx - (mother_nx*(i+ii))/nx;
          w_j[ii] = (1.0*mother_ny*(j+ii))/ny - (mother_ny*(j+ii))/ny;
          w_k[ii] = (1.0*mother_nz*(k+ii))/nz - (mother_nz*(k+ii))/nz;
          if (ii == 1) {
            if (w_i[ii] == 0) w_i[ii] = 1;
            if (w_j[ii] == 0) w_j[ii] = 1;
            if (w_k[ii] == 0) w_k[ii] = 1;
          }
        }
        for (ii = 0; ii < 2; ++ii) {
          for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
              double z =
                  w_k[kk]       * (      w_j[jj] * (      w_i[ii] * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 1, 1, 1)
                                                    + (1-w_i[ii]) * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 0, 1, 1))
                                   + (1-w_j[jj]) * (      w_i[ii] * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 1, 0, 1)
                                                    + (1-w_i[ii]) * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 0, 0, 1)))
                  + (1-w_k[kk]) * (      w_j[jj] * (      w_i[ii] * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 1, 1, 0)
                                                    + (1-w_i[ii]) * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 0, 1, 0))
                                   + (1-w_j[jj]) * (      w_i[ii] * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 1, 0, 0)
                                                    + (1-w_i[ii]) * mother_grid.GetZCorner(mother_i, mother_j, mother_k, 0, 0, 0)));
              grid_out.SetZCorner(i, j, k, ii, jj, kk, z);
            }
          }
        }
      }
    }
  }
}
