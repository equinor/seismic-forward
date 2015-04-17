// $Id: eclipsegeometry.hpp 1187 2013-06-14 13:34:45Z perroe $

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

#ifndef NRLib_GEOMETRY_GRID
#define NRLib_GEOMETRY_GRID

#include <cassert>
#include <iostream>

#include "../geometry/line.hpp"
#include "../geometry/plane.hpp"
#include "../geometry/point.hpp"
#include "../geometry/polygon.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid2d.hpp"
#include "../pointset/pointset.hpp"
#include "../surface/regularsurface.hpp"

namespace NRLib {

class EclipseGeometry {

public:
  // constructors
  EclipseGeometry();
  EclipseGeometry(size_t ni, size_t nj, size_t nk);

  size_t GetNI() const    { return ni_; }
  size_t GetNJ() const    { return nj_; }
  size_t GetNK() const    { return nk_; }
  size_t GetN()  const    { return active_.GetN(); } // = ni_*nj_*nk_

  inline void   GetIJK(size_t index, size_t& i, size_t& j, size_t& k) const;
  inline size_t GetIndex(size_t i, size_t j, size_t k) const;

  void Resize(size_t ni, size_t nj, size_t nk);

  inline bool IsActive(size_t i, size_t j, size_t k) const;
  inline bool IsActive(size_t index) const;
  inline bool IsPillarActive(size_t i, size_t j) const;

  const Line& GetPillar(size_t i, size_t j) const   { return pillars_(i, j); }

  /// Returns the z value for corner (a,b,c) in cell (i,j,k).
  inline double GetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const;

  /// Sets given input line (new_pillar) to local variable, if valid input index
  inline void SetCoordLine(size_t i, size_t j, const Line& pillar);

  /// Sets given input value (double z_in) to local variable, if valid input index
  inline void SetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c, double z);

  inline void SetActive(size_t i, size_t j, size_t k, bool is_active);

  inline void SetActive(size_t index, bool is_active);

  /// Returns the corner point for corner (a,b,c) in cell (i,j,k).
  inline Point FindCornerPoint(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const;

  /// Returns grid value in cell no (i,j,k), average of cornervalues
  /// Uses traditional ECLIPSE geometry representation.
  Point FindCellCenterPoint(size_t i, size_t j, size_t k) const;

  Point FindCellCenterPointTopBot(size_t i, size_t j, size_t k, bool top) const;

  /// Finds the cell (i,j,k) the given point (x,y,z) is inside.
  /// Uses new ECLIPSE geometry representation.
  /// Returns true if such a cell is found.
  bool FindIndex(double x, double y, double z, size_t& i, size_t& j, size_t& k) const;

  /// Finds the point on the pillar with the given z-value
  Point FindPointAtPillar(size_t i, size_t j, double z) const;

  /// Finds the point on the pillar with given z-valu. If that point is outside the grid, the closest pillar point in the grid is chosen
  Point FindPointAtPillarInsideGrid(size_t i, size_t j, double z, bool & found) const;

  /// Finds the point in the cell (i, j, k) with local coordinates (u, v, w).
  /// u, v and w must be between 0 and 1.
  Point FindPointInCell(size_t i, size_t j, size_t k, double u, double v, double w);

  /// Returns nk if all cells are inactive.
  size_t FindTopCell(size_t i, size_t j) const;
  size_t FindBottomCell(size_t i, size_t j) const;
  size_t FindTopLayer() const;
  size_t FindBottomLayer() const;

  bool IsColumnActive(size_t i, size_t j) const
  { return FindTopCell(i, j) < GetNK(); }

  /// Finds minimum area enclosing rectangle that contains all points in polymin_ and polymax_.
  /// This rectangle has a corner in (x0,y0), sides of length lx and ly, and angle is the angle between the side with length lx and the x-axis
  void FindEnclosingVolume(double& x0,
                           double & y0,
                           double& lx,
                           double& ly,
                           double& angle) const;

  ///Stores z-values in layer k for a rectangle with a corner in x0,y0 and step lengths dx and dy, angle indicates rotated angle in the xy-plane
  ///\\param lower_or_upper 0 for upper, 1 for lower
  ///\\param bilinear_else_triangles true for calulating z-coordinates inside corners by bilinear interpolation, false for calculating by intersection of plane through triangles
  void FindLayerSurface(NRLib::Grid2D<double> &z_surface,
                        size_t k,
                        int lower_or_upper,
                        double dx,
                        double dy,
                        double x0,
                        double y0,
                        double angle,
                        bool bilinear_else_triangles) const;

  void FindLayerSurfaceCornerpoint(NRLib::Grid2D<double> &z_surface,
                                   size_t k,
                                   int lower_or_upper,
                                   double dx,
                                   double dy,
                                   double x0,
                                   double y0,
                                   double angle,
                                   bool bilinear_else_triangles) const;

  void FindTopAndBotValuesOfGrid(std::vector<NRLib::Point>& top_points,
                                 std::vector<NRLib::Point>& bot_points) const;

  void FindMinAndMaxZValueAndSetPolygons();

  void WriteSpecGrid(std::ofstream& out_file) const;
  void WriteCoord   (std::ofstream& out_file) const;
  void WriteZCorn   (std::ofstream& out_file) const;
  void WriteActNum  (std::ofstream& out_file) const;
  void WriteGeometry(std::ofstream& out_file) const;

  void ReadSpecGrid(std::ifstream& in_file);
  void ReadZCorn   (std::ifstream& in_file);
  void ReadCoord   (std::ifstream& in_file);
  void ReadActNum  (std::ifstream& in_file);

  double GetDZ(size_t i, size_t j, size_t k) const;

private:
  // ----------------- PRIVATE FUNCTIONS ---------------------------
    /// Finds the z-value at the pillar on top of the grid
  double FindZTopAtPillar(size_t i, size_t j, bool & found) const;

  /// Finds the z-value at the pillar on bottom of the grid
  double FindZBotAtPillar(size_t i, size_t j, bool & found) const;

  /// Calculates the index in zcoord_-vector for cell (i, j, k), corner (a, b, c).
  inline size_t FindZCornerIndex(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const;

  /// Find column (i and j coordinate for cell), containing the point (x, y, z).
  bool GetCellIndex_ij(double x, double y, double z, size_t& i, size_t& j) const;

  /// Finds the point on the top/bottom surface
  /// \param lower_or_upper 0 for upper, 1 for lower.
  Point FindPointCellSurface(size_t i, size_t j, size_t k, int lower_or_upper,
                             double u, double v) const;

  /// Finds the local u and v coordinates when i, j, x, y and z are known.
  void FindUVCoordinates(double x, double y, double z, size_t i, size_t j,
                         double& u, double& v) const;

  /// Finds polygon around active pillars at the given level.
  Polygon FindPolygonAroundActivePillars(double z_in);

  ///Search method used in FindPolygonAroundActivePillars
  size_t SearchUp(size_t i,
                  size_t j,
                  std::vector<Point> &p,
                  NRLib::Grid2D<bool> &checked,
                  double z_in,
                  size_t start_i,
                  size_t start_j);

  ///Search method used in FindPolygonAroundActivePillars
  size_t SearchDown(size_t i,
                    size_t j,
                    std::vector<Point> &p,
                    NRLib::Grid2D<bool> &checked,
                    double z_in,
                    size_t start_i,
                    size_t start_j);

  ///Search method used in FindPolygonAroundActivePillars
  size_t SearchRight(size_t i,
                     size_t j,
                     std::vector<Point> &p,
                     NRLib::Grid2D<bool> &checked,
                     double z_in,
                     size_t start_i,
                     size_t start_j);

  ///Search method used in FindPolygonAroundActivePillars
  size_t SearchLeft(size_t i,
                    size_t j,
                    std::vector<Point> &p,
                    NRLib::Grid2D<bool> &checked,
                    double z_in,
                    size_t start_i,
                    size_t start_j);

  /// Initialize the active_pillars_ grid.
  void InitializeActivePillars();

  ///Function used by FindLayerSurface to fill in values to z_surface in the area inside the (NB) four corners (listed clockwise)
  void TriangularFillInZValuesInArea(NRLib::Grid2D<double> &z_surface,
                                     NRLib::Grid2D<int> &is_set,
                                     double x0,
                                     double y0,
                                     const std::vector<NRLib::Point>& corners,
                                     double dx,
                                     double dy) const;

  ///Function used by FindLayerSurface to fill in values to z_surface in the area inside the (NB) four corners (listed clockwise)
  void BilinearFillInZValuesInArea(NRLib::Grid2D<double> &z_surface,
                                   NRLib::Grid2D<int> &is_set,
                                   double x0,
                                   double y0,
                                   const std::vector<NRLib::Point>& corners,
                                   double dx,
                                   double dy) const;

  ///Function used by FindLayerSurface to fill in (average values of neighbour elements) to z_surface where is_set==false
  void FillInZValuesByAveraging(NRLib::Grid2D<double> &z_surface,
                                NRLib::Grid2D<int> &is_set) const;


  // --------------------- MEMBER VARIABLES -------------------------

  /// number of cells: ni*nj*nk.
  size_t ni_, nj_, nk_;

  /// Polygon around active pillars on maximum z coordinates.
  Polygon polymax_;

  /// Polygon around active pillars on minimum z coordinates.
  Polygon polymin_;

  /// Grid pillars. Size: (nx+1)*(ny+1). Each pillar is defined by a top and a bottom point.
  Grid2D<Line> pillars_;

  /// Indicator if given pillar is active or not. A pillar is inactive if all adjacent columns
  /// are inactive. Size: (nx + 1)*(ny + 1).
  Grid2D<int> active_pillars_;

  /// Z-corner values for the individual cells. Same order as in the eclipse grid files.
  /// Storage order: a, i, b, j, c, k. (Cell indexes (i, j, k, a, b, c)).
  std::vector<double> z_corners_;

  /// Indicator if given cell is active or not. True if cell is active.
  Grid<bool> active_;
};   //  class EclipseGeometry


// ===================== INLINE MEMBER FUNCTIONS ======================


void EclipseGeometry::GetIJK(size_t index, size_t& i, size_t& j, size_t& k) const
{
  assert(index < GetN());

  i = index % ni_;
  j = (index-i)/ni_ % nj_;
  k = (index - j*ni_ - i)/ni_/nj_;
}


size_t EclipseGeometry::GetIndex(size_t i, size_t j, size_t k) const
{
  assert(i < GetNI() && j < GetNI() && k < GetNK());

  return i + j*GetNI() + k*GetNI()*GetNJ();
}


bool EclipseGeometry::IsActive(size_t i, size_t j, size_t k) const
{
  return active_(i, j, k);
}


bool EclipseGeometry::IsActive(size_t index) const
{
  return active_(index);
}


bool EclipseGeometry::IsPillarActive(size_t i, size_t j) const
{
  return (active_pillars_(i, j) != 0);
}


double EclipseGeometry::GetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const
{
  size_t index = FindZCornerIndex(i, j, k, a, b, c);
  return z_corners_[index];
}


void EclipseGeometry::SetCoordLine(size_t i, size_t j, const Line& pillar)
{
  pillars_(i, j) = pillar;
}


void EclipseGeometry::SetZCorner(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c, double z)
{
  size_t index = FindZCornerIndex(i, j, k, a, b, c);
  z_corners_[index] = z;
}


void EclipseGeometry::SetActive(size_t i, size_t j, size_t k, bool is_active)
{
  active_(i, j, k) = is_active;
}


void EclipseGeometry::SetActive(size_t index, bool is_active)
{
  active_(index) = is_active;
}


Point EclipseGeometry::FindCornerPoint(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const
{
  double z = GetZCorner(i, j, k, a, b, c);
  size_t pillar_i = i + a;
  size_t pillar_j = j + b;
  return FindPointAtPillar(pillar_i, pillar_j, z);
}


size_t EclipseGeometry::FindZCornerIndex(size_t i, size_t j, size_t k, size_t a, size_t b, size_t c) const
{
  assert(i < ni_);
  assert(j < nj_);
  assert(k < nk_);
  assert(a == 0 || a == 1);
  assert(b == 0 || b == 1);
  assert(c == 0 || c == 1);

  return (2*i + a + 2*ni_*(2*j + b) + 4*ni_*nj_*(2*k + c));
}

}  //  end of namespace NRLib

#endif
