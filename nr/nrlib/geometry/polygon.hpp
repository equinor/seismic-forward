// $Id: polygon.hpp 1141 2013-02-07 10:21:39Z georgsen $

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

#ifndef NRLIB_GEOMETRY_POLYGON
#define NRLIB_GEOMETRY_POLYGON

#include <cassert>
#include <vector>

#include "../geometry/point.hpp"

namespace NRLib {

class Polygon  {

public:
  // constructors:
  Polygon();
  Polygon(size_t size_in);

  // Add input-point to end of polygon, and increase size of polygon
  void AddPoint(const Point& pt_in);
  // Returns size/number of points of given polygon
  size_t  GetSize() const;

  // Gets the polygon points
  inline std::vector<Point>& GetPoints(std::vector<Point>& pt_out, bool closed = false);
  inline const std::vector<Point>& GetPoints(std::vector<Point>& pt_out, bool closed = false) const;

  inline Point& operator()(size_t index);
  inline const Point& operator()(size_t index) const;

  Point& operator[](size_t index)                  { return (*this)(index); }
  const Point& operator[](size_t index) const      { return (*this)(index); }

  // Check if given pt_in is inside polygon, return true implies yes
  bool IsInsidePolygonXY(const Point& pt_in) const { return IsInsidePolygonXY(pt_in.x, pt_in.y); }

  // Check if given pt_in is inside polygon, return true implies yes
  bool IsInsidePolygonXY(double x, double y) const;

  /// Check if the polygon points are listed counterclockwise (neglecting z-component, 2D algorithm)
  bool IsCounterclockwise();

  /// Find the x and y limits of the polygon
  void XYExtent(double& xmin, double& xmax, double& ymin, double& ymax);

  /// Check if the polygon is convex (neglecting z-component, 2D algorithm)
  bool IsConvex();

  /// Deletes points if necessary to make polygon convex (neglecting z-component, 2D algorithm)
  //Note: Assuming the number of points that needs to be deleted is small, otherwise: GetConvexHull in PointSetSurface might be more effective.
  void MakeConvex();

  ///Finds minimum area enclosing rectangle using Rotating Calipers (neglecting z-component, 2D algorithm)
  void MinEnclosingRectangle(double &x,
                             double &y,
                             double &length1,
                             double &length2,
                             double &angle);


private:
  std::vector<Point> polygon_points_;

};   // end of class Polygon


inline Point& Polygon::operator()(size_t index)
{
  assert(index < GetSize());

  return(polygon_points_[index]);
}


inline const Point& Polygon::operator()(size_t index)  const
{
  assert(index < GetSize());

  return(polygon_points_[index]);
}

inline std::vector<Point>& Polygon::GetPoints(std::vector<Point>& pt_out, bool closed)
{
  pt_out = polygon_points_;
  if(closed && GetSize() > 0)
    pt_out.push_back(polygon_points_[0]);

  return(pt_out);
}

inline const std::vector<Point>& Polygon::GetPoints(std::vector<Point>& pt_out, bool closed) const
{
  pt_out = polygon_points_;
  if(closed && GetSize() > 0)
    pt_out.push_back(polygon_points_[0]);

  return(pt_out);
}


}   // end of namespace NRLib

#endif  // NRLib_POLYGON

