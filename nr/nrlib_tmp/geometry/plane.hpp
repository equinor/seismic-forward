// $Id: plane.hpp 882 2011-09-23 13:10:16Z perroe $

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

#ifndef NRLIB_GEOMETRY_PLANE
#define NRLIB_GEOMETRY_PLANE

#include <iostream>

#include "line.hpp"
#include "point.hpp"

namespace NRLib {

class Plane {

public:

  // Constructors:
  Plane();
  Plane(const Point& pt1_in, const Point& pt2_in, const Point& pt3_in);
  Plane(const Point& pt_in, const Point& nvec_in);

  Point  GetRefPoint() const;       // returns pt_
  Point  GetNormalVector() const;   // returns nvec1_
  double GetMin() const;  // returns min_

  /// Calc intersection point between line_in and def.plane
  /// \note Caller should make sure that line_in is not parallell to plane.
  Point  FindIntersection(const Line& line_in) const;

  /// Find intersection between plane and line.
  /// Works also for line segments and rays.
  /// \return true if intersection is found, else false.
  bool   FindIntersection(const Line& line_in, Point& intersection, double min_in = 0.0) const;

  // True if parallel directions for line_in and def.plane
  bool   IsParallel(const Line& line_in, double min_in = 0.0) const;

  // Calc projection point on def.plane from p_in
  Point  FindProjection(const Point& p_in) const;

  // Distance from p_in to projection point on plane
  double FindDistance(const Point& p_in) const;

  // True if shortest distance less than min_in value
  bool   IsPointOnPlane(const Point& p_in, double min_in = min_) const;

  // Calc line of intersection between def.plane and plane_in
  // OBS: Does not work if one of the planes are horizontal.
  Line   FindIntersection(const Plane& plane_in) const;

  /// True if point is on positive side of plane. (Direction of normal vector)
  bool   IsOnPositiveSide(const Point& p_in);

private:
  Point pt_;                 // point on plane
  Point nvec1_;              // unit normal vector (scaled to 1). From input, or calc from p1, p2, pt3)
  static const double min_;  // used instead of zero for some comparision-functions
};  // end of class Plane

}  // end of using namespace NRLib

#endif
