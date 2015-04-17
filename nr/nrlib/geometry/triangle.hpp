// $Id: triangle.hpp 882 2011-09-23 13:10:16Z perroe $

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

#ifndef NRLIB_GEOMETRY_TRIANGLE_HPP
#define NRLIB_GEOMETRY_TRIANGLE_HPP

#include "point.hpp"

namespace NRLib {

class Line;

class Triangle {
public:
  Triangle();
  Triangle(const Point& p1, const Point& p2, const Point& p3);

  inline void SetCornerPoints(const Point& p1, const Point& p2, const Point& p3);

  bool FindIntersection(const Line& line, Point& intersec_pt, bool include_edge) const;

  // Find point in triangle nearest line. Returns the distance to the line.
  double FindNearestPoint(const Line& line, Point& nearest_pt) const;

private:
  Point p1_;
  Point p2_;
  Point p3_;
};


inline void Triangle::SetCornerPoints(const NRLib::Point & p1,
                                      const NRLib::Point & p2,
                                      const NRLib::Point & p3)
{
  p1_ = p1;
  p2_ = p2;
  p3_ = p3;
}

} // namespace NRLib

#endif // NRLIB_GEOMETRY_TRIANGLE_HPP
