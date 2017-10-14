// $Id: triangle.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "triangle.hpp"

#include "line.hpp"
#include "point.hpp"

namespace NRLib {

Triangle::Triangle()
{}

Triangle::Triangle(const Point& p1, const Point& p2, const Point& p3)
  : p1_(p1),
    p2_(p2),
    p3_(p3)
{}

bool Triangle::FindIntersection(const Line& line, Point& intersec_pt, bool include_edge) const
{
  // Based on algorithm by softsurfer:
  // http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm

  Point u = p2_ - p1_;
  Point v = p3_ - p1_;

  Point n = NRLib::Cross(u, v);

  if (n.x == 0.0 && n.y == 0.0 && n.z == 0.0) {
    // Degenerated triangle. Line or point.
    return false;
  }

  Point line_dir = line.GetPt2() - line.GetPt1();
  Point w0 = line.GetPt1() - p1_;

  double a = -NRLib::Dot(n, w0);
  double b = NRLib::Dot(n, line_dir);

  if (b == 0.0) {
    // Line paralell to surface.
    // We treat this as no intersection, even if line goes trough triangle.
    return false;
  }

  // Intersection point, plane with ray.
  double r = a / b;
  if ( (line.IsEndPt1() && r < 0.0) ||
       (line.IsEndPt2() && r > 1.0) ) {
    // Intersection outside segment/ray.
    return false;
  }

  intersec_pt = line.GetPt1() + r * line_dir;

  // Check if intersec_pt is inside triangle:
  double uu = NRLib::Dot(u, u);
  double uv = NRLib::Dot(u, v);
  double vv = NRLib::Dot(v, v);

  NRLib::Point w = intersec_pt - p1_;

  double wu = NRLib::Dot(w, u);
  double wv = NRLib::Dot(w, v);

  double D = uv*uv - uu*vv;

  double s = (uv*wv - vv*wu) / D;
  if (include_edge) {
    if (s < 0.0 || s > 1.0) {
      return false;
    }
  }
  else {
    if (s <= 0.0 || s >= 1.0) {
      return false;
    }
  }

  double t = (uv*wu - uu*wv) / D;
  if (include_edge) {
    if (t < 0.0 || (s + t) > 1.0) {
      return false;
    }
  }
  else {
    if (t <= 0.0 || (s + t) >= 1.0) {
      return false;
    }
  }

  return true;
}


double Triangle::FindNearestPoint(const Line& line, Point& nearest_pt) const
{
  if (FindIntersection(line, nearest_pt, true))
    return 0.0;

  NRLib::Line edge = NRLib::Line::Segment(p1_, p2_);
  NRLib::Point dummy_point; // This point is not used.
  double min_dist = edge.FindDistance(line, nearest_pt, dummy_point);

  NRLib::Point pt;
  edge.SetPt(p2_, p3_);
  double dist = edge.FindDistance(line, pt, dummy_point);
  if (dist < min_dist) {
    min_dist = dist;
    nearest_pt = pt;
  }

  edge.SetPt(p1_, p3_);
  dist = edge.FindDistance(line, pt, dummy_point);
  if (dist < min_dist) {
    min_dist = dist;
    nearest_pt = pt;
  }

  return min_dist;
}


} // namespace NRLib
