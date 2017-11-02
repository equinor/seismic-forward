// $Id: plane.cpp 1195 2013-09-04 13:45:58Z perroe $

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

#include "plane.hpp"
#include "line.hpp"
#include "point.hpp"
#include <cmath>
#include <limits>

using namespace NRLib;

const double Plane::min_ = 0.0000001;  // def min-value, might be changed.

Plane::Plane()
{
}


Plane::Plane(const Point& pt1_in, const Point& pt2_in, const Point& pt3_in)
{
  pt_     = pt1_in;             // choose p1 as def.point
  Point a = pt1_in - pt2_in;    // line from p1 to p2
  Point b = pt1_in - pt3_in;    // line from p1 to p3
  Point cross_ab = Cross(a, b);

  assert(cross_ab.Norm() > 0); //Points are not on a line

  nvec1_ = cross_ab.GetNormalized();  // normal unit vec for plane
}


Plane::Plane(const Point& pt_in, const Point& nvec_in)
{
  pt_    = pt_in;
  nvec1_ = nvec_in.GetNormalized(); // scale input point/vector to 1.0
}


Point Plane::GetNormalVector() const
{
    return nvec1_;
}


Point Plane::GetRefPoint() const
{
   return pt_;
}


double Plane::GetMin() const
{
   return min_;
}


Point Plane::FindIntersection(const Line& line_in) const
{
  // Check whether line_in and plane are parallel -> no intersection
  assert(!IsParallel(line_in, 0.0));

  Point pt1 = line_in.GetPt1();
  Point pt2 = line_in.GetPt2();
  double t  = nvec1_.Dot(pt_ - pt1) / nvec1_.Dot(pt2 - pt1);

  Point intersec = pt1 + t*(pt2 - pt1);
  return intersec;
}


bool Plane::FindIntersection(const Line& line_in, Point& intersection, double min_in) const
{
  if (IsParallel(line_in, min_in)) {
    return false;
  }

  Point pt1 = line_in.GetPt1();
  Point pt2 = line_in.GetPt2();
  double t  = nvec1_.Dot(pt_ - pt1) / nvec1_.Dot(pt2 - pt1);

  intersection = pt1 + t*(pt2 - pt1);

  if ((line_in.IsEndPt1() && t < 0) || (line_in.IsEndPt2() && t > 1))
    return false;

  return true;
}


bool Plane::IsParallel(const Line& line_in, double min_in) const
{
  Point pt_12 = line_in.GetPt2() - line_in.GetPt1();
  double dot  = nvec1_.Dot( pt_12 );
  if (std::fabs(dot) <= min_in)
    return true;
  else
    return false;
}


Point Plane::FindProjection(const Point& p_in) const
{
  double dot = nvec1_.Dot(pt_ - p_in);
  return (p_in + dot*nvec1_);
}


double Plane::FindDistance(const Point& p_in) const
{
  return std::abs(nvec1_.Dot(p_in - pt_));
}


bool Plane::IsPointOnPlane(const Point& p_in, double min_in) const
{
  if (FindDistance(p_in) <= min_in)
    return true;
  return false;
}


Line Plane::FindIntersection(const Plane& plane_in) const
{
  Point n1    = nvec1_;                     // have been scaled to 1
  Point n2    = plane_in.GetNormalVector(); // vector have been scaled to 1
  Point cross = Cross(n1, n2);              // cross-vector n1 and n2

  if (cross.Norm() == 0) // Parallell planes.
    return Line();       // Line with no length.

  double d1     = -n1.Dot(pt_);
  double d2     = -n2.Dot(plane_in.GetRefPoint());
  double dot    = n1.Dot(n2);
  double one_minus_dotdot = 1.0 - dot*dot;
  double const1           = (d2*dot - d1) / one_minus_dotdot;
  double const2           = (d1*dot - d2) / one_minus_dotdot;
  Point inter_pt1 = const1*n1 + const2*n2;  // line when t=0

  // define the line of intersection
  Line inter_line = Line::InfiniteLine(inter_pt1, cross);
  return inter_line;
}


bool Plane::IsOnPositiveSide(const NRLib::Point& p_in)
{
  if (nvec1_.Dot(p_in - pt_) > 0)
    return true;
  return false;
}

