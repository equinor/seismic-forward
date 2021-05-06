// $Id: box.cpp 1068 2012-09-18 11:21:53Z perroe $

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

#include "box.hpp"

#include "line.hpp"
#include "plane.hpp"

#include <cassert>
#include <vector>

namespace NRLib {

Box::Box()
  : refpt_(0.0, 0.0, 0.0),
    lengthx_(0.0),
    lengthy_(0.0),
    lengthz_(0.0),
    angle_(0.0)
{}


Box::Box(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z) // for angle =0
  : refpt_(min_x, min_y, min_z),
    lengthx_(max_x - min_x),
    lengthy_(max_y - min_y),
    lengthz_(max_z - min_z),
    angle_(0.0)
{
  assert(min_x <= max_x && min_y <= max_y && min_z <= max_z);
}

Box::Box(const Point& refpt, double lengthx, double lengthy, double lengthz, double angle)
  : refpt_(refpt.x, refpt.y, refpt.z),
    lengthx_(lengthx),
    lengthy_(lengthy),
    lengthz_(lengthz),
    angle_(angle)
{
}

std::vector<Point> Box::FindIntersections(const NRLib::Line& line) const
{
  std::vector<Point> intersections;
  Point p1 = GlobalToLocal(line.GetPt1());
  Point p2 = GlobalToLocal(line.GetPt2());
  Line line_loc(p1, p2, line.IsEndPt1(), line.IsEndPt2());

  Plane side(Point(0, 0, 0), Point(1, 0, 0));
  Point p;
  if (side.FindIntersection(line_loc, p)) {
    if (p.y >= 0.0 && p.y <= lengthy_ && p.z >= 0.0 && p.z <= lengthz_) {
      intersections.push_back(p);
    }
  }

  side = Plane(Point(lengthx_, 0, 0), Point(1, 0, 0));
  if (side.FindIntersection(line_loc, p)) {
    if (p.y >= 0.0 && p.y <= lengthy_ && p.z >= 0.0 && p.z <= lengthz_) {
      intersections.push_back(p);
    }
  }

  side = Plane(Point(0, 0, 0), Point(0, 1, 0));
  if (side.FindIntersection(line_loc, p)) {
    if (p.x >= 0.0 && p.x <= lengthx_ && p.z >= 0.0 && p.z <= lengthz_) {
      intersections.push_back(p);
    }
  }

  side = Plane(Point(0, lengthy_, 0), Point(0, 1, 0));
  if (side.FindIntersection(line_loc, p)) {
    if (p.x >= 0.0 && p.x <= lengthx_ && p.z >= 0.0 && p.z <= lengthz_) {
      intersections.push_back(p);
    }
  }

  side = Plane(Point(0, 0, 0), Point(0, 0, 1));
  if (side.FindIntersection(line_loc, p)) {
    if (p.x >= 0.0 && p.x <= lengthx_ && p.y >= 0.0 && p.y <= lengthy_) {
      intersections.push_back(p);
    }
  }

  side = Plane(Point(0, 0, lengthz_), Point(0, 0, 1));
  if (side.FindIntersection(line_loc, p)) {
    if (p.x >= 0.0 && p.x <= lengthx_ && p.y >= 0.0 && p.y <= lengthy_) {
      intersections.push_back(p);
    }
  }

  assert(intersections.size() <= 2);

  // Prevent numerical instabilities:
  for (size_t i = 0; i < intersections.size(); ++i) {
    if (intersections[i].x < 0.0) {
      intersections[i].x = 0.0;
    }
    if (intersections[i].y < 0.0) {
      intersections[i].y = 0.0;
    }
    if (intersections[i].z < 0.0) {
      intersections[i].z = 0.0;
    }
    if (intersections[i].x > lengthx_) {
      intersections[i].x = lengthx_;
    }
    if (intersections[i].y > lengthy_) {
      intersections[i].y = lengthy_;
    }
    if (intersections[i].z > lengthz_) {
      intersections[i].z = lengthz_;
    }
    Point global_pt = LocalToGlobal(intersections[i]);
    intersections[i] = global_pt;
  }

  return intersections;
}


double Box::FindDistance(const Point& p) const
{
  Point p_loc = GlobalToLocal(p);
  if (IsInside(p)) {
    double distance_to_side[6];
    distance_to_side[0] = p_loc.z;
    distance_to_side[1] = lengthz_ - p_loc.z;
    distance_to_side[2] = p_loc.y;
    distance_to_side[3] = lengthy_ - p_loc.y;
    distance_to_side[4] = p_loc.x;
    distance_to_side[5] = lengthx_ - p_loc.x;

    double min_dist = distance_to_side[0];
    for (size_t i = 1; i < 6; ++i) {
      if (min_dist > distance_to_side[i])
        min_dist = distance_to_side[i];
    }
    return min_dist;
  }
  else if (p_loc.x < 0.0) {
    if (p_loc.y < 0.0) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(0.0, 0.0, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(0.0, 0.0, lengthz_));
      else
        return p_loc.GetDistance(Point(0.0, 0.0, p_loc.z));
    }
    else if (p_loc.y > lengthy_) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(0.0, lengthy_, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(0.0, lengthy_, lengthz_));
      else
        return p_loc.GetDistance(Point(0.0, lengthy_, p_loc.z));
    }
    else {
      // 0.0 < p.y < lengthy_
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(0.0, p_loc.y, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(0.0, p_loc.y, lengthz_));
      else
        return (- p_loc.x);
    }
  }
  else if (p_loc.x > lengthx_) {
    if (p_loc.y < 0.0) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(lengthx_, 0.0, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(lengthx_, 0.0, lengthz_));
      else
        return p_loc.GetDistance(Point(lengthx_, 0.0, p_loc.z));
    }
    else if (p_loc.y > lengthy_) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(lengthx_, lengthy_, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(lengthx_, lengthy_, lengthz_));
      else
        return p_loc.GetDistance(Point(lengthx_, lengthy_, p_loc.z));
    }
    else {
      // 0.0 < p.y < lengthy_
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(lengthx_, p_loc.y, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(lengthx_, p_loc.y, lengthz_));
      else
        return p_loc.x - lengthx_;
    }
  }
  else {
    // 0.0 < p.x < lengthx_
    if (p_loc.y < 0.0) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(p_loc.x, 0.0, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(p_loc.x, 0.0, lengthz_));
      else
        return (- p_loc.y);
    }
    else if (p_loc.y > lengthy_) {
      if (p_loc.z < 0.0)
        return p_loc.GetDistance(Point(p_loc.x, lengthy_, 0.0));
      else if (p_loc.z > lengthz_)
        return p_loc.GetDistance(Point(p_loc.x, lengthy_, lengthz_));
      else
        return p_loc.y - lengthy_;
    }
    else {
      // 0.0 < p.y < lengthy_
      if (p_loc.z < 0.0)
        return (- p_loc.z);
      else if (p_loc.z > lengthz_)
        return p_loc.z - lengthz_;
      else {
        // p is inside box, handled above.
        assert(0);
        return -1;
      }
    }
  }
}

Point Box::GlobalToLocal(const Point &global_pt) const
{
  Point local_pt;
  double x_rel = global_pt.x - refpt_.x;
  double y_rel = global_pt.y - refpt_.y;
  double z_rel = global_pt.z - refpt_.z;

  local_pt.x = x_rel * cos(angle_) + y_rel * sin(angle_);
  local_pt.y = - x_rel * sin(angle_) + y_rel * cos(angle_);
  local_pt.z = z_rel;
  return local_pt;
}

Point Box::LocalToGlobal(const Point &local_pt) const
{
  Point global_pt;
  global_pt.x = local_pt.x * cos(angle_) - local_pt.y * sin(angle_) + refpt_.x;
  global_pt.y = local_pt.x * sin(angle_) + local_pt.y * cos(angle_) + refpt_.y;
  global_pt.z = local_pt.z + refpt_.z;
  return global_pt;

}

Point Box::GetCorner(int i, int j, int k) const
{
  assert (i == 0 || i == 1);
  assert (j == 0 || j == 1);
  assert (k == 0 || k == 1);
  Point loc_pt;
  if(i == 0)
    loc_pt.x = 0.0;
  else // i = 1
    loc_pt.x = lengthx_;
  if(j == 0)
    loc_pt.y = 0.0;
  else // j = 1
    loc_pt.y = lengthy_;
  if(k == 0)
    loc_pt.z = 0.0;
  else // k = 1
    loc_pt.z = lengthz_;

  Point pt = LocalToGlobal(loc_pt);
  return(pt);
}

Point Box::GetCorner(bool i, bool j, bool k) const
{
  return (GetCorner(static_cast<int> (i), static_cast<int> (j), static_cast<int> (k)));
}

} // namespace NRLib
