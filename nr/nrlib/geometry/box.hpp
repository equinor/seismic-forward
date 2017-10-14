// $Id: box.hpp 1068 2012-09-18 11:21:53Z perroe $

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

#ifndef NRLIB_GEOMETRY_BOX_HPP
#define NRLIB_GEOMETRY_BOX_HPP

#include "point.hpp"

#include <vector>

namespace NRLib {

class Line;

class Box {
public:
  Box();

  Box(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z); // Box with angle 0.
  Box(const Point& refpt, double lengthx, double lengthy, double lengthz, double angle);

  double GetXRef() const { return refpt_.x; }
  double GetYRef() const { return refpt_.y; }
  double GetZRef() const { return refpt_.z; }
  double GetLX()   const { return lengthx_; }
  double GetLY()   const { return lengthy_; }
  double GetLZ()   const { return lengthz_; }
  double GetAngle() const { return angle_; }

  /// \brief Check if point is inside box.
  bool IsInside(const NRLib::Point& p) const
  { Point p_loc = GlobalToLocal(p);
    return (p_loc.x >= 0.0 && p_loc.x <= lengthx_ &&
            p_loc.y >= 0.0 && p_loc.y <= lengthy_ &&
            p_loc.z >= 0.0 && p_loc.z <= lengthz_ ); }

  /// \brief Find intersections between line and box sides.
  /// \returns 0-2 intersections.
  std::vector<NRLib::Point> FindIntersections(const NRLib::Line& line) const;

  /// \brief Find distance to box
  double FindDistance(const Point& p) const;

  Point GetCorner(int i, int j, int k) const;
  Point GetCorner(bool i, bool j, bool k) const; //Don't use this, use the int version above. This is kept until we know that it is not used any place.

  inline void Swap(Box& other) {
    refpt_.Swap(other.refpt_);
    std::swap(lengthx_, other.lengthx_);
    std::swap(lengthy_, other.lengthy_);
    std::swap(lengthz_, other.lengthz_);
    std::swap(angle_, other.angle_);
  }

private:
  Point GlobalToLocal(const NRLib::Point & p) const;
  Point LocalToGlobal(const NRLib::Point & p) const;

  Point refpt_;
  double lengthx_;
  double lengthy_;
  double lengthz_;
  double angle_;    ///< in radians
};

}

#endif // NRLIB_GEOMETRY_BOX_HPP
