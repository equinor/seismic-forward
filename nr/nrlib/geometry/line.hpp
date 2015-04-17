// $Id: line.hpp 1036 2012-06-25 14:21:49Z anner $

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

#ifndef NRLIB_GEOMETRY_LINE
#define NRLIB_GEOMETRY_LINE

#include <iostream>
#include "point.hpp"

namespace NRLib {


/// Class for handling of lines, rays (lines defined by a starting point and a
/// direction), and line segments (lines with a start and an end point).
class Line {

public:
  /// Default constructor.
  /// Constructs a line segment with zero length located at origo.
  Line();

  /// Generic constructor for all kinds of lines.
  /// In most cases the more specialized "constructors" below should be used.
  Line(const Point& p1_in, const Point &p2_in, bool end1_in = true, bool end2_in = true);

  /// Constructs a line (with infinite length).
  static Line InfiniteLine(const Point& start, const Point& dir)
  { return Line(start, start + dir, false, false); }

  static Line InfiniteLine(const Point& start, const double &angle)
  { Point end = start; end.x+=10*cos(angle); end.y+=10*sin(angle); return Line(start, end, false, false); }

  /// Constructs a ray with the given start point and direction.
  static Line Ray(const Point& start, const Point& dir)
  { return Line(start, start + dir, true, false); }

  /// Constructs a line segment from start to end;
  static Line Segment(const Point& start, const Point& end)
  { return Line(start, end, true, true); }

  // Return p1_ or p2_
  Point& GetPt1()
  { return p1_; }

  const Point& GetPt1() const
  { return p1_; }

  Point& GetPt2()
  { return p2_; }

  const Point& GetPt2() const
  { return p2_; }

  /// Returns the normalized direction vector, or 0 for a line segment with no length.
  Point GetDir() const
  { return (p2_ - p1_).GetNormalized(); }

  bool IsEndPt1() const
  { return end1_; }

  bool IsEndPt2() const
  { return end2_; }

  /// Set pt1=p1_, return value of end1_, corresp. for pt2 :
  bool GetPt1(Point& pt1) const;
  bool GetPt2(Point& pt2) const;

  void SetPt(const Point& p1_in, const Point& p2_in, bool end1_in = true, bool end2_in = true);

  /// Get length of line segment, returns infinite for lines and rays. :
  double GetLength() const;

  /// Calc shortest distance from p_in to def.line (or extention of finite line) at proj_pt
  double FindDistance(const Point& p_in) const;

  /// Calculate shortest distance between def.line and line_in.
  inline double FindDistance(const Line& line_in) const;

  /// Find shortest distance to line, and give back the two points that are nearest eachother.
  double FindDistance(const Line& line_in,
                      NRLib::Point & this_line_point,
                      NRLib::Point & in_line_point) const;

  /// Calc projection_point on def.line from p_in (or extention of finite line when proj.pt from p_in is outside)
  Point FindProjection(const Point& p_in) const;

  double FindParameter(const Point& p_in) const;

  /// Returns true if p_in is on the def.line, i.e. shortest distance to line (or extention) is less than min_in/0.0
  bool IsPointOnLine(const Point& p_in, double min_in = min_) const;

  /// Returns true when line_in and def.line are parallel (closer than min_in)
  bool IsParallel(const Line& line_in, double min_in = 0.0) const;

  /// True if def.line and line_in are overlapping (less than min_in):
  bool IsOverLapping(const Line& line_in, double min_in = min_) const;
private:
  Point p1_, p2_;            // points that define the line
  bool  end1_, end2_;        // true when p1_, p2_ are endpoints
  static const double min_;  // used instead of zero for some comparision-functions
};


// ================= INLINE FUNCTION DEFINITIONS =====================

double Line::FindDistance(const Line& line_in) const
{
  NRLib::Point p1, p2;
  return FindDistance(line_in, p1, p2);
}


} // end of namespace NRLib

#endif
